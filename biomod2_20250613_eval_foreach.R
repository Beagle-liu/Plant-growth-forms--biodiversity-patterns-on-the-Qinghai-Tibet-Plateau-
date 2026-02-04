###############################################################################################
###############################################################################################
####	
####	Species distribution modeling of plant with Biomod2
####
####	Programmed in R 2025
####
###############################################################################################
###############################################################################################
rm(list = ls())
library(biomod2)
library(raster)
library(lattice)
library(rasterVis)
library(gridExtra)
library(reshape2)
library(terra)
library(latticeExtra)
library(snowfall)
library(tidyterra)
library(mda)
library(earth)
library(maxnet)
library(randomForest)
library(gam)
library(gbm)
library(mgcv)
library(xgboost)
library(usdm)
library(TeachingDemos)
library(nlme)
library(rpart)
library(class)
library(splines)
library(plotmo)
library(plotrix)
library(foreach)
library(doParallel)



setwd("F:/LF/SDM/QTPOUT_Tree2/")
QTP_Plant_Importance <- data.frame()
QTP_Plant_Eval <- data.frame()

# Current scenario
ExplVar <- stack(list.files("F:/LF/SDM/QTP/", pattern = ".tif$", full.names = TRUE))
names(ExplVar) <- c("bio1","bio10","bio11","bio12","bio13","bio14","bio15",
                    "bio16","bio17","bio18","bio19","bio2","bio3","bio4","bio5",
                    "bio6","bio7","bio8","bio9","Clay","Elevation","Elevation_range"
                    ,"Elevation_std","gdd","Soil_TEB","Total_N","Tri")

# Load species list
species_list <- read.csv("F:/LF/SDM/QTP5/QTP_VP_GF5_Tree2.csv")
Species <- unique(species_list$species)


# Load species variables
varLIST <- read.csv("F:/LF/SDM/QTP5/selected_env_variables5_Tree2.csv")


## build species modelling wrapper ----
biomod2_wrapper <- function(spLIST){
  cat("\n> species : ", spLIST)
  
  sp_dat <- subset(species_list, species == spLIST)
  
  # Select the species-specific variables set
  
  
  # 检查变量是否存在
  explVarNames <- varLIST$selected_env[varLIST$species == spLIST]
  validVarNames <- explVarNames[explVarNames %in% names(ExplVar)]
  ExplVar_cropped <- ExplVar[[validVarNames]]
  
  # Format data for modeling
  myBiomodata <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sp_dat)), 
    resp.name = spLIST,
    expl.var = ExplVar_cropped, 
    resp.xy = sp_dat[, c("long", "lat")],
    PA.nb.rep = 1, 
    PA.nb.absences = 10000,
    PA.strategy = 'random',
    na.rm = TRUE)
  
  # Modeling options
  allModels <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                 , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
  
  # providing formated data
  sp_opt <- bm_ModelingOptions(
    data.type = 'binary',
    models = allModels,
    strategy = 'default')
  
  # Model fitting
  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodata, 
    models = c('GBM','GAM','RF', 'MAXENT.Phillips.2', 'CTA'),
    bm.options = sp_opt,
    CV.strategy = "kfold",
    CV.k = 5,
    CV.nb.rep = 5,
    metric.eval = c('TSS', 'ROC',"BOYCE"),
    var.import = 1,
    do.progress = FALSE,
    weights = NULL, # new name for Yweights
    scale.models = FALSE, # new name for rescal.all.models 
    CV.do.full.models = FALSE,
    modeling.id = spLIST
  )
  
  # Ensemble modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = 'all',
    em.by = 'all',
    em.algo = "EMca",
    metric.select = c("BOYCE"),
    metric.select.thresh = c(0.4),
    metric.eval = c('ROC',"TSS","BOYCE"),
    var.import = 1,
    do.progress = FALSE)
  
  # Variables importance
  # Variables importance
  variables_importance <- get_variables_importance(myBiomodEM)
  variables_importance <- as.data.frame(variables_importance)
  variables_importance <- variables_importance[, c("expl.var", "var.imp")]
  variables_importance$species <- paste(spLIST)
  QTP_Plant_Importance <- rbind(QTP_Plant_Importance, variables_importance)
  
  #Store evaluation metrics
  QTP_Plant_Eval[spLIST, 1] <- get_evaluations(myBiomodEM)[1, c('calibration')] # ROC
  QTP_Plant_Eval[spLIST, 2] <- get_evaluations(myBiomodEM)[2, c('sensitivity')]
  QTP_Plant_Eval[spLIST, 3] <- get_evaluations(myBiomodEM)[2, c('specificity')]
  QTP_Plant_Eval[spLIST, 4] <- get_evaluations(myBiomodEM)[2, c('calibration')] # TSS
  QTP_Plant_Eval[spLIST, 5] <- get_evaluations(myBiomodEM)[3, c('calibration')] # TSS
  
  
  
  #Project under current conditions
  myBiomodPresEF <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    new.env = ExplVar_cropped,
    models.chosen = "all",
    metric.binary = 'TSS',
    compress = TRUE,
    build.clamping.mask = FALSE,
    do.stack = FALSE,
    proj.name = "curEM.prj",
    output.format = ".tif"
  )
  
  
  QTP_Plant_Eval <- as.data.frame(QTP_Plant_Eval)
  QTP_Plant_Eval$Species <- rownames(QTP_Plant_Eval)
  
  write.table(variables_importance, paste0('QTP_Plant_Importance_', spLIST, '.csv'), sep = ";", row.names = T)
  write.table(QTP_Plant_Eval, paste0('QTP_Plant_Eval_', spLIST, '.csv'), sep = " ", row.names = F)
  
  # 显式清理内存
  rm(myBiomodata, myBiomodModelOut, myBiomodEM)
  gc()
}

## launch the spiecies modelling wrapper over species list ----
if (require(snowfall)) {
  sfInit(parallel = TRUE, cpus = 12)  # 根据实际情况调整CPU核心数
  sfExportAll()  # 导出所有变量到各节点
  # 显式加载所有需要的包到每个节点
  sfLibrary(biomod2)
  sfLibrary(raster)
  sfLibrary(terra)
  sf_out <- sfLapply(Species, biomod2_wrapper)
  sfStop()
} else { ## sequencial computation
  for (spLIST in Species){
    biomod2_wrapper(spLIST)
  }
  ## or with a lapply function in sequential model
  ## all_species_bm <- lapply(Species, biomod2_wrapper)
}