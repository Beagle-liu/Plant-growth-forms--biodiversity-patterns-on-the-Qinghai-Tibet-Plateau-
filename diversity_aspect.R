# 加载必要的包
library(phyloraster)
library(terra)
library(ape)
library(phytools)

# 设置路径和参数
setwd("F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Epiphytes/B_Correct")
tif_dir <- "F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Epiphytes/B_Correct"  # 替换为栅格文件夹路径
tree_file <- "F:/Academic_Data/QTPE_VP_GF_distribution/QTP_tree_scenario.3.newick"  # 替换为系统发育树文件路径
# --- 步骤1: 准备物种分布栅格数据 ---
# 获取所有.tif文件路径
tif_files <- list.files(pattern = "\\.tif$", full.names = TRUE)

# 读取并堆叠栅格
sp_dist <- rast(tif_files)
names(sp_dist) <- tools::file_path_sans_ext(basename(tif_files))  # 用文件名作为层名

tree <- read.tree("F:/Academic_Data/QTPE_VP_GF_distribution/QTP_tree_scenario.3.newick")  # 替换为实际文件路径



tree_pruned <- keep.tip(tree, names(sp_dist))

# 启用并行计算（在计算前设置）
phyloraster::set_parallel(TRUE, ncores = 8)  # 使用8核
# 计算物种丰富度 (SR)
sr <- rast.sr(sp_dist)
plot(we, main = "Species Richness")

# 计算系统发育多样性 (PD)
pd <- rast.pd(x = sp_dist, tree = tree_pruned)
# 计算加权特有性 (WE)
we <- rast.we(sp_dist)

# 计算系统发育特有性 (PE)
pe <- rast.pe(x = sp_dist, tree = tree_pruned)

# --- 步骤5: 保存结果 ---
writeRaster(sr, "F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Diversity_aspect/Epiphytes_SR.tif", overwrite = TRUE)
writeRaster(we, "F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Diversity_aspect/Epiphytes_WE.tif", overwrite = TRUE)
writeRaster(pd, "F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Diversity_aspect/Epiphytes_PD.tif", overwrite = TRUE)
writeRaster(pe, "F:/Academic_Data/QTPE_VP_GF_distribution/QTP/Current/GF/Diversity_aspect/Epiphytes_PE.tif", overwrite = TRUE)