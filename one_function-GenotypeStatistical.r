#01Statistical_analysis_results_Of_genotype
#2023.04.09
#Lipengle

genotype_analysis <- function(input_file) {

  packages <- c("data.table", "CMplot", "dplyr", "ggplot2")
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }

  
  
  library(data.table)
  library(CMplot)
  library(dplyr)
  library(ggplot2)
  # 读取数据
  map1 <- fread(input_file, header = TRUE)
  map2 <- map1[, 1:4]
  
  # 提取SNP信息
  mm <- map2 %>% dplyr::select(SNP = 1, Chromosome = 3, Position = 4)
  
  # 绘制SNP分布图
  CMplot(mm, plot.type = "d", bin.size = 1e6, col = c("darkgreen", "yellow", "red"),
         file = "tiff", dpi = 300, file.output = TRUE, verbose = TRUE)
  
  # 提取基因型数据
  alleles <- c(map1$alleles)
  
  # 创建组合列表
  combinations_list <- lapply(alleles, function(allele) {
    bases <- strsplit(allele, "/")[[1]]
    base1 <- bases[1]
    base2 <- bases[2]
    
    c(allele, paste0(base1, base1), paste0(base1, base2), paste0(base2, base1), paste0(base2, base2), "NN")
  })
  
  # 转换为数据框
  final_data <- as.data.frame(matrix(unlist(combinations_list), ncol = 6, byrow = TRUE))
  colnames(final_data) <- c("alleles", "AA", "aa", "Aa", "Aa1", "NN")
  
  # 合并基因型数据和基本信息
  data0 <- cbind(final_data, map1[, 5:ncol(map1)])
  
  # 计算结果矩阵
  result <- matrix(0, nrow = nrow(data0), ncol = length(data0[, 2:6]))
  
  for (i in 1:ncol(data0)) {
    row_data <- data0[i, -(1:6)]  
    
    for (j in 1:length(data0[, 2:6])) {
      count <- sum(data0[, 2:6][i, j] == row_data)
      result[i, j] <- count
    }
  }
  
  colnames(result) <- c("count_AA", "count_Aa", "count_Aa1", "count_aa", "count_NN")
  
  # 合并结果到数据框
  data1 <- cbind(final_data, result)
  
  # 计算求和列
  sum_values <- rowSums(data1[, 7:11])
  data1$sum <- sum_values
  
  # 计算频率
  data1$freHom1 <- data1$count_AA / data1$sum
  data1$freHom2 <- data1$count_aa / data1$sum
  data1$freHet <- (data1$count_Aa + data1$count_Aa1) / data1$sum
  data1$freMissing <- data1$count_NN / data1$sum
  data1$allele.freP <- (data1$count_AA + 0.5 * (data1$count_Aa + data1$count_Aa1)) / data1$sum
  data1$allele.freQ <- (data1$count_aa + 0.5 * (data1$count_Aa + data1$count_Aa1)) / data1$sum
  data1$MAF <- pmin(data1$allele.freP, data1$allele.freQ)
  data1$PIC <- 1 - (data1$freHom1^2 + data1$freHom2^2 + data1$freHet^2)  #PIC
   cat("Genotype data statistics successfully completed.\n")
  
  # Frequency distribution plot for MAF
  data1 <- na.omit(data1)
  freq_plot <- ggplot(data1, aes(x = MAF)) +
    geom_histogram(binwidth = 0.02, fill = "#1f78b4", color = "#08519c", na.rm = FALSE) +
    labs(title = "Frequency Distribution of Minor Allele Frequency (MAF)",
         x = "Minor Allele Frequency (MAF)", y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.background = element_rect(fill = "white"))
 
  ggsave("Frequency Distribution of Minor Allele Frequency (MAF).png", plot = freq_plot, width = 8, height = 6, dpi = 300, bg = "white")
  cat("Frequency distribution plot for MAF plotted successfully.\n")
  
  # 写入结果到文件
  output_file <- paste0("genotype_analysis_results_", format(Sys.Date(), "%Y%m%d"), ".csv")
  write.csv(data1, output_file, row.names = TRUE)
  
}

devtools::install("E:/data_analysis/SGGS1")
library(genotypeAnalysis)

#' Perform Genotype Analysis
#'
#' This function analyzes genotype data and generates statistics.
#'
#' @param input_file Path to the input genotype data file.
#' @return A data frame with analysis results.
#' @import data.table
#' @import dplyr
#' @import ggplot2
#' @export
genotype_analysis <- function(input_file) {
  # Your function code goes here
}

library(devtools)
document()
build()
install.packages("myGenotypePackage", repos = NULL, type = "source")



library(genotypeAnalysis)

genotype_analysis("data")







