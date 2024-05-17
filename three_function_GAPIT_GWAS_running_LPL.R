
run_GAPIT <- function(pheno_file, geno_file, snp_file) {
  source("http://zzlab.net/GAPIT/gapit_functions.txt") 
  # Read phenotype, genotype, and SNP information
  Phe=read.table(file="./data/pheno.txt",head=T)
  GD=read.table(file="./data/geno.txt", sep="\t",head=T)
  GM=read.table(file="./data/SNPinfor.txt", sep="\t",head=T)
  
  # Extract column names of phenotypes
  col_names <- names(Phe)[2:length(names(Phe))]
  setwd("~/GWAS.analysis.results")
  
  
  # Create directories for each phenotype, set working directory, and run GAPIT
  for (col in col_names) {
    dir_name <- paste("Phe_", col, sep = "")  # Create directory name
    dir.create(dir_name, showWarnings = FALSE)  # Create directory
    
    original_dir <- getwd()  # Get original working directory
    
    # Set working directory and check if data exists in the selected column
    setwd(dir_name)
    if (length(Phe[[col]]) == 0) {
      cat("No data found in column", col, "\n")
      setwd(original_dir)  # Set working directory back to original
      next  
    }
    
    # Run GAPIT
    myGAPIT_SUPER <- GAPIT(
      Y = Phe[, c(1, which(col_names == col) + 1)],  # Set phenotype data
      GD = GD,
      GM = GM,
      PCA.total = 3,
      model = c("Blink")
    )
    
    setwd(original_dir)  # Set working directory back to original
  }
}

# Example usage:
 run_GAPIT("./data/pheno.txt", "./data/geno.txt", "./data/SNPinfor.txt")
















