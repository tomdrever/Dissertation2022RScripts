# ---- READ DATA ----

cols <- c("character", "factor", "factor", "factor", "factor",
          "factor", "factor", "numeric", "numeric", "factor")

setwd("C:/Users/tomdr/Google Drive/School/Aber/Masters/Dissertation/workspace/")

covfile <- read.table("data/Covfile.txt",
                      header = TRUE, colClasses = cols, skip = 1)

lipfile <- read.table("data/Lipid raw data.txt", header = TRUE)

# ---- FILTER DATA ----

# Create matrices of data from phenotype and covariate files (without X.trait)
pheno <- as.matrix(lipfile[, 2:19])

# Lipids C8 (3), C13 (6), C20.3n.3 (14), C20.4n.3 (15) have -9 values. We can ignore sum.
exclude <- c("C8", "C13", "C20.3n.3", "C20.4n.3", "Sum") 
pheno <- pheno[, -which(colnames(pheno) %in% exclude)]

# Using data.matrix to convert chr->num AND fit to matrix in one step
cov <- data.matrix(covfile[, 2:10])

# ---- FIT MODELS ----

n <- dim(pheno)[1] # Sample size
p <- dim(pheno)[2] # Number of traits

# Create an empty list for models of each phenotype
fit <- list()
# Create an empty matrix for the residuals
pheno_residuals <- matrix(0, n, p) 

for (i in 1:p) {
  # For each trait, fit model of phenotype against covariates
  fit[[i]] <- glm(pheno[, i] ~ cov, family = Gamma()) 
  
  message(paste0("Fit for ", colnames(pheno)[[i]]))
  message(paste0("AIC:", fit[[i]]$aic))
  
  # Store residual phenotype data
  pheno_residuals[, i] <- resid(fit[[i]]) 
}

# ---- ADD EMPTY ID COLUMNS ----
# This format is needed for PLINK

# Add trait names back
pheno_out <- cbind(lipfile[, 1], pheno_residuals)

# Create vector of -9s the length of the phenotype file, attach as first column
fids <- rep(-9, dim(pheno_out)[1])
pheno_out <- cbind(fids, pheno_out)

# Add colnames
colnames(pheno_out) <- c("FID", "IID", colnames(pheno))

# Output file with header row
write.table(pheno_out, "data/phenos_excluded_residuals_gamma.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ") 
