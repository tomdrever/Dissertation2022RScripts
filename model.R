# ---- READ AND FILTER DATA ----

cols <- c("character", "factor", "factor", "factor", "factor",
          "factor", "factor", "numeric", "numeric", "factor")

setwd("C:/Users/tomdr/Google Drive/School/Aber/Masters/Dissertation/workspace/")

covfile <- read.table("data/Covfile.txt",
                      header = TRUE, colClasses = cols, skip = 1)

lipfile <- read.table("data/Lipid raw data.txt", header = TRUE)

#Fitting linear regression and saving residual data
#(from https://www.stat.purdue.edu/bigtap/online/docs/simple-association-test.html)

# Create matrices of data from phenotype and covariate files (without X.trait)

pheno <- as.matrix(lipfile[, 2:19])

# lipids C8 (3), C13 (6), C20.3n.3 (14), C20.4n.3 (15) have -9 values. We can ignore sum.
exclude <- c("C8", "C13", "C20.3n.3", "C20.4n.3", "Sum") 
pheno <- pheno[, -which(colnames(pheno) %in% exclude)]

cov <- data.matrix(covfile[, 2:10]) # (using data.matrix to convert chr->num AND fit to matrix in one step)

# ---- FIT MODELS ----

n <- dim(pheno)[1] #sample size
p <- dim(pheno)[2] #number of traits

fit <- list() # Create an empty list for models of each phenotype
pheno_residuals <- matrix(0, n, p) # Create an empty matrix for the residuals

for (i in 1:p) {
  fit[[i]] <- glm(pheno[, i] ~ cov, family = Gamma()) #for each trait, fit regression of phenotype against covariates
  message(paste0("Fit for ", colnames(pheno)[[i]]))
  message(paste0("AIC:", fit[[i]]$aic))
  pheno_residuals[, i] <- resid(fit[[i]]) #obtain residual phenotype data
}

# ---- ADD FID COLUMN AND

# add trait names back
pheno_out <- cbind(lipfile[, 1], pheno_residuals)

# create vector of -9s the length of the phenotype file, attach as first column
fids <- rep(-9, dim(pheno_out)[1])
pheno_out <- cbind(fids, pheno_out)

# add colnames
colnames(pheno_out) <- c("FID", "IID", colnames(pheno))

# output file with header row
write.table(pheno_out, "data/phenos_excluded_residuals_gamma.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ") 








