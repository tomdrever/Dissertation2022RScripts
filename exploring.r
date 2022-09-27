library(dplyr)
library(magrittr)
library(ggplot2)
library(gridExtra)

# ---- READ DATA ----

cols <- c("character", "factor", "factor", "factor", "factor",
    "factor", "factor", "numeric", "numeric", "factor")

setwd("C:/Users/tomdr/Google Drive/School/Aber/Masters/Dissertation/workspace/")

covfile <- read.table("data/Covfile.txt", header = TRUE, colClasses = cols, skip = 1)

lipfile <- read.table("data/Lipid raw data.txt", header = TRUE)


# ---- INVESTIGATE COVARIATES ---- 
# Which are important/how distributed/how usable in model (glm?)

covgraphs <- list()

# Plot histogram of each
for (cov in list("Myomax", "Breed", "Year", "BirthFlock", "LambSex", "Rtype", "KillBatch")) {
  covgraphs[[cov]] <- ggplot(covfile, aes(x=.data[[cov]])) + geom_histogram(stat="count")
}

out <- arrangeGrob(grobs = covgraphs, ncol = 4, respect = TRUE)
ggsave("../img/cov_hists.pdf", out)

# ---- INVESTIGATE PHENOTYPES ---- 
# Histograms for phenotypes / fatty acids too

lipgraphs <- list()

for (lip in colnames(lipfile)[2:19]) {
  lipgraphs[[lip]] <- ggplot(lipfile, aes(x=.data[[lip]])) + geom_histogram()
}

out1 <- arrangeGrob(grobs = lipgraphs[1:9], ncol = 3, respect = TRUE)
out2 <- arrangeGrob(grobs = lipgraphs[10:18], ncol = 3, respect = TRUE)
ggsave("../img/lipid_hists_1.pdf", out1)
ggsave("../img/lipid_hists_2.pdf", out2)

# ---- CHECKING MISSING VALUES ----

famfile <- read.table("../data/filtered50KGWAS.fam")

# Check unique traits in each list of ids
pheno_only <- setdiff(lipfile$X.Trait., famfile$V2)
length(pheno_only)
geno_only <- setdiff(famfile$V2, lipfile$X.Trait.)
length(geno_only)

# ---- INVESTIGATE PHENOTYPE CORRELATION ----

# Check correlations of all relevant pheno variables
rel_lipfile <- lipfile[-c(1, 3, 6, 14, 15, 19)]

# Big graph so try pdf
pdf("../img/pheno_correlation.pdf")
# And split in 2
pairs(rel_lipfile[1:6], cex=0.5)
pairs(rel_lipfile[7:13], cex=0.5)
dev.off()





