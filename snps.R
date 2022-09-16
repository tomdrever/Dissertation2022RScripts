# Read significant plink file
# Group by "SNP CHR BP LOOKUP" using dplyr

library(magrittr)
library(dplyr)

setwd("C:/Users/tomdr/Google Drive/School/Aber/Masters/Dissertation/")
file <- "Plink/gwas2/out/significant_snps.txt"

snps <- read.table(file, header=TRUE)

snps.grouped <- snps %>% group_by(across(all_of(c("SNP", "LOOKUP")))) %>% summarise(LIPIDS = paste(LIPID, collapse=","), PS = paste(P, collapse=","))

write.table(snps.grouped, "Plink/gwas2/out/significant_snps_grouped.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep=" ")