library(qqman)
library(magrittr)
library(dplyr)
library(wesanderson)

# ---- SETUP ----

setwd("C:/Users/tomdr/Google Drive/School/Aber/Masters/Dissertation/")

plink.dir <- "Plink/"

either.side <- 1000000 # bp either side of snp

cols <- rep(wes_palette("Darjeeling1"), length.out = 26)

processGWASoutput <- function(gwas.name, assoc.prefix, qq, bonf) {
  # Read list of .qassoc file in Plink/GWAS directory
  dir <- paste0(plink.dir, gwas.name, "/")
  files <- list.files(path=dir, pattern="\\.qassoc$", full.names=TRUE, recursive=FALSE)#
  
  # Set up snp list table
  all.snp.data <- list()
  
  # Set up pheno names list for named snp data 
  pheno.names <- list()
  
  # Set up output dir
  out.dir <- paste0(dir, "out")
  dir.create(out.dir)
  
  # For each .qassoc file
  for (file in files) {
    # ---- READING DATA ----
    
    # Cut off "Plink/" and ".qassoc" sections of filename to get pheno name
    pheno.name <- substr(file, nchar(dir)+nchar(assoc.prefix)+1, nchar(file)-7)
    
    # Read qassoc data
    qassoc.file <- read.table(file, header=TRUE)
    
    # Get gwas results in correct format, skipping nas
    qassoc.df <- na.omit(data.frame(P = qassoc.file$P, SNP = qassoc.file$SNP, CHR = qassoc.file$CHR, BP = qassoc.file$BP))
    
    # Read adjusted data
    qassoc.adj.file <- read.table(paste0(file, ".adjusted"), header=TRUE)
    
    # Get gwas results in correct format, skipping nas
    qassoc.adj.df <- data.frame(BONFP = qassoc.adj.file$BONF, SNP = qassoc.adj.file$SNP)
    
    # Get list of bonf-adjusted significant SNPs
    #snps <- qassoc.adj.df[qassoc.adj.df$BONFP < 0.05, ]$SNP
    
    # Get full SNP data from qassoc file
    #snp.data <- qassoc.df[qassoc.df$SNP %in% snps, ]
    
    # Get list of SNPs from using 0.00001 cutoff
    snp.data <- qassoc.df[qassoc.df$P < 0.00001, ]
    
    # Add lookup range for genome data viewer
    snp.data <- snp.data %>% mutate(LOOKUP = paste0("Chr",CHR, ":", BP-either.side, "-", BP+either.side))
    
    # ---- PLOTING AND SAVING----
    
    if (bonf) {
      # JPEG for manhattan plot with bonf-adjusted data
      jpeg(file=paste0(out.dir, "/", pheno.name, "_bonf.jpeg"), 
           width = 6, height = 6, units = 'in', res = 300)
      
      # Attach bonf-adjusted data to qassoc data
      merged.df <- merge(qassoc.df, qassoc.adj.df, by="SNP")
      
      manhattan(merged.df, p = "BONFP")
      
      # Add pheno name as title
      mtext(paste(pheno.name, "with bonferroni-adjusted P values"), side = 3, line = -3, outer = TRUE)
      
      dev.off()
    }
    
    if (qq) {
      # JPEG for qq file
      jpeg(file=paste0(out.dir, "/", pheno.name, "_qq.jpeg"), 
           width = 6, height = 6, units = 'in', res = 300)
      
      
      qq(qassoc.df$P)
      
      # Add pheno name as title
      mtext(paste("QQ plot of P values for", pheno.name), side = 3, line = -3, outer = TRUE)
      
      dev.off()
    }
    
    # JPEG for manhattan plot
    jpeg(file=paste0(out.dir, "/", pheno.name, ".jpeg"), 
         width = 6, height = 6, units = 'in', res = 300)
    
    # Manhattan
    manhattan(qassoc.df, genomewideline = FALSE, col = cols, annotatePval = 0.00001)
    
    # Add pheno name as title
    mtext(paste0("P-Values by position in genome for SNPs for ", pheno.name), side = 3, line = -3, outer = TRUE)
    
    # Reset for next JPEG
    dev.off()
    
    # Add snp df to list
    all.snp.data <- c(all.snp.data, list(snp.data))
    pheno.names <- c(pheno.names, pheno.name)

    # Save individual snp data file
    #write.table(snp.data, paste0(out.dir, "/", pheno.name, "_snps.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ") 
  }
  
  # Merge snp df list into single df
  names(all.snp.data) <- pheno.names
  merged.snp.data <- bind_rows(all.snp.data, .id = "LIPID")
  
  # Save file
  write.table(merged.snp.data, paste0(out.dir, "/significant_snps.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")
}

processGWASoutput("gwas2", "g2assoc.", qq = FALSE, bonf = FALSE)



