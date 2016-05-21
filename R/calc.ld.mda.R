# This is Peter's script for estimating LD decay in the panel of
# inbred strains genotyped using the Mouse Diversity Array (MDA).
#
# These calculations use the genotype data at for 30 inbred lab
# strains. We use only the 297,329 SNPs on chromosomes 1-19 with no
# missing genotypes that are polymorphic in the 30 inbred
# strains. These data were downloaded from the Mouse Diversity
# Genotyping Array web resource at the Jackson Labs Center for Genome
# Dynamics. All genomic positions are based on Mouse Genome Assembly
# 37 from the NCBI database (mm9).
#
# We use the --r2 command in PLINK v1.9beta3.33 to compute the r^2
# measure of LD (see Pritchard and Przeworski, 2001) based on
# maximum-likelihood estimates of the haplotype frequences.
#
library(data.table)

# SCRIPT PARAMETERS
# -----------------
# Which chromosome to analyze.
# which.chr <- 1

# PLINK file containing the genotype data.
geno.file <- "/tmp/pcarbo/data/mda/mda"

# These are the bins used for the LD calculations.
bins <- seq(0,5e6,1e5)

set.seed(1)

# LOAD SNP INFO
# -------------
cat("Loading SNP data.\n")
map        <- read.table("/tmp/pcarbo/data/mda/mda.bim",
                         stringsAsFactors = FALSE)
names(map) <- c("chr","id","dist","pos","A1","A2")
map        <- subset(map,chr == which.chr)
cat("Loaded SNP data for",nrow(map),"markers.\n")

# LOAD SNP ALLELE FREQUENCIES
# ---------------------------
cat("Loading SNP allele frequencies")
system(sprintf(paste("~/bin/plink2 --bfile %s --freq --chr %d",
                     "--out /tmp/pcarbo/plink/mda"),
               geno.file,which.chr))
out <- read.table("/tmp/pcarbo/plink/mda.frq",header = TRUE,
                  stringsAsFactors = FALSE)
map <- cbind(map,data.frame(maf = out$MAF))
rm(out)

# Sample LD estimates within each bin from this many randomly selected
# SNP pairs.
ns <- round(nrow(map)/10)

# COMPUTE LD USING PLINK
# ----------------------
ld.matrix <- vector("list",19)
cat("Computing LD using PLINK:\n")
out.file <- sprintf("/tmp/pcarbo/plink/chr%d",which.chr)
system(sprintf("~/bin/plink2 --bfile %s --chr %d --r2 square --out %s",
               geno.file,which.chr,out.file),ignore.stdout = FALSE)
out <- fread(paste0(out.file,".ld"),verbose = FALSE,showProgress = TRUE,
             sep = "\t",colClasses = rep("numeric",sum(map$chr == which.chr)),
             header = FALSE)
class(out)    <- "data.frame"
out           <- as.matrix(out)
rownames(out) <- map$id
colnames(out) <- map$id
ld.matrix[[which.chr]] <- out
cat("\n")
rm(out.file,out)

# SAMPLE LD WITHIN BASE-PAIR INTERVALS
# ------------------------------------
source("sample.ld.R")

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(list = c("bins","ld","map"),
     file = paste0("LDdecay.MDA.chr",which.chr,".RData"))
