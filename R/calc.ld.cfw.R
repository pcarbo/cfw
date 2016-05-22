# Script for estimating LD decay in the CFW mice. These calculations
# use the genotype data for 1,150 male CFW mice at 92,734 SNPs on
# autosomal chromosomes. The base-pair positions are based on NCBI
# release 38 of the mouse genome assembly. The --r2 command in PLINK
# v1.9beta3.33 is used to compute the r^2 measure of LD (see Pritchard
# and Przeworski, 2001) based on maximum-likelihood estimates of the
# haplotype frequencies.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
# PLINK .bed/.bim/.fam files containing the genotype data.
geno.file <- "cfw"

# These are the bins used for the LD calculations.
bins <- seq(0,5e6,1e5)

# Sample LD estimates within each bin from this many randomly selected
# SNP pairs.
ns <- 1e4

set.seed(1)

# LOAD SNP INFO
# -------------
cat("Loading SNP data.\n")
map <- read.table("cfw.bim",stringsAsFactors = FALSE)
names(map) <- c("chr","id","dist","pos","A1","A2")
cat("Loaded SNP data for",nrow(map),"markers.\n")

# LOAD SNP ALLELE FREQUENCIES
# ---------------------------
cat("Loading SNP allele frequencies")
system("mkdir -p out_plink")
system(sprintf("plink2 --bfile %s --freq --out out_plink/cfw",geno.file))
out <- read.table("out_plink/cfw.frq",header = TRUE,
                  stringsAsFactors = FALSE)
map <- cbind(map,data.frame(maf = out$MAF))
rm(out)

# COMPUTE LD USING PLINK
# ----------------------
ld.matrix <- vector("list",19)
cat("Computing LD using PLINK: ")
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep = "")
for (i in 1:19) {
  caterase(sprintf("chromosome %d",i))
  out.file <- sprintf("out_plink/chr%d",i)
  system(sprintf("plink2 --bfile %s --chr %d --r2 square --out %s",
                 geno.file,i,out.file),ignore.stdout = TRUE)
  r <- fread(paste0(out.file,".ld"),verbose = FALSE,showProgress = FALSE)
  class(r)    <- "data.frame"
  r           <- as.matrix(r)
  rownames(r) <- subset(map,chr == i)$id
  colnames(r) <- subset(map,chr == i)$id
  ld.matrix[[i]] <- r
}
cat("\n")
rm(out.file,i,r)

# SAMPLE LD WITHIN BASE-PAIR INTERVALS
# ------------------------------------
source("sample.ld.R")

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(list = c("bins","ld","map"),
     file = "LDdecay.CFW.RData")
