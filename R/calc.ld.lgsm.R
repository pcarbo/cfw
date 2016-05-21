# Script for estimating LD decay in mice from the LG x SM advanced
# intercross line. These calculations use the genotype data at 4,524
# SNPs on chromosomes 1-19 for 687 F34 LG/J x SM/J crosses (353 males
# and 334 females). All SNP information and genomic positions are
# based on Mouse Genome Assembly 37 (mm9, July 2007). The --r2 command
# in PLINK v1.9beta3.33 is used to compute the r^2 measure of LD (see
# Pritchard and Przeworski, 2001) based on maximum-likelihood
# estimates of the haplotype frequences.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
# PLINK .bed/.bim/.fam files containing the genotype data.
geno.file <- "lgsm"

# These are the bins used for the LD calculations.
bins <- seq(0,5e6,1e5)

# Sample LD estimates within each bin from this many randomly selected
# SNP pairs.
ns <- 1e4

# Initialize the pseudorandom number generator.
set.seed(1)

# LOAD SNP INFO
# -------------
cat("Loading SNP data.\n")
map <- read.table("lgsm.bim",stringsAsFactors = FALSE)
names(map) <- c("chr","id","dist","pos","A1","A2")
cat("Loaded SNP data for",nrow(map),"markers.\n")

# LOAD SNP ALLELE FREQUENCIES
# ---------------------------
cat("Loading SNP allele frequencies")
system("mkdir -p out_plink")
system(sprintf("plink2 --bfile %s --freq --out out_plink/lgsm",geno.file))
out <- read.table("out_plink/lgsm.frq",header = TRUE,stringsAsFactors = FALSE)
map <- cbind(map,data.frame(maf = out$MAF))
rm(out)

# Keep only SNPs on chromosomes 1-19.
map <- subset(map,is.element(chr,1:19))

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
save(list = c("bins","ld.lgsm","map"),
     file = "LDdecay.LGSM.RData")
