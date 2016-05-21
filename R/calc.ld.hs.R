# Script for estimating LD decay in the heterogeneous stock (HS)
# mice. Relevant publications for the HS data are Valdar et al, Nature
# Genetics, 2006; Shifman et al, PLoS Biology, 2006. These
# calculations use the genotype data for 1,940 HS mice downloaded from
# http://mus.well.ox.ac.uk/mouse/HS. Base-pair positions used in these
# calculations are from the "sex-averaged" genetic map downloaded from
# the same webpage. Positions for 9,416 SNPs (out of 11745 SNPs total)
# on autosomal chromosomes are provided in this genetic map. The
# base-pair positions are based on NCBI release 34 of the mouse genome
# assembly. The --r2 command in PLINK v1.9beta3.33 is used to compute
# the r^2 measure of LD (see Pritchard and Przeworski, 2001) based on
# maximum-likelihood estimates of the haplotype frequences.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
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
map <- read.table("chrall_sexaveraged.bim",stringsAsFactors = FALSE)
names(map) <- c("chr","id","dist","pos","A1","A2")
cat("Loaded SNP data for",nrow(map),"markers.\n")

# LOAD SNP ALLELE FREQUENCIES
# ---------------------------
cat("Loading SNP allele frequencies")
out <- NULL
system("mkdir -p out_plink")
for (i in 1:19) {
  geno.file <- sprintf("chr%d",i)
  system(sprintf(paste("plink2 --bfile %s --freq --make-founders",
                       "--out out_plink/hs"),geno.file))
  out <- rbind(out,read.table("out_plink/hs.frq",header = TRUE,
                              stringsAsFactors = FALSE))
}
rownames(out) <- out$SNP
map <- cbind(map,data.frame(maf = out[map$id,"MAF"]))
rm(out)

# COMPUTE LD USING PLINK
# ----------------------
cat("Computing LD using PLINK: ")
ld.matrix <- vector("list",19)
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep = "")
for (i in 1:19) {
  caterase(sprintf("chromosome %d",i))
  geno.file <- sprintf("chr%d",i)
  out.file  <- sprintf("out_plink/chr%d",i)
 system(sprintf("plink2 --bfile %s --r2 square --make-founders --out %s",
                geno.file,out.file),ignore.stdout = TRUE)
  r           <- fread(paste0(out.file,".ld"))
  class(r)    <- "data.frame"
  r           <- as.matrix(r)
  rownames(r) <- subset(map,chr == i)$id
  colnames(r) <- subset(map,chr == i)$id
  ld.matrix[[i]] <- r
}
cat("\n")
rm(geno.file,out.file,i,r)

# SAMPLE LD WITHIN BASE-PAIR INTERVALS
# ------------------------------------
# Retain only the SNPs with base-pair positions.
map <- subset(map,!is.na(pos))
source("sample.ld.R")

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(list = c("bins","ld","map"),
     file = "LDdecay.HS.RData")
