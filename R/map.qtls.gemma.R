# Map QTLs for phenotypes measured in CFW outbred mice using the linear
# mixed model (LMM) analysis implemented in GEMMA.
library(qtl)
library(data.table)
source("misc.R")
source("gemma.R")
source("read.data.R")
source("data.manip.R")
source("qtl.analyses.R")

# SCRIPT PARAMETERS
# -----------------
which.analysis <- "EDL2"
chromosomes    <- NULL
incl.markers   <- NULL
gemmadir       <- "gemma_out"
resultdir      <- "/Users/pcarbo"
gemma.exe      <- "~/shyamg/bin/gemma"
num.perms      <- 0
seed           <- 1

# Get the phenotype and covariates used in the QTL mapping.
analysis   <- analyses[[which.analysis]]
phenotype  <- analysis$pheno
covariates <- analysis$cov
outliers   <- analysis$outliers

# LOAD PHENOTYPE DATA
# -------------------
# Load the phenotype data, and discard outlying phenotype values. I
# create binary covariates from some of the categorical phenotypes.
cat("Loading phenotype data.\n")
pheno <- read.pheno("pheno.csv")
pheno <- prepare.pheno(pheno)
pheno <- cbind(pheno,
               binary.from.categorical(pheno$FCbox,paste0("FCbox",1:4)),
               binary.from.categorical(pheno$PPIbox,paste0("PPIbox",1:5)),
               binary.from.categorical(pheno$methcage,
                                       paste0("methcage",1:12)),
               binary.from.categorical(pheno$round,paste0("SW",1:25)))
stop()
if (!is.null(outliers))
  pheno <- remove.outliers(pheno,phenotype,covariates,outliers)

# Only analyze samples (i.e. rows of the genotype and phenotype
# matrices) for which the phenotype and all the covariates are
# observed.
pheno <- pheno[which(none.missing.row(pheno[c(phenotype,covariates)])),]

# LOAD GENOTYPE DATA
# ------------------
# Load the "mean genotypes", or the mean alternative allele counts.
cat("Loading genotype data.\n")
load("../data/geno.dosage.RData")

# Discard genotype samples from mislabeled flowcell samples.
X <- X[which(discard == "no"),]

# Align the phenotypes and genotypes
ids   <- intersect(pheno$id,rownames(X))
pheno <- pheno[match(ids,pheno$id),]
X     <- X[match(ids,rownames(X)),]

# Discard SNPs with low "imputation quality" assessed by inspecting
# the genotype probabilities. Retain SNPs for which: (1) at least 95%
# of the samples have a maximum probability genotype greater than than
# 0.5; (2) the minor allele frequency is greater than 2%.
f       <- apply(X,2,compute.maf)
markers <- which(quality > 0.95 & f > 0.02)
map     <- map[markers,]
X       <- X[,markers]

# If requested, include SNPs as covariates in the linear regression of
# the phenotype.
if (length(incl.markers) > 0) {
  covariates <- c(covariates,incl.markers)
  markers    <- match(incl.markers,map$snp)
  d          <- data.frame(X[,markers])
  names(d)   <- incl.markers
  pheno      <- cbind(pheno,d)
}

# Initialize the random number generator.
set.seed(seed)

if (num.perms == 0) {

  # MAP QTLs 
  # --------
  # Calculate p-values using GEMMA.
  gwscan.gemma <- run.gemma(phenotype,covariates,pheno,X,map,
                            gemmadir,gemma.exe,chromosomes)

  # Save results to file.
  cat("Saving results to file.\n")
  save(list = c("analysis","gwscan.gemma"),
       file = paste0(resultdir,"/",analysis$file))
} else {

  # ESTIMATE NULL DISTRIBUTION USING PERMUTATIONS
  # ---------------------------------------------
  # Set the local directory to the location of the GEMMA files.
  setwd(gemmadir)
  
  # Create text files containing the mean genotypes and map information for
  # all markers in the format used by GEMMA.
  cat("Writing SNP and genotype data to separate files.\n")
  write.gemma.geno(paste0(gemmadir,"/geno.txt"),X,map)
  write.gemma.map(paste0(gemmadir,"/map.txt"),map)

  # Write out the kinship matrix to file.
  cat("Writing identity kinship matrix to file.\n");
  write.table(diag(nrow(pheno)),paste0(gemmadir,"/kinship.txt"),
              sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
  
  # Initialize the table containing the minimum p-values.
  perms.gemma           <- matrix(NA,num.perms,1)
  colnames(perms.gemma) <- which.analysis
  
  # Repeat for each permutation.
  for (i in 1:num.perms) {

    # Permute the rows of the phenotype table.
    cat("Permutation #",i,'\n',sep="")
    rows <- sample(nrow(pheno))
      
    # Write the phenotype and covariate data to separate text files.
    cat(" * Writing phenotype and covariate data to file.\n")
    write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno[rows,])
    write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),
                           covariates,pheno[rows,])

    # Now we are ready to calculate the p-values for all markers using
    # GEMMA.
    cat(" * Computing p-values for all markers using GEMMA.\n")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)

    # Load the results of the GEMMA association analysis, and get the
    # minimum p-value.
    gwscan <- read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
    perms.gemma[i] <- max(gwscan$log10p)
  }
  
  # Save results to file.
  cat("Saving results to file.\n")
  class(perms.gemma) <- c("scanoneperm","matrix")
  save(list = c("seed","analysis","perms.gemma"),
       file = paste0(resultdir,"/perms.",analysis$file))
}
