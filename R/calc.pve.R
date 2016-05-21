# Estimate the proportion of variance in the phenotype, after removing
# linear effects of covariates, that is explained by the available
# genetic variants.
library(data.table)
source("misc.R")
source("read.data.R")
source("data.manip.R")
source("qtl.analyses.R")
source("polygenic.R")

# SCRIPT PARAMETERS
# -----------------
which.analysis <- "tibia"
incl.markers   <- NULL

# Get the phenotype and covariates used in the QTL mapping, and the
# pathname of the file for saving the results.
analysis   <- analyses[[which.analysis]]
phenotype  <- analysis$pheno
covariates <- analysis$cov
outliers   <- analysis$outliers
h          <- seq(0.01,0.99,0.01)

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
if (!is.null(outliers))
  pheno <- remove.outliers(pheno,phenotype,covariates,outliers)

# Convert the "abnormal" BMD phenotype from a factor to a number.
pheno <- transform(pheno,abBMD = binfactor2num(abBMD))

# Only analyze samples (i.e. rows of the genotype and phenotype
# matrices) for which the phenotype and all the covariates are
# observed.
pheno <- pheno[which(none.missing.row(pheno[c(phenotype,covariates)])),]

# LOAD SNP DATA
# -------------
# Load the marker info and the mean alternative allele counts ("dosages").
cat("Loading SNP data.\n")
map     <- read.map("map.txt")
out     <- read.geno.dosage("geno.txt",nrow(map))
discard <- out$discard
X       <- out$geno
rm(out)

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
markers <- which(map$quality > 0.95 & f > 0.02)
map     <- map[markers,]
X       <- X[,markers]

# If requested, include SNPs as covariates in the linear regression of
# the phenotype.
if (length(incl.markers) > 0) {
  covariates <- c(covariates,incl.markers)
  markers    <- match(incl.markers,map$snp)
  d          <- data.frame(X[,markers])
  names(d)   <- incl.markers
  X          <- X[,-markers]
  map        <- map[-markers,]
  pheno      <- cbind(pheno,d)
}

# Get the phenotype data.
y <- pheno[[phenotype]]
n <- length(y)

# Get the covariate data.
Z <- pheno[covariates]
for (col in covariates)
  if (is.factor(Z[[col]]))
    Z[[col]] <- binfactor2num(Z[[col]])
Z <- as.matrix(cbind(data.frame(intercept = rep(1,n)),Z))

# Adjust the genotypes and phenotypes so that the linear effects of
# the covariates are removed. This is equivalent to integrating out
# the regression coefficients corresponding to the covariates with
# respect to an improper, uniform prior; see Chipman, George and
# McCulloch, "The Practical Implementation of Bayesian Model
# Selection," 2001. The equivalent expressions in MATLAB are  
#
#   y = y - Z*((Z'*Z)\(Z'*y))
#   X = X - Z*((Z'*Z)\(Z'*X))  
#
# Note that this should give the same result as centering the
# columns of X and subtracting the mean from y when we have only
# one covariate, the intercept.
y <- y - c(Z %*% solve(crossprod(Z),t(y %*% Z)))
X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)

# COMPUTE POLYGENIC MODEL ESTIMATES
# ---------------------------------
# Give summary of the analysis.
cat("Evaluating polygenic model for",phenotype,"in",n,"mice with",
    nrow(map),"SNPs,\n")
if (!is.null(covariates)) {
  cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
} else {
  cat("with no covariates included.\n")
}
out   <- polygenic.model(X,y,h)
logw  <- out$logw
sigma <- out$sigma
rm(out)

# Compute the normalized importance weights.
w <- normalizelogweights(logw)
