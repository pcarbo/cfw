# This R script is used in the calc.ld.*.R scripts to sample LD
# estimates within given base-pair intervals, using "allele frequency
# matching". See the supplementary methods of the Nature Genetics for
# details.
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep = "")

# Initialize storage for the results.
p            <- nrow(map)
nk           <- length(bins) - 1
ld           <- matrix(0,ns,nk)
colnames(ld) <- bins[1:nk]/1e6
markers      <- which(map$maf > 0.2)
    
# Repeat for each base-pair interval.
for (k in 1:nk) {
  a <- bins[k]
  b <- bins[k+1]
  cat(sprintf("Sampling LD for bp distances in %0.2f-%0.2f Mb: ",a/1e6,b/1e6))

  # Sample SNP pairs within the selected base-pair interval.
  iter <- 0
  while (iter < ns) {
    caterase(iter)

    # Get a random pair of SNPs.
    i <- sample(markers,1)
    j <- with(map,which(chr == chr[i] &
                        abs(maf - maf[i]) < 0.05 &
                        abs(pos - pos[i]) > a &
                        abs(pos - pos[i]) < b))
    if (length(j) > 0) {
      iter <- iter + 1
      if (length(j) > 1)
        j <- sample(j,1)
      ld[iter,k] <- with(map,ld.matrix[[chr[i]]][id[i],id[j]])
    }
  }
  cat("\n")
}
rm(p,nk,k,a,b,iter,i,j)
