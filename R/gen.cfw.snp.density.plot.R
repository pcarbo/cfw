# Script to plot distribution of GBS SNPs across chromosomes 1-19.
library(lattice)

# LOAD GBS SNP DATA
# -----------------
# Load the SNP data genotyped in the CFW mouse cohort using
# genotyping-by-sequencing (GBS).
cat("Loading GBS SNP data.\n")
map <- read.csv("map.csv",stringsAsFactors = FALSE,quote = "",
                comment.char = "#")
map <- transform(map,
                 chr = factor(chr),
                 ref = factor(ref),
                 alt = factor(alt))

# CREATE SNP-DENSITY PLOT
# -----------------------
chromosomes <- list(10:1,19:11)
density     <- vector("list",19)
for (k in 1:2) {
  trellis.device(height = 6,width = 6)
  trellis.par.set(par.ylab.text = list(cex = 0.75),
                  axis.text = list(cex = 0.6),
                  axis.line = list(col = "white"))
  bins <- seq(0,200,1)
  ypos <- (-0.05)
  for (i in chromosomes[[k]]) {
    r <- table(cut(subset(map,chr == i)$pos/1e6,bins))
    density[[i]] <- r
    r[r>100] <- 100
    print(barchart(r,horizontal = FALSE,box.width = 1,
                   col = "dodgerblue",border = "dodgerblue",
                   scales = list(x = list(limits = c(0,201),at = seq(0,200,50),
                                          labels = paste(seq(0,200,50),"Mb")),
                                 y = list(limits = c(0,100),at = c(0,50,100),
                                           labels = c(0,50,"100+"))),
                   xlab = "",ylab = paste("chr",i)),
          position = c(0,ypos,1,ypos + 0.2),
          more = TRUE)
    ypos <- ypos + 0.09
  }
}

