# Script to compare the GBS SNP data against the "gold standard"
# Wellcome Trust SNP data.
library(lattice)

# LOAD GBS SNP DATA
# -----------------
# Load the SNP data genotyped in the CFW mouse cohort using
# genotyping-by-sequencing (GBS).
cat("Loading GBS SNP data.\n")
map <- read.csv("../data/map.csv",stringsAsFactors = FALSE,quote = "",
                comment.char = "#")
map <- transform(map,
                 chr = factor(chr),
                 ref = factor(ref),
                 alt = factor(alt))

# LOAD WELLCOME-TRUST SNP DATA
# ----------------------------
# cat("Loading Wellcome-Trust SNP data.\n")
# map.ref        <- read.table("WT.snps.pos.txt.gz",sep = " ",quote = "")
# names(map.ref) <- c("chr","pos")
# map.ref        <- transform(map.ref,chr = factor(chr))

# CREATE SCATTERPLOT
# ------------------
# Count the number of markers within each 1-Mb bin in both the
# Wellcome-Trust and GBS panels.
# x <- NULL
# y <- NULL
# for (i in 1:19) {
#   r    <- subset(map.ref,chr == i)
#   bins <- seq(min(r$pos),max(r$pos),1e6)
#   x    <- c(x,table(cut(r$pos,bins)))
#   y    <- c(y,table(cut(subset(map,chr == i)$pos,bins)))
# }

# # Create the scatterplot using the lattice library.
# trellis.device(height = 4,width = 4)
# print(xyplot(y ~ x,data.frame(x = x,y = y),pch = 20,col = "dodgerblue",
#              cex = 0.65,xlab = "Wellcome-Trust SNP count",
#              ylab = "GBS SNP count",
#              scales = list(x = list(limits = c(-500,2e4)),
#                            y = list(limits = c(-10,300)))))
# rm(i,r)

# CREATE SNP DENSITY HISTOGRAMS
# -----------------------------
# trellis.device(height = 5,width = 4)
# trellis.par.set(par.main.text = list(cex = 1,font = 1))
# print(barchart(table(cut(log10(x+1),seq(0,5,5/30))),horizontal = FALSE,
#                col = "black",border = "black",xlab = "log10 SNP count",
#                ylab = "num. 1-Mb windows",main = "Wellcome-Trust SNP panel",
#                scales = list(x = list(at = c(1,30)))),
#       split = c(1,1,1,2),
#       more = TRUE)
# print(barchart(table(cut(log10(y+1),seq(0,3,0.1))),horizontal = FALSE,
#                col = "darkorange",border = "darkorange",
#                xlab = "log10 SNP count",ylab = "num. 1-Mb windows",
#                main = "GBS SNP panel",
#                scales = list(x = list(at = c(1,30)))),
#       split = c(1,2,1,2))

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
rm(bins,ypos,i,r,chromosomes)
