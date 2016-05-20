# Script to plot distribution of MegaMUGA SNPs that are polymorphic in
# CFW mice across chromosomes 1-19.
library(lattice)

# LOAD MegaMUGA SNP DATA
# ----------------------
# Load the SNP data genotyped in the CFW mouse cohort using MegaMUGA
# array.
cat("Loading MegaMUGA SNP data.\n")
map <- read.csv("../data/megamuga.txt.gz",stringsAsFactors = FALSE,
                quote = "",sep = " ",header = TRUE)
names(map) <- c("id","chr","pos")
map <- transform(map,chr = factor(chr,paste0("chr",c(1:19,"X","Y"))))
map <- transform(map,chr = as.integer(chr))
map <- subset(map,chr < 20)
map <- map[order(map$chr,map$pos),]

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
    print(barchart(r,horizontal = FALSE,box.width = 1,col = "red",
                   border = "red",
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

