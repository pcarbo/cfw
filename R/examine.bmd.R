# A small script to compare the distribution of bone-mineral density
# (BMD) against BMD data from other studies.
library(lattice)
source("read.data.R")
source("plotting.tools.R")

# Load the phenotype data from the CFW mice. So that I can compare
# these BMD measurements to the BMD data from the HMDP mouse study, I
# do not remove outliers or transform the data.
pheno <- read.pheno("pheno.csv")

# Load the phenotype data from the HMDP mice.
hmdp.pheno <- transform(read.csv("../data/farber2011.TableS1.csv",
                                 header = TRUE,stringsAsFactors = FALSE,
                                 comment.char = "#",check.names = FALSE),
                        sex       = factor(sex,c("M","F")),
                        totalbody = 1000 * totalbody,
                        femur     = 1000 * femur,
                        spine     = 1000 * spine)

# Set up the graphics device, and set up the location of each panel in
# the grid.
trellis.device(height = 4,width = 5)
trellis.par.set(list(fontsize = list(text = 9),
                     layout.widths = list(right.padding = -1),
                     layout.heights = list(top.padding = -1,
                                           bottom.padding = 0)))
panel.layout <- create.grid.layout(2,2,c("BMD","totalbody","femur","spine"))

# Show the distribution of BMD in the CFW mice.
panel.col <- "darkorange"
print(histogram(~BMD,pheno,col = panel.col,border = panel.col,nint = 24,
                xlab = "BMD",ylab = "% of samples",xlim = c(40,160),
                xaxp = c(40,160,6),ylim = c(0,20)),
      split = grid.split(panel.layout,"BMD"),more = TRUE)
                      
# Show the distribution of total body BMD in the HMDP mice.
panel.col <- "dodgerblue"
print(histogram(~totalbody,hmdp.pheno,col = panel.col,border = panel.col,
                nint = 24,xlab = "Total Body BMD",ylab = "% of samples"),
      split = grid.split(panel.layout,"totalbody"),more = TRUE)

# Show the distribution of femur BMD in the HMDP mice.
print(histogram(~femur,hmdp.pheno,col = panel.col,border = panel.col,
                nint = 18,xlab = "Femur BMD",ylab = "% of samples",
                xlim = c(40,160),xaxp = c(40,160,6),ylim = c(0,20)),
      split = grid.split(panel.layout,"femur"),more = TRUE)

# Show the distribution of lumbar spine BMD in the HMDP mice.
print(histogram(~spine,hmdp.pheno,col = panel.col,border = panel.col,
                nint = 24,xlab = "Lumbar Spine BMD",ylab = "% of samples"),
      split = grid.split(panel.layout,"spine"))
