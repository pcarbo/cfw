# SUMMARY
# -------
# Some functions for creating plots to summarize results. Here is an
# overview of the functions defined in this file:
#
#   create.grid.layout(nrow,ncol,panels)
#   grid.split(panel.layout,panel)
#   effect.plot(pheno,geno,map,phenotype,marker)
#   odds.plot(pheno,geno,map,phenotype,marker)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Return a data structure detailing a grid layout for figure panels
# with the specified number of rows and columns. List elements "row"
# and "col" give the row and column number of each of the figure
# panels; each column is filled from left to right. The figure panels
# are named according to the input "panels".
create.grid.layout <- function (nrow, ncol, panels = NULL) {
  out <- list(nrow = nrow,
              ncol = ncol,
              row  = rep(1:nrow,times = ncol),
              col  = rep(1:ncol,each  = nrow))
  names(out$row) <- panels
  names(out$col) <- panels
  return(out)
}

# ----------------------------------------------------------------------
# Return the "split" argument for the "print.trellis" function that
# corresponds to the specified panel and grid layout.
grid.split <- function (panel.layout, panel)
  with(panel.layout,c(col[panel],row[panel],ncol,nrow))

# ----------------------------------------------------------------------
# Produce a box-percentile plot to visualize the relationship between
# phenotype and genotype.
effect.plot <- function (pheno, geno, map, phenotype, marker) {

  # Convert the genotype counts to genotypes.
  A            <- map$ref
  B            <- map$alt
  genotypes    <- c(paste0(A,A),paste0(A,B),paste0(B,B))
  levels(geno) <- genotypes
    
  # Create a single table containing the genotype and phenotype data.
  d        <- as.data.frame(cbind(pheno,geno))
  names(d) <- c("pheno","geno")
  
  # Output the figure.
  ylabels <- paste0(genotypes," (",table(geno),")")
  return(bwplot(formula("geno ~ pheno"),d,
                probs = c(0.01,0.05,0.125,0.25,0.375),
                panel = panel.bpplot,means = FALSE,nout = 0.003,
                scat1d.opts = list(lwd = 4,tfrac = 0.01,col = "slategray"),
                xlab = phenotype,ylab = "genotype",
                main = paste0(marker," (chr ",map$chr,", ",
                              round(map$pos/1e6,digits = 2)," Mb)"),
                par.settings = list(par.main.text = list(cex = 1,font = 1),
                                    layout.heights = list(axis.top = 0)),
                scales = list(x = list(tck = 0.5),
                              y = list(labels = ylabels,tck = 0))))
}

# ----------------------------------------------------------------------
# Produce a bar chart to visualize the relationship between the
# genotype and a binary trait or phenotype.
odds.plot <- function (pheno, geno, map, phenotype, marker) {

  # Convert the binary trait to a factor.
  pheno <- factor(pheno)
    
  # Convert the genotype counts to genotypes.
  A            <- map$ref
  B            <- map$alt
  genotypes    <- c(paste0(A,A),paste0(A,B),paste0(B,B))
  levels(geno) <- genotypes

  # Create a 3 x 2 contingency table of the genotype (rows) versus
  # the phenotype (columns).
  counts <- table(geno,factor(as.numeric(pheno) - 1,1:0))

  # Divide by the total count in each row so that the row entries add
  # up to 1.
  for (r in 1:3)
    counts[r,] <- counts[r,] / sum(counts[r,])
  
  # Output the figure.
  xlabels <- paste0(genotypes,"\n(",table(geno),")")
  return(barchart(counts,horizontal = FALSE,lwd = 1,box.ratio = 1,
                  col = c("black","white"),border = "black",
                  xlab = "genotype",ylab = paste0("p(",phenotype," = 1)"),
                  scales = list(x = list(labels = xlabels,tck = 0),
                                y = list(tck = 0.5)),
                  par.settings = list(par.main.text = list(cex = 1,font = 1),
                                      layout.heights = list(axis.top = 0)),
                  main = paste0(marker," (chr ",map$chr,", ",
                                round(map$pos)," Mb)")))
}
