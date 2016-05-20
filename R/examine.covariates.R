# A script to show scatterplots of phenotype versus covariate, and to
# calculate the proportion of variance in a phenotype interest that is
# explained by various candidate covariates (e.g., body weight).
library(lattice)
source("misc.R")
source("read.data.R")
source("data.manip.R")
source("plotting.tools.R")

# SCRIPT PARAMETERS
# -----------------
# I organize the analyses into "clusters" of phenotypes: muscle and
# bone traits ("muscle+bone"), body weights ("bw"), other
# physiological traits ("physio"), fear conditioning phenotypes
# ("fc"), methamphetamine sensitivity phenotypes ("meth"), prepulse
# inhibition phenotypes ("ppi"), and others (see below).
cluster <- "muscle+bone"

if (cluster == "muscle+bone") {

  # MUSCLE + BONE TRAITS
  # --------------------
  # Examine the muscle weight phenotypes (TA, EDL, soleus, plantaris,
  # gastroc), the first muscle weight principal component ("mwpc"),
  # tibia length ("tibia"), and bone-mineral density ("BMD").
  panels <-
    list(TA1      = list(pheno="TA",       cov="sacweight",col="dodgerblue"),
         TA2      = list(pheno="TA",       cov="tibia",    col="dodgerblue"),
         EDL1     = list(pheno="EDL",      cov="sacweight",col="dodgerblue"),
         EDL2     = list(pheno="EDL",      cov="tibia",    col="dodgerblue"),
         soleus1  = list(pheno="soleus",   cov="sacweight",col="dodgerblue"),
         soleus2  = list(pheno="soleus",   cov="tibia",    col="dodgerblue"),
         plant1   = list(pheno="plantaris",cov="sacweight",col="dodgerblue"),
         plant2   = list(pheno="plantaris",cov="tibia",    col="dodgerblue"),
         gastroc1 = list(pheno="gastroc",  cov="sacweight",col="dodgerblue"),
         gastroc2 = list(pheno="gastroc",  cov="tibia",    col="dodgerblue"),
         mwpc1    = list(pheno="mwpc",     cov="sacweight",col="darkblue"),
         mwpc2    = list(pheno="mwpc",     cov="tibia",    col="darkblue"),
         BMD1     = list(pheno="BMD",      cov="sacweight",col="darkorange"),
         BMD2     = list(pheno="BMD",      cov="tibia",    col="darkorange"),
         BMD3     = list(pheno="BMD",      cov="gastroc",  col="darkorange"),
         tibia    = list(pheno="tibia",    cov="sacweight",col="darkorange"))

  # The panels are arranged on a 4 x 4 grid.
  trellis.device(height = 6,width = 6.25)
  panel.layout <- list(nrow = 4,ncol = 4)
} else if (cluster == "bw") {

  # LOOK MORE CLOSELY AT BODY WEIGHTS
  # ---------------------------------
  panels <-
    list(panel1 = list(pheno="bw1",      cov="bw0",      col="darkblue"),
         panel2 = list(pheno="bw2",      cov="bw1",      col="darkblue"),
         panel3 = list(pheno="bw3",      cov="bw2",      col="darkblue"),
         panel4 = list(pheno="sacweight",cov="bw3",      col="darkblue"),
         panel5 = list(pheno="PPIweight",cov="sacweight",col="darkblue"))

  # The panels are arranged on a 2 x 3 grid.
  trellis.device(height = 3,width = 4.75)
  panel.layout <- list(nrow = 2,ncol = 3)
} else if (cluster == "physio") {

  # OTHER PHYSIOLOGICAL TRAITS
  # --------------------------
  # Examine fasting glucose levels ("fastglucose"), body weights,
  # testes weight, and tail length.
  panels <-
  list(bw0      = list(pheno="bw0",         cov="glucoseage",col="darkgreen"),
       bw1      = list(pheno="bw1",         cov="methage",   col="darkgreen"),
       ppiwt    = list(pheno="PPIweight",   cov="PPIage",    col="darkgreen"),
       testes   = list(pheno="testesweight",cov="sacweight", col="firebrick"),
       glucose1 = list(pheno="fastglucose", cov="bw0",       col="firebrick"),
       glucose2 = list(pheno="fastglucose", cov="glucoseage",col="firebrick"),
       tail1    = list(pheno="taillength",  cov="glucoseage",col="darkorange"),
       tail2    = list(pheno="taillength",  cov="bw0",       col="darkorange"))

  # The panels are arranged on a 3 x 3 grid.
  trellis.device(height = 5,width = 5)
  panel.layout <- list(nrow = 3,ncol = 3)
} else if (cluster == "fc") {

  # FEAR CONDITIONING TRAITS
  # ------------------------
  panels <- list(
AvContextD2a = list(pheno="AvContextD2",   cov="PreTrainD1",col="dodgerblue"),
AvContextD3a = list(pheno="AvAltContextD3",cov="PreTrainD1",col="dodgerblue"),
AvToneD3a    = list(pheno="AvToneD3",      cov="PreTrainD1",col="dodgerblue"),
AvContextD2b = list(pheno="AvContextD2",   cov="AvToneD1",col="darkblue"),
AvContextD3b = list(pheno="AvAltContextD3",cov="AvToneD1",col="darkblue"),
AvToneD3b    = list(pheno="AvToneD3",      cov="AvToneD1",col="darkblue"),
extinction   = list(pheno="D3.360",        cov="D3.180",col="forestgreen"))
  
  # The panels are arranged on a 3 x 2 grid.
  trellis.device(height = 4.5,width = 5)
  panel.layout <- list(nrow = 3,ncol = 3)
} else if (cluster == "meth") {

  # METHAMPHETAMINE SENSITIVITY TRAITS
  # ----------------------------------
  panels <- list(d3.0to15a   = list(pheno = "D3totaldist0to15",
                                    cov   = "D1totaldist0to15",
                                    col   = "darkblue"),
                 d3.15to30a  = list(pheno = "D3totaldist15to30",
                                    cov   = "D1totaldist15to30",
                                    col   = "darkblue"),
                 d3.0to30a   = list(pheno = "D3totaldist0to30",
                                    cov   = "D1totaldist0to30",
                                    col   = "darkblue"),
                 d3.0to15b   = list(pheno = "D3totaldist0to15",
                                    cov   = "D2totaldist0to15",
                                    col   = "dodgerblue"),
                 d3.15to30b  = list(pheno = "D3totaldist15to30",
                                    cov   = "D2totaldist15to30",
                                    col   = "dodgerblue"),
                 d3.0to30b   = list(pheno = "D3totaldist0to30",
                                    cov   = "D2totaldist0to30",
                                    col   = "dodgerblue"))

  # The panels are arranged on a 3 x 2 grid.
  trellis.device(height = 4.75,width = 3.5)
  panel.layout <- list(nrow = 3,ncol = 2)
} else if (cluster == "meth2") {

  # MORE PHENOTYPES FROM METHAMPHETAMINE SENSITIVITY TESTS
  # ------------------------------------------------------
  panels <- list(
    D1hact0to30 = list(pheno = "D1hact0to30",cov="bw1",col="darkblue"),
    D2hact0to30 = list(pheno = "D2hact0to30",cov="bw1",col="darkblue"),
    D3hact0to30 = list(pheno = "D3hact0to30",cov="bw1",col="darkblue"),
    D1vact0to30 = list(pheno = "D1vact0to30",cov="bw1",col="magenta"),
    D2vact0to30 = list(pheno = "D2vact0to30",cov="bw1",col="magenta"),
    D3vact0to30 = list(pheno = "D3vact0to30",cov="bw1",col="magenta"))

  # The panels are arranged on a 3 x 2 grid.
  trellis.device(height = 4.75,width = 3.5)
  panel.layout <- list(nrow = 3,ncol = 2)
} else if (cluster == "ppi") {

  # PREPULSE INHIBITION (PPI) PHENOTYPES.
  panels <- list(
    pp3avg     = list(pheno="pp3avg",    cov="PPIweight",col="indianred"),
    pp6avg     = list(pheno="pp6avg",    cov="PPIweight",col="indianred"),
    pp12avg    = list(pheno="pp12avg",   cov="PPIweight",col="indianred"),
    pp3PPIavg  = list(pheno="pp3PPIavg", cov="PPIweight",col="darkblue"),
    pp6PPIavg  = list(pheno="pp6PPIavg", cov="PPIweight",col="darkblue"),
    pp12PPIavg = list(pheno="pp12PPIavg",cov="PPIweight",col="darkblue"),
    startle    = list(pheno="startle",   cov="PPIweight",col="dodgerblue"),
    avgnostim  = list(pheno="avgnostim", cov="PPIweight",col="darkgreen"),
    p120b4     = list(pheno="p120b4",    cov="PPIweight",col="orangered"),
    habit      = list(pheno="p120b4",    cov="p120b1",   col="orangered"))    

  # The panels are arranged on a 3 x 4 grid.
  trellis.device(height = 4.5,width = 6.25)
  panel.layout <- list(nrow = 3,ncol = 4)
} else if (cluster == "bmi") {
    
    # MORE PHENOTYPES FROM METHAMPHETAMINE SENSITIVITY TESTS
    # ------------------------------------------------------
    panels <- list(
    sac.bmi.tibia1=list(pheno="sacwt.bmi.tibia",cov="sacweight",col="blue"),
    sac.bmi.tibia2=list(pheno="sacwt.bmi.tibia",cov="tibia",col="blue"),
    sac.bmi.tibia4=list(pheno="sacwt.bmi.tibia",cov="bw0",col="blue"))
    
    # The panels are arranged on a 2 x 2 grid.
    trellis.device(height = 4.5,width = 2)
    panel.layout <- list(nrow = 3,ncol = 1)
}

# Set up the graphics device.
trellis.par.set(list(fontsize = list(text = 8),
                     layout.widths = list(right.padding = -1),
                     layout.heights = list(top.padding = -1,
                                           bottom.padding = 0)))

# Set up the location of each panel on the grid.
panel.layout <- with(panel.layout,
                     create.grid.layout(nrow,ncol,names(panels)))

# Load the phenotype data.
phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- read.pheno("pheno.csv")
pheno      <- prepare.pheno(pheno)

# Jitter glucoseage for the scatterplots.
n     <- nrow(pheno)
pheno <- transform(pheno,glucoseage = glucoseage + rnorm(n,sd = 0.5))

# Show the histogram plots.
for (panel in names(panels)) {

  # Get the phenotype and covariate to investigate, and the colour for
  # drawing the points in the scatterplot.
  r         <- panels[[panel]]
  phenotype <- r$pheno
  covariate <- r$cov
  panel.col <- r$col

  # Get the phenotype (Y) and covariate (X) data.
  data        <- pheno[c(phenotype,covariate)]
  names(data) <- c("y","x")
  
  # Fit the linear regression for the phenotype given the covariate,
  # and get the proportion of variance in the phenotype explained by
  # the covariate.
  model <- lm(y ~ x,data)
  mu    <- coef(model)["(Intercept)"]
  b     <- coef(model)["x"]
  pve   <- summary(model)$r.squared
  
  # Show the scatterplot of covariate vs phenotype, with the best
  # fit regression line.
  print(xyplot(y ~ x,data,ylab = phenotype,
               xlab=paste0(covariate," (PVE = ",round(100*pve,digits=2),"%)"),
               scales = list(x = list(tck = 0.6),y = list(tck = 0.6)),
               panel = function (...) {
                 panel.xyplot(...,pch = 20,cex = 0.4,col = panel.col);
                 panel.abline(mu,b,lwd = 2,col = "limegreen");
             }),
        split = grid.split(panel.layout,panel),
        more = TRUE)
}
