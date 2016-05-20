# A script to show scatterplots of phenotype vs covariate, and
# calculate the proportion of variance in a phenotype of interest that
# is explained by variance candidate binary covariates (e.g., testing
# apparatus).
library(lattice)
library(Hmisc)
source("misc.R")
source("read.data.R")
source("data.manip.R")
source("plotting.tools.R")

# SCRIPT PARAMETERS
# -----------------
# I organize the analyses into "clusters" of phenotypes: muscle and
# bone traits ("muscle+bone"), other physiological traits ("physio"),
# fear conditioning phenotypes ("fc"), methamphetamine sensitivity
# phenotypes ("meth"), and prepulse inhibition phenotypes ("ppi"), and
# others (see below).
cluster <- "muscle+bone"

if (cluster == "muscle+bone") {

  # MUSCLE + BONE TRAITS
  # --------------------
  # Examine the muscle weight phenotypes (TA, EDL, soleus, plantaris,
  # gastroc), the first muscle weight principal component ("mwpc"),
  # tibia length ("tibia"), and bone-mineral density ("BMD").
  panels <-
    list(TA      = list(pheno="TA",       cov="SW16"),
         EDL     = list(pheno="EDL",      cov="SW16"),
         soleus  = list(pheno="soleus",   cov="SW16"),
         plant   = list(pheno="plantaris",cov="SW16"),
         gastroc = list(pheno="gastroc",  cov="SW16"),
         mwpc    = list(pheno="mwpc",     cov="SW16"),
         BMD     = list(pheno="BMD",      cov="SW16"),
         tibia1  = list(pheno="tibia",    cov="SW16"),
         tibia2  = list(pheno="tibia",    cov="SW6"))

  # The panels are arranged on a 3 x 3 grid.
  trellis.device(height = 4.25,width = 4)
  panel.layout <- list(nrow = 3,ncol = 3)
} else if (cluster == "physio") {

  # OTHER PHYSIOLOGICAL TRAITS
  # --------------------------
  # Examine fasting glucose levels ("fastglucose"), body weights
  # testes weight, and tail length.
  panels <-
    list(sacweight    = list(pheno="bw0",   cov="SW17"),
         fastglucose1 = list(pheno="fastglucose", cov="SW1"),
         fastglucose2 = list(pheno="fastglucose", cov="SW11"),
         taillength1  = list(pheno="taillength",  cov="SW3"),
         taillength2  = list(pheno="taillength",  cov="SW4"),
         taillength3  = list(pheno="taillength",  cov="SW19"),
         taillength4  = list(pheno="taillength",  cov="SW20"),
         taillength5  = list(pheno="taillength",  cov="SW22"),
         taillength6  = list(pheno="taillength",  cov="SW24"))
    
  # The panels are arranged on a 3 x 3 grid.
  trellis.device(height = 4.5,width = 4.5)
  panel.layout <- list(nrow = 3,ncol = 3)
} else if (cluster == "fc") {

  # FEAR CONDITIONING TRAITS
  # ------------------------
  panels <-
    list(
         # Covariate = FC box #1.
         d1pretrain1 = list(pheno="PreTrainD1",    cov="FCbox1"),
         d1tone1     = list(pheno="AvToneD1",      cov="FCbox1"),
         d2ctxt1     = list(pheno="AvContextD2",   cov="FCbox1"),
         d3altctxt1  = list(pheno="AvAltContextD3",cov="FCbox1"),
         d3tone1     = list(pheno="AvToneD3",      cov="FCbox1"),

         # Covariate = FC box #2.
         d1pretrain2 = list(pheno="PreTrainD1",    cov="FCbox2"),
         d1tone2     = list(pheno="AvToneD1",      cov="FCbox2"),
         d2ctxt2     = list(pheno="AvContextD2",   cov="FCbox2"),
         d3altctxt2  = list(pheno="AvAltContextD3",cov="FCbox2"),
         d3tone2     = list(pheno="AvToneD3",      cov="FCbox2"),

         # Covariate = FC box #3.
         d1pretrain3 = list(pheno="PreTrainD1",    cov="FCbox3"),
         d1tone3     = list(pheno="AvToneD1",      cov="FCbox3"),
         d2ctxt3     = list(pheno="AvContextD2",   cov="FCbox3"),
         d3altctxt3  = list(pheno="AvAltContextD3",cov="FCbox3"),
         d3tone3     = list(pheno="AvToneD3",      cov="FCbox3"),

         # Covariate = FC box #4.
         d1pretrain4 = list(pheno="PreTrainD1",    cov="FCbox4"),
         d1tone4     = list(pheno="AvToneD1",      cov="FCbox4"),
         d2ctxt4     = list(pheno="AvContextD2",   cov="FCbox4"),
         d3altctxt4  = list(pheno="AvAltContextD3",cov="FCbox4"),
         d3tone4     = list(pheno="AvToneD3",      cov="FCbox4"),

         # Covariate = round SW17.
         d1pretrain5 = list(pheno="PreTrainD1",    cov="SW17"),
         d1tone5     = list(pheno="AvToneD1",      cov="SW17"),
         d2ctxt5     = list(pheno="AvContextD2",   cov="SW17"),
         d3altctxt5  = list(pheno="AvAltContextD3",cov="SW17"),
         d3tone5     = list(pheno="AvToneD3",      cov="SW17"),
        
         # Day 2 + 3 traits freezing, controlling for Day 1 measurements.
         d2ctxt6    = list(pheno="AvContextD2",   cov="pretrainbin"),
         d3altctxt6 = list(pheno="AvAltContextD3",cov="pretrainbin"),
         d3tone6    = list(pheno="AvToneD3",      cov="pretrainbin"),
         d2ctxt7    = list(pheno="AvContextD2",   cov="d1tonebin"),
         d3altctxt7 = list(pheno="AvAltContextD3",cov="d1tonebin"),
         d3tone7    = list(pheno="AvToneD3",      cov="d1tonebin"))

  # The panels are arranged on a 5 x 7 grid.
  trellis.device(height = 6.5,width = 8.5)
  panel.layout <- list(nrow = 5,ncol = 7)
} else if (cluster == "meth") {

  # METHAMPHETAMINE SENSITIVITY TRAITS
  # ----------------------------------
  panels <-
    list(d1.0to15a  = list(pheno="D1totaldist0to15", cov="methcage7"),
         d1.15to30a = list(pheno="D1totaldist15to30",cov="methcage7"),
         d2.0to15a  = list(pheno="D2totaldist0to15", cov="methcage7"),
         d2.15to30a = list(pheno="D2totaldist15to30",cov="methcage7"),
         d3.0to15a  = list(pheno="D3totaldist0to15", cov="methcage7"),
         d3.15to30a = list(pheno="D3totaldist15to30",cov="methcage7"))

  # The panels are arranged on a 2 x 3 grid.
  trellis.device(height = 3,width = 4.5)
  panel.layout <- list(nrow = 2,ncol = 3)
} else if (cluster == "meth2") {

  # MORE PHENOTYPES FROM METHAMPHETAMINE SENSITIVITY TESTS
  # ------------------------------------------------------
  panels <- list(
    D1hact0to30 = list(pheno = "D1hact0to30",cov = "methcage8"),
    D2hact0to30 = list(pheno = "D2hact0to30",cov = "methcage8"),
    D3hact0to30 = list(pheno = "D3hact0to30",cov = "methcage8"),
    D1vact0to30 = list(pheno = "D1vact0to30",cov = "methcage8"),
    D2vact0to30 = list(pheno = "D2vact0to30",cov = "methcage8"),
    D3vact0to30 = list(pheno = "D3vact0to30",cov = "methcage8"))
  
  # The panels are arranged on a 2 x 3 grid.
  trellis.device(height = 3,width = 4.5)
  panel.layout <- list(nrow = 2,ncol = 3)
} else if (cluster == "ppi") {

  # PREPULSE INHIBITION (PPI) PHENOTYPES
  # ------------------------------------
  panels <-
    list(
         pp3avg      = list(pheno="pp3avg",    cov="PPIbox3"),
         pp6avg      = list(pheno="pp6avg",    cov="PPIbox3"),
         pp12avg     = list(pheno="pp12avg",   cov="PPIbox3"),
         pp3PPIavg   = list(pheno="pp3PPIavg", cov="PPIbox3"),
         pp6PPIavg   = list(pheno="pp6PPIavg", cov="PPIbox3"),
         pp12PPIavg  = list(pheno="pp12PPIavg",cov="PPIbox3"),
         startle     = list(pheno="startle",   cov="PPIbox3"),
         avgnostim   = list(pheno="avgnostim", cov="PPIbox3"),
         p120b1      = list(pheno="p120b1",    cov="PPIbox3"),
         p120b4      = list(pheno="p120b4",    cov="PPIbox3"))

  # The panels are arranged on a 3 x 3 grid.
  trellis.device(height = 4,width = 5)
  panel.layout <- list(nrow = 3,ncol = 4)
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

# Create binary covariates from some of the categorical phenotypes.
pheno <-
  cbind(pheno,
        binary.from.categorical(pheno$FCbox,paste0("FCbox",1:4)),
        binary.from.categorical(pheno$methcage,paste0("methcage",1:12)),
        binary.from.categorical(pheno$PPIbox,paste0("PPIbox",1:5)),
        binary.from.categorical(pheno$round,paste0("SW",1:25)))

# Show the box-percentile plots.
for (panel in names(panels)) {

  # Get the phenotype and covariate to investigate.
  r         <- panels[[panel]]
  phenotype <- r$pheno
  covariate <- r$cov

  # Get the phenotype (Y) and covariate (X) data, then convert the
  # binary covariate from a factor to an integer (with possible values
  # of 0 and 1) so that we can use it in a linear model of the
  # phenotype.
  data        <- pheno[c(phenotype,covariate)]
  names(data) <- c("y","x")
  data        <- transform(data,x = binfactor2num(x))

  # Fit the linear regression for the phenotype given the covariate,
  # and get the proportion of variance in the phenotype explained by
  # the covariate.
  model <- lm("y ~ x",data)
  pve   <- summary(model)$r.squared
      
  # Show the distribution of the phenotype conditioned on the binary
  # covariate using a "box-percentile" plot.
  print(bwplot(formula(paste(covariate,"~",phenotype)),pheno,
               probs = c(0.01,0.05,0.125,0.25,0.375),
               means = FALSE,nout = 0.003,panel = panel.bpplot,
               scales = list(x = list(tck = 0.5),y = list(tck = 0)),
               scat1d.opts = list(lwd = 4,tfrac = 0.01,col = "orangered"),
               xlab=paste0(phenotype," (PVE = ",round(100*pve,digits=2),"%)"),
               ylab = covariate),
        split = grid.split(panel.layout,panel),
        more = TRUE)
}
