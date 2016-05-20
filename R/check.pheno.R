# This is a small script to check the whether the observed quantiles
# for each phenotype, conditioned on different sets of covariates,
# match what we would expect under the normal distribution.
library(lattice)
source("misc.R")
source("read.data.R")
source("data.manip.R")

# Load the phenotype data.
pheno <- read.pheno("pheno.csv")
pheno <- prepare.pheno(pheno)

# Create binary covariates from some of the categorical phenotypes.
pheno <- cbind(pheno,
               binary.from.categorical(pheno$FCbox,paste0("FCbox",1:4)),
               binary.from.categorical(pheno$methcage,paste0("methcage",1:12)),
               binary.from.categorical(pheno$PPIbox,paste0("PPIbox",1:5)),
               binary.from.categorical(pheno$round,paste0("SW",1:25)))

# PHENOTYPE = TA
# --------------
# TA with no additional covariates.
trellis.device(height = 2.5,width = 3)
trellis.par.set(list(fontsize = list(text = 9)))
r <- resid(lm(TA ~ SW16,pheno,na.action = na.exclude))
r[r < (-20) | r > 20] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of TA | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# TA given tibia length.
r <- resid(lm(TA ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-18)] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of TA | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# TA given tibia length and body weight.
r <- resid(lm(TA ~ SW16 + tibia + sacweight,pheno,na.action = na.exclude))
r[r < (-15) | r > 15] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of TA | SW16 + tibia + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = EDL
# ---------------
# EDL with no additional covariates.
r <- resid(lm(EDL ~ SW16,pheno,na.action = na.exclude))
r[r < (-4) | r > 4] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of EDL | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# EDL given tibia length.
r <- resid(lm(EDL ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-5) | r > 4] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of EDL | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# EDL given tibia length and body weight.
r <- resid(lm(EDL ~ SW16 + tibia + sacweight,pheno,na.action = na.exclude))
r[r < (-4) | r > 4] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of EDL | SW16 + tibia + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = SOLEUS
# ------------------
# Soleus with no additional covariates.
r <- resid(lm(soleus ~ SW16,pheno,na.action = na.exclude))
r[r < (-4) | r > 4] <- NA
print(histogram(r,col = "gold",border = "gold",nint = 24,
                xlab = "residual of soleus | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Soleus given tibia length.
r <- resid(lm(soleus ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-4) | r > 4] <- NA
print(histogram(r,col = "gold",border = "gold",nint = 24,
                xlab = "residual of soleus | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Soleus given tibia length and body weight.
r <- resid(lm(soleus ~ SW16 + tibia + sacweight,pheno,na.action = na.exclude))
r[r < (-3.5) | r > 3.5] <- NA
print(histogram(r,col = "gold",border = "gold",nint = 24,
                xlab = "residual of soleus | SW16 + tibia + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PLANTARIS
# ---------------------
# Plantaris with no additional covariates.
r <- resid(lm(plantaris ~ SW16,pheno,na.action = na.exclude))
r[r < (-8) | r > 8] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 24,
                xlab = "residual of plantaris | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Plantaris given tibia length.
r <- resid(lm(plantaris ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-9) | r > 8] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 24,
                xlab = "residual of plantaris | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Plantaris given tibia length and body weight.
r <- resid(lm(plantaris ~ SW16 + tibia + sacweight,pheno,
              na.action = na.exclude))
r[r < (-8) | r > 7] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 24,
                xlab = "residual of plantaris | SW16 + tibia + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = GASTROC
# -------------------
# Gastroc with no additional covariates.
r <- resid(lm(gastroc ~ SW16,pheno,na.action = na.exclude))
r[r < (-50) | r > 50] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of gastroc | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Gastroc given tibia length.
r <- resid(lm(gastroc ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-40) | r > 50] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of gastroc | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Gastroc given tibia length and body weight.
r <- resid(lm(gastroc ~ SW16 + tibia + sacweight,pheno,
              na.action = na.exclude))
r[r < (-40) | r > 40] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of gastroc | SW16 + tibia + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = MWPC
# ----------------
# MWPC with no additional covariates.
r <- resid(lm(mwpc ~ SW16,pheno,na.action = na.exclude))
r[r < (-6) | r > 5] <- NA
print(histogram(r,col = "darkorchid",border = "darkorchid",nint = 24,
                xlab = "residual of mwpc | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# MWPC given body weight.
r <- resid(lm(mwpc ~ SW16 + sacweight,pheno,na.action = na.exclude))
r[r < (-4) | r > 4.5] <- NA
print(histogram(r,col = "darkorchid",border = "darkorchid",nint = 24,
                xlab = "residual of mwpc | SW16 + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# MWPC given tibia length.
r <- resid(lm(mwpc ~ SW16 + tibia,pheno,na.action = na.exclude))
r[r < (-5) | r > 5] <- NA
print(histogram(r,col = "darkorchid",border = "darkorchid",nint = 24,
                xlab = "residual of mwpc | SW16 + tibia"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = BMD
# ---------------
# BMD with no additional covariates.
r <- resid(lm(BMD ~ SW16,pheno,na.action = na.exclude))
r[r > 0.14] <- NA
print(histogram(r,col = "cornflowerblue",border = "cornflowerblue",
                nint = 24,xlab = "residual of BMD | SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = TIBIA
# -----------------
# Tibia with no additional covariates.
r <- resid(lm(tibia ~ SW6 + SW16,pheno,na.action = na.exclude))
r[r < (-1.4)] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 24,
                xlab = "residual of tibia | SW6 + SW16"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Tibia given body weight.
r <- resid(lm(tibia ~ SW6 + SW16 + sacweight,pheno,na.action = na.exclude))
r[r < (-1.5)] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 24,
                xlab = "residual of tibia | SW6 + SW16 + sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = SACWEIGHT
# ---------------------
# Sacweight with no additional covariates.
r <- resid(lm(sacweight ~ SW17,pheno,na.action = na.exclude))
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 24,
                xlab = "residual of sacweight | SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = BW0
# ---------------
# bw0 given glucose age.
r <- resid(lm(bw0 ~ glucoseage + SW17,pheno,na.action = na.exclude))
r[r < (-8.5) | r > 8.5] <- NA
print(histogram(r,col = "yellowgreen",border = "yellowgreen",nint = 24,
                xlab = "residual of bw0 | glucoseage + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = BW1
# ---------------
r <- resid(lm(bw1 ~ methage + SW17,pheno,na.action = na.exclude))
r[r < (-9) | r > 10] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 24,
                xlab = "residual of bw1 | methage + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PPIWEIGHT
# ---------------------
r <- resid(lm(PPIweight ~ SW17,pheno,na.action = na.exclude))
print(histogram(r,col = "limegreen",border = "limegreen",nint = 24,
                xlab = "residual of PPIweight | SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = FASTGLUCOSE
# -----------------------
# Fasting glucose levels given bw0.
r <- resid(lm(fastglucose ~ SW1 + SW11 + bw0,pheno,na.action = na.exclude))
r[r < (-60) | r > 60] <- NA
print(histogram(r,col = "indianred",border = "indianred",nint = 24,
                xlab = "residual of fastglucose | SW1 + SW11 + bw0"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = TESTESWEIGHT
# ------------------------
# Testes weight given sacweight.
r <- resid(lm(testesweight ~ sacweight,pheno,na.action = na.exclude))
r[r < (-0.075)] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of testesweight | sacweight"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = TAILLENGTH
# ----------------------
# Tail length given age, body weight, and binary indicators for rounds
# SW3, SW4, SW19, SW20, SW22 and SW44.
r <- resid(lm(taillength ~ bw0 + glucoseage + SW3 + SW4 + SW19 +
              SW20 + SW22 + SW24,pheno,na.action = na.exclude))
print(histogram(r,col = "darkorange",border = "darkorange",nint = 20,
                xlab = "residual of taillength | bw0 + glucoseage + round"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PRETRAIND1
# ----------------------
r <- resid(lm(PreTrainD1 ~ FCbox1+FCbox2+FCbox3 + SW10+SW16+SW17+SW20,pheno,
              na.action = na.exclude))
print(histogram(r,col = "darkcyan",border = "darkcyan",nint = 24,
                xlab = "PreTrainD1 | FCbox + SW--"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = AVTONED1
# ----------------------
r <- resid(lm(AvToneD1 ~ FCbox1 + FCbox2 + FCbox3 + SW10 + SW7 + SW14 + SW20,
              pheno,na.action = na.exclude))
print(histogram(r,col = "darkcyan",border = "darkcyan",nint = 24,
                xlab = "AvToneD1 | FCbox + SW--"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = AVCONTEXTD2
# -----------------------
# Freezing on Day 2 in response to test chamber conditioned on which
# FC box is used in the tests, and a binary indicator for round SW17.
r <- resid(lm(AvContextD2 ~ FCbox1 + FCbox2 + FCbox3 + SW17,pheno,
              na.action = na.exclude))
print(histogram(r,col = "darkcyan",border = "darkcyan",nint = 24,
                xlab = "resid(AvContextD2) | FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Freezing on Day 2 in response to test chamber ("AvContextD2")
# conditioned on freezing on Day 1 after exposure to conditioned
# stimulus ("AvToneD1"), and conditioned on which FC box is used in
# the tests, plus a binary indicator for round SW17.
r <- resid(lm(AvContextD2 ~ AvToneD1 + FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
print(histogram(r,col = "darkcyan",border = "darkcyan",nint = 24,
                xlab = "resid(AvContextD2) | AvToneD1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Freezing on Day 2 in response to test chamber ("AvContextD2")
# conditioned on Day 1 freezing measures ("PreTrainD1" and
# "AvToneD1"), and conditioned on which FC box is used in the tests,
# plus a binary indicator for round SW17.
r <- resid(lm(AvContextD2 ~ PreTrainD1 + AvToneD1 +
              FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
print(histogram(r,col = "darkcyan",border = "darkcyan",nint = 24,
                xlab = "resid(AvContextD2) | D1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = AVALTCONTEXTD3
# --------------------------
# Proportion of freezing in the altered context on Day 3
# ("AvAltContextD3") conditioned on freezing on Day 1 after exposure to
# conditioned stimulus ("AvToneD1"), and conditioned on which FC box
# is used in the tests, plus a binary indicator for round SW17.
r <- resid(lm(AvAltContextD3 ~ AvToneD1 + FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
r[r > 1] <- NA
print(histogram(r,col = "darkred",border = "darkred",nint = 24,
                xlab = "resid(AvAltContextD3) | AvToneD1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Proportion of freezing in the altered context on Day 3
# ("AvAltContextD3") conditioned on Day 1 freezing measurements
# ("PreTrainD1" and "AvToneD1"), and conditioned on which FC box is
# used in the tests, plus a binary indicator for round SW17.
r <- resid(lm(AvAltContextD3 ~ PreTrainD1 + AvToneD1 +
              FCbox1 + FCbox2 + FCbox3 + SW17, 
              pheno,na.action = na.exclude))
r[r > 0.95] <- NA
print(histogram(r,col = "darkred",border = "darkred",nint = 24,
                xlab = "resid(AvAltContextD3) | D1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = AVTONED3
# --------------------
# Freezing on Day 3 after receiving the tones conditioned on which FC
# box is used in the tests, and a binary indicator for round SW17.
r <- resid(lm(AvToneD3 ~ FCbox1 + FCbox2 + FCbox3 + SW17,pheno,
              na.action = na.exclude))
print(histogram(r,col = "tomato",border = "tomato",nint = 24,
                xlab = "residual of AvToneD3 | FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Freezing on Day 3 after receiving the tones ("AvToneD3") conditioned
# on freezing on Day 1 after exposure to conditioned stimulus
# ("AvToneD1"), and conditioned on which FC box is used in the tests,
# plus a binary indicator for round SW17.
r <- resid(lm(AvToneD3 ~ AvToneD1 + FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
print(histogram(r,col = "tomato",border = "tomato",nint = 24,
                xlab = "resid(AvToneD3) | AvToneD1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Freezing on Day 3 after receiving the tones ("AvToneD3") conditioned
# on Day 1 freezing measures ("PreTrainD1" and "AvToneD1"), and
# conditioned on which FC box is used in the tests, plus a binary
# indicator for round SW17.
r <- resid(lm(AvToneD3 ~ PreTrainD1 + AvToneD1 +
              FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
print(histogram(r,col = "tomato",border = "tomato",nint = 24,
                xlab = "resid(AvToneD3) | D1 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3.360
# ------------------
# Proportion of freezing on day 3 after the 360-second tone,
# conditioned on proportion of freezing after the 180-second tone, and
# conditioned on which FC box is used in the tests, plus a binary
# indicator for round SW17.
r <- resid(lm(D3.360 ~ D3.180 + FCbox1 + FCbox2 + FCbox3 + SW17,
              pheno,na.action = na.exclude))
print(histogram(r,col = "yellowgreen",border = "yellowgreen",nint = 24,
                xlab = "residual of D3.360 | D3.180 + FCbox + SW17"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTALDIST0TO15
# ----------------------------
# Locomotor response on Day 1 in 0-15 minute interval with no
# additional covariates.
r <- resid(lm(D1totaldist0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-2000) | r > 2200] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of D1totaldist0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTALDIST15TO30
# -----------------------------
# Locomotor response on Day 1 in 15-30 minute interval with no
# additional covariates.
r <- resid(lm(D1totaldist15to30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 1900] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of D1totaldist15to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTALDIST0TO30
# ----------------------------
# Locomotor response on Day 1 in 0-30 minute interval with no
# additional covariates.
r <- resid(lm(D1totaldist0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-3500) | r > 4000] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of D1totaldist0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTALDIST0TO15
# ----------------------------
# Locomotor response on Day 2 in 0-15 minute interval with no
# additional covariates.
r <- resid(lm(D2totaldist0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-2000) | r > 2500] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of D2totaldist0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTALDIST15TO30
# -----------------------------
# Locomotor response on Day 2 in 15-30 minute interval with no
# additional covariates.
r <- resid(lm(D2totaldist15to30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 2000] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of D2totaldist15to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTALDIST0TO30
# ----------------------------
# Locomotor response on Day 2 in 0-30 minute interval with no
# additional covariates.
r <- resid(lm(D2totaldist0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 4500] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of D2totaldist0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTALDIST0TO15
# ----------------------------
# Locomotor response on Day 3 in 0-15 minute interval with no
# additional covariates.
r <- resid(lm(D3totaldist0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r > 8500] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of D3totaldist0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Locomotor response on Day 3 in 0-15 minute interval conditioned on
# the locomotor response on Day 2 over the same interval.
r <- resid(lm(D3totaldist0to15 ~ D2totaldist0to15 + methcage7,pheno,
              na.action = na.exclude))
r[r < (-7000) | r > 8000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "resid(D3totaldist0to15) | D2 + methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTALDIST15TO30
# -----------------------------
# Locomotor response on Day 3 in 15-30 minute interval with no
# additional covariates.
r <- resid(lm(D3totaldist15to30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 12000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of D3totaldist15to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Locomotor response on Day 3 in 15-30 minute interval conditioned on
# the locomotor response on Day 2 over the same interval.
r <- resid(lm(D3totaldist15to30 ~ D2totaldist15to30 + methcage7,pheno,
              na.action = na.exclude))
r[r > 10000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "resid(D3totaldist15to30) | D2 + methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTALDIST0TO30
# ----------------------------
# Locomotor response on Day 3 in 0-30 minute interval with no
# additional covariates.
r <- resid(lm(D3totaldist0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 20000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "residual of D3totaldist0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# Locomotor response on Day 3 in 0-30 minute interval conditioned on
# the locomotor response on Day 2 over the same interval.
r <- resid(lm(D3totaldist0to30 ~ D2totaldist0to30 + methcage7,pheno,
              na.action = na.exclude))
r[r > 17000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 24,
                xlab = "resid(D3totaldist0to30) | D2 + methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST5
# ----------------------
# Locomotor response on Day 1 over 0-5 minute interval.
r <- resid(lm(D1TOTDIST5 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-1000) | r > 1000] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist5 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST10
# -----------------------
# Locomotor response on Day 1 over 5-10 minute interval.
r <- resid(lm(D1TOTDIST10 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-750) | r > 750] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist10 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST15
# -----------------------
# Locomotor response on Day 1 over 10-15 minute interval.
r <- resid(lm(D1TOTDIST15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-750) | r > 750] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST20
# -----------------------
# Locomotor response on Day 1 over 15-20 minute interval.
r <- resid(lm(D1TOTDIST20 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-750) | r > 750] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist20 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST25
# -----------------------
# Locomotor response on Day 1 over 20-25 minute interval.
r <- resid(lm(D1TOTDIST25 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-750) | r > 750] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist25 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1TOTDIST30
# -----------------------
# Locomotor response on Day 1 over 25-30 minute interval.
r <- resid(lm(D1TOTDIST30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 700] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of d1totdist30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST5
# ----------------------
# Locomotor response on Day 2 over 0-5 minute interval.
r <- resid(lm(D2TOTDIST5 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-1250) | r > 1250] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist5 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST10
# -----------------------
# Locomotor response on Day 2 over 5-10 minute interval.
r <- resid(lm(D2TOTDIST10 ~ methcage7,pheno,na.action = na.exclude))
r[r > 1000] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist10 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST15
# -----------------------
# Locomotor response on Day 2 over 10-15 minute interval.
r <- resid(lm(D2TOTDIST15 ~ methcage7,pheno,na.action = na.exclude))
r[r > 850] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST20
# -----------------------
# Locomotor response on Day 2 over 15-20 minute interval.
r <- resid(lm(D2TOTDIST20 ~ methcage7,pheno,na.action = na.exclude))
r[r > 1000] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist20 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST25
# -----------------------
# Locomotor response on Day 2 over 20-25 minute interval.
r <- resid(lm(D2TOTDIST25 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist25 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2TOTDIST30
# -----------------------
# Locomotor response on Day 2 over 25-30 minute interval.
r <- resid(lm(D2TOTDIST30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 900] <- NA
print(histogram(r,col = "firebrick",border = "firebrick",nint = 20,
                xlab = "residual of d2totdist30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST5
# ----------------------
# Locomotor response on Day 3 over 0-5 minute interval.
r <- resid(lm(D3TOTDIST5 ~ methcage7,pheno,na.action = na.exclude))
r[r > 2000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist5 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST10
# -----------------------
# Locomotor response on Day 3 over 5-10 minute interval.
r <- resid(lm(D3TOTDIST10 ~ methcage7,pheno,na.action = na.exclude))
r[r > 4000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist10 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST15
# -----------------------
# Locomotor response on Day 3 over 10-15 minute interval.
r <- resid(lm(D3TOTDIST15 ~ methcage7,pheno,na.action = na.exclude))
r[r > 4000] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST20
# -----------------------
# Locomotor response on Day 3 over 15-20 minute interval.
r <- resid(lm(D3TOTDIST20 ~ methcage7,pheno,na.action = na.exclude))
r[r > 4500] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist20 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST25
# -----------------------
# Locomotor response on Day 3 over 20-25 minute interval.
r <- resid(lm(D3TOTDIST25 ~ methcage7,pheno,na.action = na.exclude))
r[r > 4500] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist25 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3TOTDIST30
# -----------------------
# Locomotor response on Day 3 over 25-30 minute interval.
r <- resid(lm(D3TOTDIST30 ~ methcage7,pheno,na.action = na.exclude))
r[r > 3750] <- NA
print(histogram(r,col = "dodgerblue",border = "dodgerblue",nint = 20,
                xlab = "residual of d3totdist30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1CTRTIME0TO15
# --------------------------
# Proportion of time spent in the center of the arena during Day 1 of
# the MA sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D1ctrtime0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.5)] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 20,
                xlab = "residual of D1ctrtime0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2CTRTIME0TO15
# --------------------------
# Proportion of time spent in the center of the arena during Day 2 of
# the MA sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D2ctrtime0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.75)] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 20,
                xlab = "residual of D2ctrtime0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3CTRTIME0TO15
# --------------------------
# Proportion of time spent in the center of the arena during Day 3 of
# the MA sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D3ctrtime0to15 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.75)] <- NA
print(histogram(r,col = "forestgreen",border = "forestgreen",nint = 20,
                xlab = "residual of D3ctrtime0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1CTRTIME0TO30
# --------------------------
# Proportion of time spent in the center of the arena during Day 1 of
# the MA sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D1ctrtime0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.6)] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D1ctrtime0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2CTRTIME0TO30
# --------------------------
# Proportion of time spent in the center of the arena during Day 2 of
# the MA sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D2ctrtime0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.75)] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D2ctrtime0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3CTRTIME0TO30
# --------------------------
# Proportion of time spent in the center of the arena during Day 3 of
# the MA sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D3ctrtime0to30 ~ methcage7,pheno,na.action = na.exclude))
r[r < (-0.85)] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D3ctrtime0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1HACT0TO15
# -----------------------
# Horizontal activity phenotype measured on Day 1 of the MA
# sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D1hact0to15 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "royalblue",border = "royalblue",nint = 20,
                xlab = "residual of D1hact0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2HACT0TO15
# -----------------------
# Horizontal activity phenotype measured on Day 2 of the MA
# sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D2hact0to15 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "royalblue",border = "royalblue",nint = 20,
                xlab = "residual of D2hact0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3HACT0TO15
# -----------------------
# Horizontal activity phenotype measured on Day 3 of the MA
# sensitivity tests over the 0-15 minute interval.
r <- resid(lm(D3hact0to15 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "royalblue",border = "royalblue",nint = 20,
                xlab = "residual of D3hact0to15 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1HACT0TO30
# -----------------------
# Horizontal activity phenotype measured on Day 1 of the MA
# sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D1hact0to30 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of D1hact0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2HACT0TO30
# -----------------------
# Horizontal activity phenotype measured on Day 2 of the MA
# sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D2hact0to30 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of D2hact0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = DHACT0TO30
# -----------------------
# Horizontal activity phenotype measured on Day 3 of the MA
# sensitivity tests over the 0-30 minute interval.
r <- resid(lm(D3hact0to30 ~ methcage7,pheno,na.action = na.exclude))
print(histogram(r,col = "darkblue",border = "darkblue",nint = 20,
                xlab = "residual of D3hact0to30 | methcage7"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1VACT0TO15
# -----------------------
# Vertical activity phenotype measured on Day 1 of the MA sensitivity
# tests over the 0-15 minute interval.
r <- resid(lm(D1vact0to15 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r < (-0.85) | r > 0.85] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D1vact0to15 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2VACT0TO15
# -----------------------
# Vertical activity phenotype measured on Day 2 of the MA sensitivity
# tests over the 0-15 minute interval.
r <- resid(lm(D2vact0to15 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r < (-1) | r > 1] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D2vact0to15 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3VACT0TO15
# -----------------------
# Vertical activity phenotype measured on Day 3 of the MA sensitivity
# tests over the 0-15 minute interval.
r <- resid(lm(D3vact0to15 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r < (-1.25) | r > 1.25] <- NA
print(histogram(r,col = "limegreen",border = "limegreen",nint = 20,
                xlab = "residual of D3vact0to15 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D1VACT0TO30
# -----------------------
# Vertical activity phenotype measured on Day 1 of the MA sensitivity
# tests over the 0-30 minute interval.
r <- resid(lm(D1vact0to30 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r < (-1) | r > 1] <- NA
print(histogram(r,col = "orangered",border = "orangered",nint = 20,
                xlab = "residual of D1vact0to30 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D2VACT0TO30
# -----------------------
# Vertical activity phenotype measured on Day 2 of the MA sensitivity
# tests over the 0-30 minute interval.
r <- resid(lm(D2vact0to30 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r < (-1) | r > 1] <- NA
print(histogram(r,col = "orangered",border = "orangered",nint = 20,
                xlab = "residual of D2vact0to30 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = D3VACT0TO30
# -----------------------
# Vertical activity phenotype measured on Day 3 of the MA sensitivity
# tests over the 0-30 minute interval.
r <- resid(lm(D3vact0to30 ~ methcage7 + methcage8 + methcage9 + methcage10 +
              methcage11 + methcage12,pheno,na.action = na.exclude))
r[r > 1.5] <- NA
print(histogram(r,col = "orangered",border = "orangered",nint = 20,
                xlab = "residual of D3vact0to30 | methcage"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PP3PPIAVG
# ---------------------
# Prepulse inhibition during the 3 dB trials given the box used for
# the PPI tests.
r <- resid(lm(pp3PPIavg ~ PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r < (-0.9)] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of pp3PPIavg | PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PP6PPIAVG
# ---------------------
# Prepulse inhibition during the 6 dB trials given the box used for
# the PPI tests.
r <- resid(lm(pp6PPIavg ~ PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r < (-1.1)] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of pp6PPIavg | PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PP12PPIAVG
# ----------------------
# Prepulse inhibition during the 6 dB trials given the box used for
# the PPI tests.
r <- resid(lm(pp12PPIavg ~ PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r < (-1)] <- NA
print(histogram(r,col = "darkblue",border = "darkblue",nint = 24,
                xlab = "residual of pp12PPIavg | PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = PPIAVG
# ------------------
# Prepulse inhibition averaged over all the trials given the box used
# for the PPI tests.
r <- resid(lm(PPIavg ~ PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r < (-1)] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "residual of PPIavg | PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = P120B4
# ------------------
# Startle response during the fourth block of pulse-alone trials
# (p120b1) given the startle response during the first block of
# pulse-along trials (p120b4), and given the box used for the PPI
# tests. This is a measure of habituation to the pulses.
r <- resid(lm(p120b4 ~ p120b1 + PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r < (-140) | r > (140)] <- NA
print(histogram(r,col = "orangered",border = "orangered",nint = 24,
                xlab = "residual of p120b4 | p120b1 + PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))
readLines("stdin",n = 1)

# PHENOTYPE = STARTLE
# -------------------
# Average startle response in the pulse-only trials given the box used
# for the PPI tests.
r <- resid(lm(startle ~ PPIbox1 + PPIbox2 + PPIbox3 + PPIbox4,
              pheno,na.action = na.exclude))
r[r > 250] <- NA
print(histogram(r,col = "darkorange",border = "darkorange",nint = 24,
                xlab = "PPI startle | PPIbox"))
cat("num. phenotype samples =",sum(!is.na(r)),"\n")
print(check.normal.quantiles(r))

