# Script for generating the LD decay plot.
library(lattice)
library(latticeExtra)

# Load the LD calculations in the 6 mouse populations.
load("LDdecay.HS.RData")
ld.hs <- ld
load("LDdecay.CFW.RData")
ld.cfw <- ld
load("LDdecay.MDA.RData")
ld.mda <- ld
load("LDdecay.LGSM.RData")
ld.lgsm <- ld
load("LDdecay.HMDP.RData")
ld.hmdp <- ld
load("LDdecay.DO.RData")
ld.do <- ld

# Plot LD (r^2) against base-pair position.
trellis.device(height = 3.5,width = 3.75)
n   <- length(bins)
pos <- bins[1:(n-1)]/1e6
print(xyplot(ld ~ pos,data.frame(pos = pos,ld = colMeans(ld.hs,na.rm = TRUE)),
             type = "l",lwd = 2,col = "dodgerblue",xlab = "position (Mb)",
             ylab = "r^2",
             scales = list(y = list(limits=c(-0.02,0.85),at=seq(0,1,0.1)))) +
    as.layer(xyplot(ld~pos,data.frame(pos=pos,ld=colMeans(ld.mda,na.rm=TRUE)),
                    type = "l",lwd = 2,col = "midnightblue")) +
    as.layer(xyplot(ld~pos,data.frame(pos=pos,ld=colMeans(ld.lgsm,na.rm=TRUE)),
                     type = "l",lwd = 2,col = "forestgreen")) +
    as.layer(xyplot(ld~pos,data.frame(pos=pos,ld=colMeans(ld.hmdp,na.rm=TRUE)),
                    type = "l",lwd = 2,col = "limegreen")) +
    as.layer(xyplot(ld~pos,data.frame(pos=pos,ld=colMeans(ld.do,na.rm=TRUE)),
                    type = "l",lwd = 2,col = "cyan")) +
    as.layer(xyplot(ld~pos,data.frame(pos=pos,ld=colMeans(ld.cfw,na.rm=TRUE)),
                    type = "l",lwd = 2,col = "darkorange")))

