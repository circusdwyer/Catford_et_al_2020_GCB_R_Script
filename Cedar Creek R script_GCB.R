rm(list=ls())

## packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(multcomp)
library(multcompView)
library(factoextra)

## colours for warming treatments
control_col <- "blue3"
low_col <-"gold"
high_col <- "red"


## load custom functions
source("R/JMD_R_functions.r")

## read in and compile the trait data
source("R/Compile_trait_data.R")

## summarise to subplot-scale data from 'species_in_subplot_data'
## also join with trait data
source("R/Compile_subplot_data.R")

## read in and compile trait data for the larger species pool
source("R/Compile_pool_data.R")

##update species names using names provided by Jane Catford
alt_names_exp<-e93_names$USDA[match(species_in_subplot_data$Species, e93_names$CDR)]
species_in_subplot_data$Species<-ifelse(is.na(alt_names_exp), species_in_subplot_data$Species, alt_names_exp)
species_in_subplot_data$Species<-ifelse(species_in_subplot_data$Species=="Dichanthelium villosissimum (Nash) Freckmann var. praecocius (Hitchc. & Chase)", 
                                        "Dichanthelium villosissimum", species_in_subplot_data$Species)
species_in_subplot_data$Species<-ifelse(species_in_subplot_data$Species=="Plantago (purshii) patagonica", 
                                        "Plantago patagonica", species_in_subplot_data$Species)
species_in_subplot_data$Species<-ifelse(species_in_subplot_data$Species=="Tragopogon dubius (major)", 
                                        "Tragopogon dubius", species_in_subplot_data$Species)
species_in_subplot_data$perennial<-ifelse(species_in_subplot_data$Duration == "Perennial", 1,0)

## replace names in traits dataframe too
alt_names_traits<-e93_names$USDA[match(traits$Species, e93_names$CDR)]
traits$Species<-ifelse(is.na(alt_names_traits), traits$Species, alt_names_traits)
traits$Species<-ifelse(traits$Species=="Dichanthelium villosissimum (Nash) Freckmann var. praecocius (Hitchc. & Chase)", 
                            "Dichanthelium villosissimum", traits$Species)
traits$Species<-ifelse(traits$Species=="Plantago (purshii) patagonica", 
                            "Plantago patagonica", traits$Species)
traits$Species<-ifelse(traits$Species=="Tragopogon dubius (major)", 
                            "Tragopogon dubius", traits$Species)


## summarise to plot-scale data from 'species_in_subplot_data' 
## this also creates 'species_in_plot_data'
source("R/Compile_plot_data.R")



################
### FIGURE 2####
################

## first run the models
subplot_data$log.sown.richness<- log(subplot_data$sown.richness)
subplot_data$log.subplot.sown.Mass.g.m2<- log(subplot_data$subplot.sown.Mass.g.m2)

## richness
m.rich.rich1<-lmer(subplot.colonist.richness ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
plot(fitted(m.rich.rich1), resid(m.rich.rich1))
summary(m.rich.rich1)
anova(m.rich.rich1)
## remove interaction
m.rich.rich2<-lmer(subplot.colonist.richness ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.rich.rich2)


## biomass
m.bio.rich1<-lmer(log(subplot.colonist.Mass.g.m2) ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
plot(fitted(m.bio.rich1), resid(m.bio.rich1))
summary(m.bio.rich1)
anova(m.bio.rich1)
m.bio.rich2<-lmer(log(subplot.colonist.Mass.g.m2) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.bio.rich2)
anova(m.bio.rich2)

m.bio.rich3<-lmer(log(subplot.colonist.Mass.g.m2) ~ log.sown.richness + (1|Plot), subplot_data)
summary(m.bio.rich3)


pdf(file=paste("Outputs/Figures/Figure_2.pdf", sep=""),height=5.5,width=10.5,useDingbats=F)   
par(mfrow=c(1,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

plot.x<-seq.func(log(subplot_data$sown.richness))

m.rich.rich.pred.control<-lmer.predict(mod=m.rich.rich2, newdat=data.frame(rep(1,100), plot.x,0 ,0), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.richness ~ jitter(log.sown.richness, amount=0.05), ylab="Colonist richness in subplot", xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=m.rich.rich.pred.control$y, upper=m.rich.rich.pred.control$phi, lower=m.rich.rich.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=1)
m.rich.rich.pred.low<-lmer.predict(mod=m.rich.rich2, newdat=data.frame(rep(1,100), plot.x,0 ,1), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.richness ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=m.rich.rich.pred.low$y, upper=m.rich.rich.pred.low$phi, lower=m.rich.rich.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=1)
m.rich.rich.pred.high<-lmer.predict(mod=m.rich.rich2, newdat=data.frame(rep(1,100), plot.x,1 ,0), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.richness ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.rich.rich.pred.high$y, upper=m.rich.rich.pred.high$phi, lower=m.rich.rich.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=1)
mtext("(a)", line=0.5, side=3, adj=0, cex=1.2, font=2)
text(rep(log(4.5), 3), c(37, 35.25, 33.5), labels=c("Control (no warming)", "Low", "High"), col=c(control_col, low_col, high_col), pos=4)

m.bio.rich.pred<-lmer.predict(mod=m.bio.rich3, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), log="y", ylab=expression(Total~colonist~biomass~'in'~subplot~(g~m^{-2})), 
                                      xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.bio.rich.pred$y, upper=m.bio.rich.pred$phi, lower=m.bio.rich.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(b)", line=0.5, side=3, adj=0, cex=1.2, font=2)
dev.off()



################
### FIGURE 3####
################

### Colonists
pdf(file=paste("Outputs/Figures/Figure_3.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

plot.x<-seq.func(subplot_data$log.sown.richness)

## c4 grass
mc4full<-lmer(log(subplot.colonist.c4.grass.Mass.g.m2+0.01) ~ log.subplot.sown.Mass.g.m2 + (log.sown.richness + I(log.sown.richness^2)) * Heat.treatment + (1|Plot), subplot_data)
summary(mc4full)
mc4b<-lmer(log(subplot.colonist.c4.grass.Mass.g.m2+0.01) ~ log.sown.richness +I(log.sown.richness^2) + Heat.treatment + (1|Plot), subplot_data)
summary(mc4b)
anova(mc4b)
mc4<-lmer(log(subplot.colonist.c4.grass.Mass.g.m2+0.01) ~ log.sown.richness +I(log.sown.richness^2) + (1|Plot), subplot_data)
summary(mc4)

mc4.pred<-lmer.predict(mod=mc4, newdat=data.frame(rep(1,100), plot.x, plot.x^2), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.c4.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonizing~C4~grass~biomass~(g~m^{-2})), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.c4.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.c4.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mc4.pred$y, upper=mc4.pred$phi, lower=mc4.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(a)", line=0.5, side=3, adj=0, cex=1, font=2)


## c3 grass
mc3full<-lmer(log(subplot.colonist.c3.grass.Mass.g.m2 + 0.01) ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment+ (1|Plot), subplot_data)
summary(mc3full)

mc3<-lmer(log(subplot.colonist.c3.grass.Mass.g.m2 + 0.01) ~ log.sown.richness * Heat.treatment+ (1|Plot), subplot_data)


mc3.pred.control<-lmer.predict(mod=mc3, newdat=data.frame(rep(1,100), plot.x,0 ,0,plot.x*0,plot.x*0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.c3.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonizing~C3~grass~biomass~(g~m^{-2})), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(0.01, 50)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=mc3.pred.control$y, upper=mc3.pred.control$phi, lower=mc3.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=1)
mc3.pred.low<-lmer.predict(mod=mc3, newdat=data.frame(rep(1,100), plot.x,0 ,1,plot.x*0,plot.x*1), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.c3.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=mc3.pred.low$y, upper=mc3.pred.low$phi, lower=mc3.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=1)
mc3.pred.high<-lmer.predict(mod=mc3, newdat=data.frame(rep(1,100), plot.x,1 ,0,plot.x*1,plot.x*0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.c3.grass.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mc3.pred.high$y, upper=mc3.pred.high$phi, lower=mc3.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=1)
mtext("(b)", line=0.5, side=3, adj=0, cex=1, font=2)
text(rep(log(4.5), 3), c(50,30,18), labels=c("Control (no warming)", "Low", "High"), col=c(control_col, low_col, high_col), pos=4)


## legumes
mlegumefull<-lmer(log(subplot.colonist.legume.Mass.g.m2 + 0.01) ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment+ (1|Plot), subplot_data)
summary(mlegumefull)
mlegumeb<-lmer(log(subplot.colonist.legume.Mass.g.m2+0.01) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(mlegumeb)
anova(mlegumeb)
mlegume<-lmer(log(subplot.colonist.legume.Mass.g.m2+0.01) ~ log.sown.richness + (1|Plot), subplot_data)
summary(mlegume)

mlegume.pred<-lmer.predict(mod=mlegume, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.legume.Mass.g.m2+0.01 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonizing~legume~biomass~(g~m^{-2})), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(0.01, 20)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.legume.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.legume.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mlegume.pred$y, upper=mlegume.pred$phi, lower=mlegume.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(c)", line=0.5, side=3, adj=0, cex=1, font=2)


## Forbs
mforbfull<-lmer(log(subplot.colonist.forb.Mass.g.m2 + 0.01) ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment+ (1|Plot), subplot_data)
summary(mforbfull)
mforbb<-lmer(log(subplot.colonist.forb.Mass.g.m2+0.01) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(mforbb)
anova(mforbb)
mforb<-lmer(log(subplot.colonist.forb.Mass.g.m2+0.01) ~ log.sown.richness + (1|Plot), subplot_data)
summary(mforb)

mforb.pred<-lmer.predict(mod=mforb, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=T)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonizing~forb~biomass~(g~m^{-2})), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(0.1,200)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mforb.pred$y, upper=mforb.pred$phi, lower=mforb.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(d)", line=0.5, side=3, adj=0, cex=1, font=2)

dev.off()


################
### FIGURE 4####
################

pdf(file=paste("Outputs/Figures/Figure_4.pdf", sep=""),height=9,width=8.5, useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

plot.x<-seq.func(subplot_data$log.sown.richness)

## SLA
m.sla1<-lmer(subplot.colonist.log.sla.cwm ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
summary(m.sla1)
anova(m.sla1)

m.sla2<-lmer(subplot.colonist.log.sla.cwm ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.sla2)

m.sla3<-lmer(subplot.colonist.log.sla.cwm ~ log.sown.richness + (1|Plot), subplot_data)
summary(m.sla3)

m.sla.pred<-exp(lmer.predict(mod=m.sla3, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=F))
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(exp(subplot.colonist.log.sla.cwm) ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonist~SLA~(mm^{2}~mg^{-1})~CWM), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(exp(subplot.colonist.log.sla.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(exp(subplot.colonist.log.sla.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.sla.pred$y, upper=m.sla.pred$phi, lower=m.sla.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(a)", line=0.5, side=3, adj=0, cex=1, font=2)

## seed mass
m.sm1<-lmer(subplot.colonist.log.sm.cwm ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
summary(m.sm1)
anova(m.sm1)
m.sm2<-lmer(subplot.colonist.log.sm.cwm ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.sm2)
anova(m.sm2)

m.sm.pred.control<-exp(lmer.predict(mod=m.sm2, newdat=data.frame(rep(1,100), plot.x,0 ,0), se.mult=1.96, binom=F, poisson=F))
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(exp(subplot.colonist.log.sm.cwm) ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonist~Seed~Mass~(g)~CWM), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.control$y, upper=m.sm.pred.control$phi, lower=m.sm.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=1)
m.sm.pred.low<-exp(lmer.predict(mod=m.sm2, newdat=data.frame(rep(1,100), plot.x,0 ,1), se.mult=1.96, binom=F, poisson=F))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(exp(subplot.colonist.log.sm.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.low$y, upper=m.sm.pred.low$phi, lower=m.sm.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=1)
m.sm.pred.high<-exp(lmer.predict(mod=m.sm2, newdat=data.frame(rep(1,100), plot.x,1 ,0), se.mult=1.96, binom=F, poisson=F))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(exp(subplot.colonist.log.sm.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.high$y, upper=m.sm.pred.high$phi, lower=m.sm.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=1)
mtext("(b)", line=0.5, side=3, adj=0, cex=1, font=2)
text(rep(log(1), 3), c(-4.3, -4.6, -4.9), labels=c("Control (no warming)", "Low", "High"), col=c(control_col, low_col, high_col), pos=4)


## mh
m.mh1<-lmer(subplot.colonist.sqrt.height.cwm ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
summary(m.mh1)
anova(m.mh1)
m.mh2<-lmer(subplot.colonist.sqrt.height.cwm ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.mh2)
m.mh3<-lmer(subplot.colonist.sqrt.height.cwm ~ log.sown.richness + (1|Plot), subplot_data)
summary(m.mh3)

m.mh.pred<-lmer.predict(mod=m.mh3, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.sqrt.height.cwm ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonist~Maximum~Height~(cm)~CWM), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", yaxt="n"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
axis(side=2, at=sqrt(c(9, 16, 25, 36, 47, 64, 81)), labels=c(9, 16, 25, 36, 47, 64, 81))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.sqrt.height.cwm ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.sqrt.height.cwm ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.mh.pred$y, upper=m.mh.pred$phi, lower=m.mh.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=2)
mtext("(c)", line=0.5, side=3, adj=0, cex=1, font=2)

## LDMC
m.ldmc1<-lmer(subplot.colonist.log.ldmc.cwm ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
summary(m.ldmc1)
anova(m.ldmc1)
m.ldmc2<-lmer(subplot.colonist.log.ldmc.cwm ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.ldmc2)
m.ldmc3<-lmer(subplot.colonist.log.ldmc.cwm ~ log.sown.richness + (1|Plot), subplot_data)
summary(m.ldmc3)

m.ldmc.pred<-exp(lmer.predict(mod=m.ldmc3, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=F))
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(exp(subplot.colonist.log.ldmc.cwm) ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Colonist~LDMC~(mg~g^{-1})~CWM), 
                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(exp(subplot.colonist.log.ldmc.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data[subplot_data$Heat.treatment=="High",], points(exp(subplot.colonist.log.ldmc.cwm) ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.ldmc.pred$y, upper=m.ldmc.pred$phi, lower=m.ldmc.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=2)
mtext("(d)", line=0.5, side=3, adj=0, cex=1, font=2)

dev.off()



## does the seed mass result hold if we remove legumes?
m.sm1.no.legumes<-lmer(subplot.colonist.log.sm.cwm.no.legumes ~ log.subplot.sown.Mass.g.m2 + log.sown.richness * Heat.treatment + (1|Plot), subplot_data)
summary(m.sm1.no.legumes)
anova(m.sm1.no.legumes)
m.sm2.no.legumes<-lmer(subplot.colonist.log.sm.cwm.no.legumes ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data)
summary(m.sm2.no.legumes)
anova(m.sm2.no.legumes)

m.sm.pred.control<-lmer.predict(mod=m.sm2.no.legumes, newdat=data.frame(rep(1,100), plot.x,0 ,0), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="Control",], plot(subplot.colonist.log.sm.cwm.no.legumes ~ jitter(log.sown.richness, amount=0.05), ylab="Colonist log(Seed Mass) CWM (Legumes excluded)", xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.control$y, upper=m.sm.pred.control$phi, lower=m.sm.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=1)
m.sm.pred.low<-lmer.predict(mod=m.sm2.no.legumes, newdat=data.frame(rep(1,100), plot.x,0 ,1), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="Low",], points(subplot.colonist.log.sm.cwm.no.legumes ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.low$y, upper=m.sm.pred.low$phi, lower=m.sm.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=1)
m.sm.pred.high<-lmer.predict(mod=m.sm2.no.legumes, newdat=data.frame(rep(1,100), plot.x,1 ,0), se.mult=1.96, binom=F, poisson=F)
with(subplot_data[subplot_data$Heat.treatment=="High",], points(subplot.colonist.log.sm.cwm.no.legumes ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=m.sm.pred.high$y, upper=m.sm.pred.high$phi, lower=m.sm.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=1)
text(rep(log(1), 3), c(0.1, -0.2, -0.5), labels=c("Control (no warming)", "Low", "High"), col=c(control_col, low_col, high_col), pos=4)




## export fully compiled table (Table S3)
table_apply_function<-function (mod_tab){
  vars<-rownames(mod_tab)
  rownames(mod_tab)<-NULL
  step.1<-round(mod_tab, 3)[,1:2]
  df<-paste(round(mod_tab, 0)[,3], round(mod_tab, 0)[,4], sep=";")
  step.3<-round(mod_tab, 3)[,5:6]
  return(cbind(vars,step.1, df, step.3))
}

Table_S3_list<-list(table_apply_function(anova(m.rich.rich1)), table_apply_function(anova(m.bio.rich1)),
                    table_apply_function(anova(mc4full)), table_apply_function(anova(mc3full)),
                    table_apply_function(anova(mlegumefull)), table_apply_function(anova(mforbfull)),
                    table_apply_function(anova(m.sla1)), table_apply_function(anova(m.sm1)),
                    table_apply_function(anova(m.mh1)), table_apply_function(anova(m.ldmc1)))

names(Table_S3_list) <- c("Richness of invading assemblage",
                          "Biomass of invading assemblage",
                          "log(Biomass of C4 grass colonists",
                          "log(Biomass of C3 grass colonists",
                          "log(Biomass of legume grass colonists",
                          "log(Biomass of forb colonists",
                          "log(SLA) CWM of invading assemblage",
                          "log(Seed Mass) CWM of invading assemblage",
                          "sqrt(Max. Height) CWM of invading assemblage",
                          "log(LDMC) CWM of invading assemblage")

write.csv(do.call('rbind', Table_S3_list), "Outputs/Supps/Table_S3_list.csv")




################
### FIGURE 5####
################

pdf(file=paste("Outputs/Figures/Figure_5.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

plot.x<-seq.func(log(subplot_data$sown.richness))

## c4 grass
subplot_data.sown.c4<-subset(subplot_data, subplot.sown.c4.grass.Mass.g.m2>0)
mc4full<-lmer(log(subplot.sown.c4.grass.Mass.g.m2) ~ log.sown.richness * Heat.treatment + (1|Plot), subplot_data.sown.c4)

mc4<-lmer(log(subplot.sown.c4.grass.Mass.g.m2) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data.sown.c4)
summary(mc4)
anova(mc4)

mc4.pred.control<-lmer.predict(mod=mc4, newdat=data.frame(rep(1,100), plot.x,0 ,0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.c4[subplot_data.sown.c4$Heat.treatment=="Control",], plot(subplot.sown.c4.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Sown~C4~grass~biomass~(g~m^{-2})), 
                                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y"))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=mc4.pred.control$y, upper=mc4.pred.control$phi, lower=mc4.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=1)
mc4.pred.low<-lmer.predict(mod=mc4, newdat=data.frame(rep(1,100), plot.x,0 ,1), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.c4[subplot_data.sown.c4$Heat.treatment=="Low",], points(subplot.sown.c4.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=mc4.pred.low$y, upper=mc4.pred.low$phi, lower=mc4.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=1)
mc4.pred.high<-lmer.predict(mod=mc4, newdat=data.frame(rep(1,100), plot.x,1 ,0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.c4[subplot_data.sown.c4$Heat.treatment=="High",], points(subplot.sown.c4.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mc4.pred.high$y, upper=mc4.pred.high$phi, lower=mc4.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=1)
mtext("(a)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      Positive warming effect (P<0.001)", line=0.5, side=3, adj=0, cex=1)

text(rep(log(1), 3), c(300,245,200), labels=c("Control (no warming)", "Low", "High"), col=c(control_col, low_col, high_col), pos=4)


## c3 grass
subplot_data.sown.c3<-subset(subplot_data, subplot.sown.c3.grass.Mass.g.m2>1)
mc3full<-lmer(log(subplot.sown.c3.grass.Mass.g.m2) ~ log.sown.richness * Heat.treatment + (1|Plot), subplot_data.sown.c3)
mc3b<-lmer(log(subplot.sown.c3.grass.Mass.g.m2) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data.sown.c3)
summary(mc3b)
anova(mc3b)
mc3<-lmer(log(subplot.sown.c3.grass.Mass.g.m2) ~ log.sown.richness + (1|Plot), subplot_data.sown.c3)
summary(mc3)

mc3.pred<-lmer.predict(mod=mc3, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.c3[subplot_data.sown.c3$Heat.treatment=="Control",], plot(subplot.sown.c3.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Sown~C3~grass~biomass~(g~m^{-2})), 
                                                                                 xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(2,120)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data.sown.c3[subplot_data.sown.c3$Heat.treatment=="Low",], points(subplot.sown.c3.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data.sown.c3[subplot_data.sown.c3$Heat.treatment=="High",], points(subplot.sown.c3.grass.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mc3.pred$y, upper=mc3.pred$phi, lower=mc3.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=2)
mtext("(b)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      No warming effect (P=0.07)", line=0.5, side=3, adj=0, cex=1)

## legumes
subplot_data.sown.legume<-subset(subplot_data, subplot.sown.legume.Mass.g.m2>0)
mlegumefull<-lmer(log(subplot.sown.legume.Mass.g.m2) ~ log.sown.richness * Heat.treatment + (1|Plot), subplot_data.sown.legume)
mlegume<-lmer(log(subplot.sown.legume.Mass.g.m2) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data.sown.legume)
summary(mlegume)
anova(mlegume)

mlegume.pred.control<-lmer.predict(mod=mlegume, newdat=data.frame(rep(1,100), plot.x,0 ,0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.legume[subplot_data.sown.legume$Heat.treatment=="Control",], plot(subplot.sown.legume.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Sown~legume~biomass~(g~m^{-2})), 
                                                                                         xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(20,700)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
plot.CI.func(x.for.plot=plot.x, pred=mlegume.pred.control$y, upper=mlegume.pred.control$phi, lower=mlegume.pred.control$plo, env.colour=control_col, env.trans=40, line.colour=control_col, line.weight=3, line.type=2)
mlegume.pred.low<-lmer.predict(mod=mlegume, newdat=data.frame(rep(1,100), plot.x,0 ,1), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.legume[subplot_data.sown.legume$Heat.treatment=="Low",], points(subplot.sown.legume.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
plot.CI.func(x.for.plot=plot.x, pred=mlegume.pred.low$y, upper=mlegume.pred.low$phi, lower=mlegume.pred.low$plo, env.colour=low_col, env.trans=40, line.colour=low_col, line.weight=3, line.type=2)
mlegume.pred.high<-lmer.predict(mod=mlegume, newdat=data.frame(rep(1,100), plot.x,1 ,0), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.legume[subplot_data.sown.legume$Heat.treatment=="High",], points(subplot.sown.legume.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mlegume.pred.high$y, upper=mlegume.pred.high$phi, lower=mlegume.pred.high$plo, env.colour=high_col, env.trans=40, line.colour=high_col, line.weight=3, line.type=2)
mtext("(c)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      Positive warming effect (P=0.004)", line=0.5, side=3, adj=0, cex=1)

## Forbs
subplot_data.sown.forb<-subset(subplot_data, subplot.sown.forb.Mass.g.m2>0)
mforbfull<-lmer(log(subplot.sown.forb.Mass.g.m2) ~ log.sown.richness * Heat.treatment + (1|Plot), subplot_data.sown.forb)
mforbb<-lmer(log(subplot.sown.forb.Mass.g.m2) ~ log.sown.richness + Heat.treatment + (1|Plot), subplot_data.sown.forb)
summary(mforbb)
anova(mforbb)

mforb<-lmer(log(subplot.sown.forb.Mass.g.m2) ~ log.sown.richness + (1|Plot), subplot_data.sown.forb)
summary(mforb)

mforb.pred<-lmer.predict(mod=mforb, newdat=data.frame(rep(1,100), plot.x), se.mult=1.96, binom=F, poisson=T)
with(subplot_data.sown.forb[subplot_data.sown.forb$Heat.treatment=="Control",], plot(subplot.sown.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), ylab=expression(Sown~forb~biomass~(g~m^{-2})), 
                                                                                     xlab="Number of sown species (per plot)", pch=19, col=control_col, xaxt="n", log="y", ylim=c(5,1000)))
axis(side=1, at=log(c(1,2,4,16)), labels=c(1,2,4,16))
with(subplot_data.sown.forb[subplot_data.sown.forb$Heat.treatment=="Low",], points(subplot.sown.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=low_col))
with(subplot_data.sown.forb[subplot_data.sown.forb$Heat.treatment=="High",], points(subplot.sown.forb.Mass.g.m2 ~ jitter(log.sown.richness, amount=0.05), pch=19, col=high_col))
plot.CI.func(x.for.plot=plot.x, pred=mforb.pred$y, upper=mforb.pred$phi, lower=mforb.pred$plo, env.colour="black", env.trans=40, line.colour="black", line.weight=3, line.type=1)
mtext("(d)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      No warming effect (P=0.13)", line=0.5, side=3, adj=0, cex=1)
dev.off()



################
### FIGURE 6####
################

## create relative trait distance variables (and standardise them)
## first at subplot scale
species_in_subplot_data<-left_join(species_in_subplot_data, subplot_no_groups, by="subplot")
species_in_subplot_data$std.subplot.log.sla.rel.dist<-scale(species_in_subplot_data$log.sla - species_in_subplot_data$subplot.sown.log.sla.cwm)[,1]
species_in_subplot_data$std.subplot.log.sm.rel.dist<-scale(species_in_subplot_data$log.sm - species_in_subplot_data$subplot.sown.log.sm.cwm)[,1]
species_in_subplot_data$std.subplot.sqrt.height.rel.dist<-scale(species_in_subplot_data$sqrt.height - species_in_subplot_data$subplot.sown.sqrt.height.cwm)[,1]
species_in_subplot_data$std.subplot.log.ldmc.rel.dist<-scale(species_in_subplot_data$log.ldmc - species_in_subplot_data$subplot.sown.log.ldmc.cwm)[,1]

## plot scale
species_in_plot_data$plot.std.log.sla.rel.dist<-scale(species_in_plot_data$log.sla - species_in_plot_data$plot.sown.log.sla.cwm)[,1]
species_in_plot_data$plot.std.log.sm.rel.dist<-scale(species_in_plot_data$log.sm - species_in_plot_data$plot.sown.log.sm.cwm)[,1]
species_in_plot_data$plot.std.sqrt.height.rel.dist<-scale(species_in_plot_data$sqrt.height - species_in_plot_data$plot.sown.sqrt.height.cwm)[,1]
species_in_plot_data$plot.std.log.ldmc.rel.dist<-scale(species_in_plot_data$log.ldmc - species_in_plot_data$plot.sown.log.ldmc.cwm)[,1]

species_in_subplot_data_colonists <-droplevels(subset(species_in_subplot_data, Taxa_status=="colonist"))
species_in_subplot_data_colonists$std.log.sown.richness <- scale(log(species_in_subplot_data_colonists$sown.richness))[,1]


mEXP1<-lmer(log(Mass.g.m2) ~ (std.log.sla + std.log.sm + std.sqrt.height + std.log.ldmc) * (std.log.sown.richness + Heat.treatment) + (1|Plot) + (1|subplot:Plot) + (std.log.sown.richness + Heat.treatment|Species), data = species_in_subplot_data_colonists, REML=F)
plot(mEXP1)
summary(mEXP1)$coef
write.csv(round(summary(mEXP1)$coef, 3), "Outputs/Supps/Table_S2_unordered.csv")

mEXP2<-lmer(log(Mass.g.m2) ~ (std.log.sla + std.log.sm + std.sqrt.height + std.log.ldmc) * std.log.sown.richness + (1|Plot) + (1|subplot:Plot) + (std.log.sown.richness|Species), data = species_in_subplot_data_colonists, REML=F)
summary(mEXP2)


species_in_plot_data_colonists <-droplevels(subset(species_in_plot_data, Taxa_status=="colonist"))

## standardise the covariates
species_in_plot_data_colonists$plot.std.log.sown.richness<-scale(log(species_in_plot_data_colonists$sown.richness))[,1]
species_in_plot_data_colonists$plot.std.log.sown.Mass.g.m2<-scale(log(species_in_plot_data_colonists$plot.sown.Mass.g.m2))[,1]

traits.model<-lmer(log(Mass.g.m2) ~ (std.log.sla + std.log.sm + std.sqrt.height + std.log.ldmc) * plot.std.log.sown.richness  + (1|Plot) + (plot.std.log.sown.richness|Species), data = species_in_plot_data_colonists, REML=T)
summary(traits.model)
write.csv(round(summary(traits.model)$coef, 3), "Outputs/Supps/Table_S4_model_1.csv")
r.squaredGLMM(traits.model)

## compare raw traits and trait distances using AIC
rel.dist.model<-lmer(log(Mass.g.m2) ~ (plot.std.log.sla.rel.dist + plot.std.log.sm.rel.dist + plot.std.sqrt.height.rel.dist + plot.std.log.ldmc.rel.dist) * plot.std.log.sown.richness  + (1|Plot) + (plot.std.log.sown.richness|Species), data = species_in_plot_data_colonists, REML=T)
summary(rel.dist.model)
write.csv(round(summary(rel.dist.model)$coef, 3), "Outputs/Supps/Table_S4_model_2.csv")
r.squaredGLMM(rel.dist.model)

AICc(rel.dist.model, traits.model)


## do some plots
## this bit grabs the random slopes for each species as well as their associated SEs
varfix <- vcov(traits.model)[2,2]
re <- ranef(traits.model,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm


traits_for_model<-traits[match(rownames(coef(traits.model)$Species), traits$Species),c(6:9)]
traits_for_model$slopes<- rowSums(cbind(1, traits_for_model) * as.matrix(coef(traits.model)$Species[,c(6:10)]))
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
traits_for_model$Species<-rownames(coef(traits.model)$Species)

## species to highlight in the plots
highlight_species<-c("Andropogon gerardii", "Baptisia alba", "Berteroa incana", "Chamaesyce glyptosperma", "Crepis tectorum", 
"Dalea purpurea", "Lupinus perennis", "Mollugo verticillata", "Physalis hispida", "Polygonum tenue", 
"Rosa arkansana", "Rumex acetosella", "Silene antirrhina", "Sorghastrum nutans")

high_div_colonist_col <-  "blue4" 
low_div_colonist_col <- "red"
  
pdf(file=paste("Outputs/Figures/Figure_6.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 3, 2), mgp=c(2.6,0.8,0)) 

## SLA
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(species_in_plot_data_colonists$std.log.sla))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(traits.model)[c(6,7)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(traits.model)[c(6,7), c(6,7)], as.matrix(data.for.sla.slope.regression))))
with(traits_for_model, plot(slopes ~ std.log.sla, ylab= "Slope of colonist biomass - sown community richness relationship", xlab="log(SLA) (standardised)", ylim=c(-2.2,2), type="n"))
plot.CI.func(x.for.plot=seq.func(species_in_plot_data_colonists$std.log.sla), pred=pred.sla.slope.regression, 
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs, 
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
with(traits_for_model, points(std.log.sla, slopes, col="grey70", pch=19))
with(traits_for_model, arrows(std.log.sla, slopes+slope.SEs, std.log.sla, slopes-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], points(std.log.sla, slopes, 
            col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), pch=19))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], arrows(std.log.sla, slopes+slope.SEs, std.log.sla, slopes-slope.SEs, length = 0, 
            angle = 30, code = 2, col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], text(std.log.sla, slopes, labels=highlight_species, cex=0.5, font=3))
curve(cbind(1,x)%*%fixef(traits.model)[c(6,7)], add=T, lwd=3) 
abline(h=0, lty=3, lwd=2)
mtext("(a)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      log(SLA)*", line=0.5, side=3, adj=0, cex=1)

## seed mass
data.for.sm.slope.regression<-cbind(int=rep(1,100), x=seq.func(species_in_plot_data_colonists$std.log.sm))
pred.sm.slope.regression<- data.for.sm.slope.regression%*%fixef(traits.model)[c(6,8)]
sm.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sm.slope.regression) %*% tcrossprod(vcov(traits.model)[c(6,8), c(6,8)], as.matrix(data.for.sm.slope.regression))))
with(traits_for_model, plot(slopes ~ std.log.sm, ylab= "Slope of colonist biomass - sown community richness relationship", xlab="log(Seed mass) (standardised)", ylim=c(-2.2,2), type="n"))
plot.CI.func(x.for.plot=seq.func(species_in_plot_data_colonists$std.log.sm), pred=pred.sm.slope.regression, 
             upper=pred.sm.slope.regression + 1.96*sm.slope.regression.SEs,
             lower=pred.sm.slope.regression - 1.96*sm.slope.regression.SEs, 
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
with(traits_for_model, points(std.log.sm, slopes, col="grey70", pch=19))
with(traits_for_model, arrows(std.log.sm, slopes+slope.SEs, std.log.sm, slopes-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], points(std.log.sm, slopes, 
                                                                               col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), pch=19))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], arrows(std.log.sm, slopes+slope.SEs, std.log.sm, slopes-slope.SEs, length = 0, 
                                                                               angle = 30, code = 2, col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], text(std.log.sm, slopes, labels=highlight_species, cex=0.5, font=3))
curve(cbind(1,x)%*%fixef(traits.model)[c(6,8)], add=T, lwd=3) 
abline(h=0, lty=3, lwd=2)
mtext("(b)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      log(Seed mass)*", line=0.5, side=3, adj=0, cex=1)

## height
data.for.mh.slope.regression<-cbind(int=rep(1,100), x=seq.func(species_in_plot_data_colonists$std.sqrt.height))
pred.mh.slope.regression<- data.for.mh.slope.regression%*%fixef(traits.model)[c(6,9)]
sm.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.mh.slope.regression) %*% tcrossprod(vcov(traits.model)[c(6,9), c(6,9)], as.matrix(data.for.mh.slope.regression))))
with(traits_for_model, plot(slopes ~ std.sqrt.height, ylab= "Slope of colonist biomass - sown community richness relationship", xlab="sqrt(Max. height) (standardised)", ylim=c(-2.2,2), type="n"))
plot.CI.func(x.for.plot=seq.func(species_in_plot_data_colonists$std.sqrt.height), pred=pred.mh.slope.regression, 
             upper=pred.mh.slope.regression + 1.96*sm.slope.regression.SEs,
             lower=pred.mh.slope.regression - 1.96*sm.slope.regression.SEs, 
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
with(traits_for_model, points(std.sqrt.height, slopes, col="grey70", pch=19))
with(traits_for_model, arrows(std.sqrt.height, slopes+slope.SEs, std.sqrt.height, slopes-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], points(std.sqrt.height, slopes, 
                                                                               col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), pch=19))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], arrows(std.sqrt.height, slopes+slope.SEs, std.sqrt.height, slopes-slope.SEs, length = 0, 
                                                                               angle = 30, code = 2, col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], text(std.sqrt.height, slopes, labels=highlight_species, cex=0.5, font=3))
curve(cbind(1,x)%*%fixef(traits.model)[c(6,9)], add=T, lwd=3) 
abline(h=0, lty=3, lwd=2)
mtext("(c)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      sqrt(Max. height)*", line=0.5, side=3, adj=0, cex=1)

## LDMC
data.for.ldmc.slope.regression<-cbind(int=rep(1,100), x=seq.func(species_in_plot_data_colonists$std.log.ldmc))
pred.ldmc.slope.regression<- data.for.ldmc.slope.regression%*%fixef(traits.model)[c(6,10)]
sm.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.ldmc.slope.regression) %*% tcrossprod(vcov(traits.model)[c(6,10), c(6,10)], as.matrix(data.for.ldmc.slope.regression))))
with(traits_for_model, plot(slopes ~ std.log.ldmc, ylab= "Slope of colonist biomass - sown community richness relationship", xlab="log(LDMC) (standardised)", ylim=c(-2.2,2), type="n"))
plot.CI.func(x.for.plot=seq.func(species_in_plot_data_colonists$std.log.ldmc), pred=pred.ldmc.slope.regression, 
             upper=pred.ldmc.slope.regression + 1.96*sm.slope.regression.SEs,
             lower=pred.ldmc.slope.regression - 1.96*sm.slope.regression.SEs, 
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=2)
with(traits_for_model, points(std.log.ldmc, slopes, col="grey70", pch=19))
with(traits_for_model, arrows(std.log.ldmc, slopes+slope.SEs, std.log.ldmc, slopes-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], points(std.log.ldmc, slopes, 
                                                                               col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), pch=19))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], arrows(std.log.ldmc, slopes+slope.SEs, std.log.ldmc, slopes-slope.SEs, length = 0, 
                                                                               angle = 30, code = 2, col=ifelse(slopes>0, high_div_colonist_col, low_div_colonist_col), lwd = 1))
with(traits_for_model[traits_for_model$Species %in%highlight_species,], text(std.log.ldmc, slopes, labels=highlight_species, cex=0.5, font=3))
curve(cbind(1,x)%*%fixef(traits.model)[c(6,10)], add=T, lwd=3, lty=2) 
abline(h=0, lty=3, lwd=2)
mtext("(d)", line=0.5, side=3, adj=0, cex=1, font = 2)
mtext("      log(LDMC) N.S.", line=0.5, side=3, adj=0, cex=1)

dev.off()




################
### FIGURE S1###
################

## how do continuous traits vary by functional group?
traits_for_plot<-traits %>%
  filter(tilman.groups %in% c("C4 Grass", "C3 Grass", "Forb", "Legume"))

pdf(file=paste("Outputs/Supps/Figure_S1_traits_by_Functional_group.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

plot(traits_for_plot$log.sla ~ as.factor(traits_for_plot$tilman.groups), ylab="log(SLA)", xlab="Functional group")
mtext("(a)", line=0.5, side=3, adj=0, cex=1)
plot(traits_for_plot$log.sm ~ as.factor(traits_for_plot$tilman.groups), ylab="log(Seed mass)", xlab="Functional group")
mtext("(b)", line=0.5, side=3, adj=0, cex=1)
plot(traits_for_plot$sqrt.height ~ as.factor(traits_for_plot$tilman.groups), ylab="sqrt(Max. height)", xlab="Functional group")
mtext("(c)", line=0.5, side=3, adj=0, cex=1)
plot(traits_for_plot$log.ldmc ~ as.factor(traits_for_plot$tilman.groups), ylab="log(LDMC)", xlab="Functional group")
mtext("(d)", line=0.5, side=3, adj=0, cex=1)

dev.off()

################
### FIGURE S2###
################

fg.by.lifespan<-as.factor(paste(traits_for_plot$Duration, traits_for_plot$tilman.groups, sep=" "))

pdf(file=paste("Outputs/Supps/Figure_S2_traits_by_FG_&_lifespan.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

duration_labels1<-c("Annual", "Annual", "Perennial", "Perennial")
group_labels1<-c("C4 Grass", "Legume", "C3 Grass", "Forb")
duration_labels2<-c("Annual", "Biennial", "Perennial", "Perennial")
group_labels2<-c("Forb", "Forb", "C4 Grass", "Legume")


plot(traits_for_plot$log.sla ~ fg.by.lifespan, ylab="log(SLA)", xlab="", xaxt="n")
axis(side=1, at=c(1,3,5,7), labels = rep("", 4), tck= -0.03)
axis(side=1, at=c(2,4,6,8), labels = rep("", 4), tck= -0.12)
mtext(at=c(1,3,5,7), text=duration_labels1, side=1, line=0.75, cex=0.75)
mtext(at=c(1,3,5,7), text=group_labels1, side=1, line=1.5, cex=0.75)
mtext(at=c(2,4,6,8), text=duration_labels2, side=1, line=2.5, cex=0.75)
mtext(at=c(2,4,6,8), text=group_labels2, side=1, line=3.25, cex=0.75)
mtext("(a)", line=0.5, side=3, adj=0, cex=1)
      
plot(traits_for_plot$log.sm ~ fg.by.lifespan, ylab="log(Seed mass)", xlab="", xaxt="n")
axis(side=1, at=c(1,3,5,7), labels = rep("", 4), tck= -0.03)
axis(side=1, at=c(2,4,6,8), labels = rep("", 4), tck= -0.12)
mtext(at=c(1,3,5,7), text=duration_labels1, side=1, line=0.75, cex=0.75)
mtext(at=c(1,3,5,7), text=group_labels1, side=1, line=1.5, cex=0.75)
mtext(at=c(2,4,6,8), text=duration_labels2, side=1, line=2.5, cex=0.75)
mtext(at=c(2,4,6,8), text=group_labels2, side=1, line=3.25, cex=0.75)
mtext("(b)", line=0.5, side=3, adj=0, cex=1)

plot(traits_for_plot$sqrt.height ~ fg.by.lifespan, ylab="sqrt(Max. height)", xlab="", xaxt="n")
axis(side=1, at=c(1,3,5,7), labels = rep("", 4), tck= -0.03)
axis(side=1, at=c(2,4,6,8), labels = rep("", 4), tck= -0.12)
mtext(at=c(1,3,5,7), text=duration_labels1, side=1, line=0.75, cex=0.75)
mtext(at=c(1,3,5,7), text=group_labels1, side=1, line=1.5, cex=0.75)
mtext(at=c(2,4,6,8), text=duration_labels2, side=1, line=2.5, cex=0.75)
mtext(at=c(2,4,6,8), text=group_labels2, side=1, line=3.25, cex=0.75)
mtext("(c)", line=0.5, side=3, adj=0, cex=1)

plot(traits_for_plot$log.ldmc ~ fg.by.lifespan, ylab="log(LDMC)", xlab="", xaxt="n")
axis(side=1, at=c(1,3,5,7), labels = rep("", 4), tck= -0.03)
axis(side=1, at=c(2,4,6,8), labels = rep("", 4), tck= -0.12)
mtext(at=c(1,3,5,7), text=duration_labels1, side=1, line=0.75, cex=0.75)
mtext(at=c(1,3,5,7), text=group_labels1, side=1, line=1.5, cex=0.75)
mtext(at=c(2,4,6,8), text=duration_labels2, side=1, line=2.5, cex=0.75)
mtext(at=c(2,4,6,8), text=group_labels2, side=1, line=3.25, cex=0.75)
mtext("(d)", line=0.5, side=3, adj=0, cex=1)

dev.off()



#############################################
####### FIGURE S3 - SPECIES POOL PCA ########
#############################################
## run the PCA
pool_trait_pca<-princomp(pool_traits[,2:5], cor=T)
summary(pool_trait_pca) ##contribution of each PC axis
loadings(pool_trait_pca) ## trait loadings on each PC axis

## make a nice version - colour-code species pool groups
colonist_col <- "red"
sown_col <- "blue"
other_col <- "grey"

pdf(file=paste0("Outputs/Supps/Figure_S3.pdf"), height=10, width=10,useDingbats=F)   
print(fviz_pca_biplot(pool_trait_pca, axes = c(1, 2), col.var="black",geom = c("point", "text"), 
                      habillage = pool_traits$pool_category,
                      title = "PC1 & PC2",
                      palette = c(colonist_col, other_col, sown_col)))
dev.off()


##################################################
####### FIGURE S4 - SPECIES POOL BW PLOTS ########
##################################################
pool_traits$pool_cat_intro<-as.factor(paste(pool_traits$pool_category, ifelse(pool_traits$introduced==1, "intro", "native"), sep="_"))           
x_axis_labels_top<-c("Non-native", "Native", "Non-native", "Native", "Sown")
x_axis_labels_bottom<-c("colonist", "colonist", "non-colonist", "non-colonist", "")

pdf(file=paste0("Outputs/Supps/Figure_S4.pdf"), height=11, width=11, useDingbats=F)
par(mfrow=c(2,2))
with(pool_traits, plot(exp(log.sla) ~ pool_cat_intro, log="y", ylab="SLA", xlab="", xaxt="n", ylim=c(5,75)))
axis(side=1, at=1:5, x_axis_labels_top, cex.axis=1)
mtext(at=1:5, text=x_axis_labels_bottom, side=1, line=2, cex=0.8)
msla<-glht(aov(log.sla ~ pool_cat_intro, pool_traits), mcp(pool_cat_intro = "Tukey"))
sla_letters<-cld(msla,level=0.05,decreasing=TRUE)$mcletters$Letters
text(1:5, rep(75,5), sla_letters, font=2, cex=1.5)
mtext("(a)", line=0.5, side=3, adj=0, cex=1.5)

with(pool_traits, plot(exp(log.sm) ~ pool_cat_intro, log="y", ylab="Seed mass (g)", xlab="", xaxt="n", ylim=c(0.00001,0.1)))
axis(side=1, at=1:5, x_axis_labels_top, cex.axis=1)
mtext(at=1:5, text=x_axis_labels_bottom, side=1, line=2, cex=0.8)
msm<-glht(aov(log.sm ~ pool_cat_intro, pool_traits), mcp(pool_cat_intro = "Tukey"))
sm_letters<-cld(msm,level=0.05,decreasing=TRUE)$mcletters$Letters
text(1:5, rep(0.1,5), sm_letters, font=2, cex=1.5)
mtext("(b)", line=0.5, side=3, adj=0, cex=1.5)

with(pool_traits, plot(sqrt.height ~ pool_cat_intro, ylab="Max. height (cm)", xlab="", xaxt="n", yaxt="n", ylim=c(1,13)))
axis(side=1, at=1:5, x_axis_labels_top, cex.axis=1)
mtext(at=1:5, text=x_axis_labels_bottom, side=1, line=2, cex=0.8)
axis(side=2, c(2,4,6,8,10,12), labels=c(2,4,6,8,10,12)^2)
mmh<-glht(aov(sqrt.height ~ pool_cat_intro, pool_traits), mcp(pool_cat_intro = "Tukey"))
mh_letters<-cld(mmh,level=0.05,decreasing=TRUE)$mcletters$Letters
text(1:5, rep(13,5), mh_letters, font=2, cex=1.5)
mtext("(c)", line=0.5, side=3, adj=0, cex=1.5)

with(pool_traits, plot(exp(log.ldmc) ~ pool_cat_intro, log="y", ylab="LDMC", xlab="", xaxt="n", ylim=c(100,700)))
axis(side=1, at=1:5, x_axis_labels_top, cex.axis=1)
mtext(at=1:5, text=x_axis_labels_bottom, side=1, line=2, cex=0.8)
mldmc<-glht(aov(log.ldmc ~ pool_cat_intro, pool_traits), mcp(pool_cat_intro = "Tukey"))
ldmc_letters<-cld(mldmc,level=0.05,decreasing=TRUE)$mcletters$Letters
text(1:5, rep(700,5), ldmc_letters, font=2, cex=1.5)
mtext("(d)", line=0.5, side=3, adj=0, cex=1.5)

dev.off()


################
### FIGURE S5###
################

### DO sown CWMS CHANGE WITH WARMING? ###
## Proper multiple comparisons to test if sown CWMs differ between warming treatments (within each diversity level)
## make a contrast matrix indicating which combinations we want to compare
K<-matrix(c(1,0,0,-1,0,0,0,0,0,
            0,1,0,0,-1,0,0,0,0,
            0,0,1,0,0,-1,0,0,0,
            1,0,0,0,0,0,-1,0,0,
            0,1,0,0,0,0,0,-1,0,
            0,0,1,0,0,0,0,0,-1,
            0,0,0,1,0,0,-1,0,0,
            0,0,0,0,1,0,0,-1,0,
            0,0,0,0,0,1,0,0,-1), byrow=T, nrow=9, ncol=9,
          dimnames=list(paste(rep(c(1,4,16), 3), rep(c("Control - Low", "Control - High", "Low - High"), each=3), sep=":"), paste(rep(c(1,4,16), 3), rep(c("Control", "low", "high"), each=3), sep=":")))

## make an interaction term and associated 'cell' model so the glht gives us the contrasts we want
subplot_data$div.by.warming <- with(subplot_data, interaction(as.factor(subplot.original.NumSp), Heat.treatment))

## SLA
cell.sla.m <- lmer(subplot.sown.log.sla.cwm~-1+div.by.warming+(1|Plot),data=subplot_data)
summary(glht(cell.sla.m, linfct = K)) 

## Seed Mass
cell.sm.m <- lmer(subplot.sown.log.sm.cwm~-1+div.by.warming+(1|Plot),data=subplot_data)
summary(glht(cell.sm.m, linfct = K)) 

## Max. Height
cell.mh.m <- lmer(subplot.sown.sqrt.height.cwm~-1+div.by.warming+(1|Plot),data=subplot_data)
summary(glht(cell.mh.m, linfct = K)) 

##LDMC
cell.ldmc.m <- lmer(subplot.sown.log.ldmc.cwm~-1+div.by.warming+(1|Plot),data=subplot_data)
summary(glht(cell.ldmc.m, linfct = K)) 


## now make the plots
## create a warming treatment that will order as we want it on the plot
subplot_data$plot_heat<-ifelse(subplot_data$Heat.treatment=="Control", 0, ifelse(subplot_data$Heat.treatment=="Low", 1, 2))
## same for sown diversity
subplot_data$plot_div<-ifelse(subplot_data$subplot.original.NumSp==1, 0, ifelse(subplot_data$subplot.original.NumSp==4, 1, 2))
## combine into a unique identifier for the diversity x heat combinations to make some boxplots
subplot_data$div_heat<-as.factor(with(subplot_data, paste(plot_div, plot_heat, sep="_")))

pdf(file="Outputs/Supps/Figure_S5_sown_CWMs_by_treatment.pdf", height=9,width=8.5, useDingbats=F)
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

with(subplot_data, plot(subplot.sown.log.sla.cwm ~ div_heat, cex.lab=1.25, cex.axis=1, xaxt="n", xlab="", ylab="Sown log(SLA) CWM", ylim=c(2.25, 3.1)))
axis(side=1, at=c(1:9), labels=rep("", times=9))
mtext(at=c(1:9), text=rep(c("Control", "Low", "High"), 3), side=1, line=0.8, cex=0.7)
# Second level lines
axis(side=1, at=c(1:3), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(4:6), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(7:9), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
# Second level labels
mtext(at=c(2,5,8), text=c("1 Species","4 Species", "16 Species"), side=1, line=3.2, cex=1)
mtext("(a)", line=0.5, side=3, adj=0, cex=1)
arrows(x0 = c(1,4,7), y0 = rep(3.05,3), x1 = c(3,6,9), y1 = rep(3.05,3), length = 0.05, angle = 90, code = 3)
text(c(2,5,8), rep(3.05,3), labels=c("N.S.", "N.S.", "N.S."), pos=3)

with(subplot_data, plot(subplot.sown.log.sm.cwm ~ div_heat, cex.lab=1.25, cex.axis=1, xaxt="n", xlab="", ylab="Sown log(Seed Mass) CWM", ylim=c(-9, -3.1)))
axis(side=1, at=c(1:9), labels=rep("", times=9))
mtext(at=c(1:9), text=rep(c("Control", "Low", "High"), 3), side=1, line=0.8, cex=0.7)
# Second level lines
axis(side=1, at=c(1:3), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(4:6), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(7:9), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
# Second level labels
mtext(at=c(2,5,8), text=c("1 Species","4 Species", "16 Species"), side=1, line=3.2, cex=1)
mtext("(b)", line=0.5, side=3, adj=0, cex=1)
arrows(x0 = c(1,4,7), y0 = rep(-3.45,3), x1 = c(3,6,9), y1 = rep(-3.45,3), length = 0.05, angle = 90, code = 3)
text(c(2,5,8), rep(-3.45, 3), labels=c("N.S.", "N.S.", "p=0.007"), pos=3)

with(subplot_data, plot(subplot.sown.sqrt.height.cwm ~ div_heat, cex.lab=1.25, cex.axis=1, xaxt="n", xlab="", ylab="Sown sqrt(Max. Height) CWM", ylim=c(3, 13.1)))
axis(side=1, at=c(1:9), labels=rep("", times=9))
mtext(at=c(1:9), text=rep(c("Control", "Low", "High"), 3), side=1, line=0.8, cex=0.7)
# Second level lines
axis(side=1, at=c(1:3), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(4:6), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(7:9), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
# Second level labels
mtext(at=c(2,5,8), text=c("1 Species","4 Species", "16 Species"), side=1, line=3.2, cex=1)
mtext("(c)", line=0.5, side=3, adj=0, cex=1)
arrows(x0 = c(1,4,7), y0 = rep(12.5,3), x1 = c(3,6,9), y1 = rep(12.5,3), length = 0.05, angle = 90, code = 3)
text(c(2,5,8), rep(12.5,3), labels=c("N.S.", "N.S.", "p=0.069"), pos=3)

with(subplot_data, plot(subplot.sown.log.ldmc.cwm ~ div_heat, cex.lab=1.25, cex.axis=1, xaxt="n", xlab="", ylab="Sown log(LDMC) CWM", ylim=c(4.9, 6.4)))
axis(side=1, at=c(1:9), labels=rep("", times=9))
mtext(at=c(1:9), text=rep(c("Control", "Low", "High"), 3), side=1, line=0.8, cex=0.7)
# Second level lines
axis(side=1, at=c(1:3), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(4:6), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
axis(side=1, at=c(7:9), xlim=c(1:3), labels=c("","",""), line=2.65, tck=0.025, lwd=1.25)
# Second level labels
mtext(at=c(2,5,8), text=c("1 Species","4 Species", "16 Species"), side=1, line=3.2, cex=1)
mtext("(d)", line=0.5, side=3, adj=0, cex=1)
arrows(x0 = c(1,4,7), y0 = rep(6.325,3), x1 = c(3,6,9), y1 = rep(6.325,3), length = 0.05, angle = 90, code = 3)
text(c(2,5,8), rep(6.325,3), labels=c("N.S.", "N.S.", "N.S."), pos=3)

dev.off()## must close pdf plot comman so that it saves to file



################
### FIGURE S6###
################
## color blind safe
c4_col <- "#56B4E9"
c3_col <- "#D55E00"  
forb_col <- "#F0E442"
legume_col <- "#009E73"

pdf(file=paste("Outputs/Supps/Figure_S6.pdf", sep=""),height=9,width=8.5,useDingbats=F)   
par(mfrow=c(2,2), mar=c(4.5, 4, 4, 2), mgp=c(2.6,0.8,0)) 

mc4.subplot.sown<-lmer(log(subplot.colonist.c4.grass.Mass.g.m2+1) ~ sqrt(subplot.sown.c4.grass.Mass.g.m2) + sqrt(subplot.sown.c3.grass.Mass.g.m2) + sqrt(subplot.sown.forb.Mass.g.m2) + sqrt(subplot.sown.legume.Mass.g.m2) + (1|Plot), subplot_data)
summary(mc4.subplot.sown)
with(subplot_data, plot(subplot.colonist.c4.grass.Mass.g.m2+1 ~ jitter(sqrt(subplot.sown.c4.grass.Mass.g.m2), amount=0.5), pch=19, col=c4_col, log="y",
                        ylab="C4 grass colonist biomass + 1", xlab="sqrt(sown biomass)", ylim=c(1,150)))
with(subplot_data, points(jitter(sqrt(subplot.sown.c3.grass.Mass.g.m2), amount=0.5), subplot.colonist.c4.grass.Mass.g.m2+1, pch=19, col=c3_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.forb.Mass.g.m2), amount=0.5), subplot.colonist.c4.grass.Mass.g.m2+1, pch=19, col=forb_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.legume.Mass.g.m2), amount=0.5), subplot.colonist.c4.grass.Mass.g.m2+1, pch=19, col=legume_col))
text(c(17,17,17,17), exp(c(5,4.75,4.5,4.25)), labels=c("sown C4 grass biomass", "sown C3 grass biomass", "sown forb biomass", "sown legume biomass*"), col=c(c4_col, c3_col, forb_col, legume_col))

curve(exp(cbind(1, x, mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc4.subplot.sown)), add=T, col=c4_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), x,
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc4.subplot.sown)), 
      from=0, to=13, add=T, col=c3_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                x, mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc4.subplot.sown)), 
      from=0, to=19, add=T, col=forb_col, lwd=3)								

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), x)%*%fixef(mc4.subplot.sown)), add=T, col=legume_col, lwd=3)
mtext("(a) C4 grass colonist biomass", line=0.5, side=3, adj=0, cex=1)


mc3.subplot.sown<-lmer(log(subplot.colonist.c3.grass.Mass.g.m2+1) ~ sqrt(subplot.sown.c4.grass.Mass.g.m2) + sqrt(subplot.sown.c3.grass.Mass.g.m2) + sqrt(subplot.sown.forb.Mass.g.m2) + sqrt(subplot.sown.legume.Mass.g.m2) + (1|Plot), subplot_data)
summary(mc3.subplot.sown)
with(subplot_data, plot(subplot.colonist.c3.grass.Mass.g.m2+1 ~ jitter(sqrt(subplot.sown.c4.grass.Mass.g.m2), amount=0.5), pch=19, col=c4_col, log="y",
                        ylab="C3 grass colonist biomass + 1", xlab="sqrt(sown biomass)", ylim=c(1,60)))
with(subplot_data, points(jitter(sqrt(subplot.sown.c3.grass.Mass.g.m2), amount=0.5), subplot.colonist.c3.grass.Mass.g.m2+1, pch=19, col=c3_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.forb.Mass.g.m2), amount=0.5), subplot.colonist.c3.grass.Mass.g.m2+1, pch=19, col=forb_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.legume.Mass.g.m2), amount=0.5), subplot.colonist.c3.grass.Mass.g.m2+1, pch=19, col=legume_col))
text(c(17,17,17,17), exp(c(4,3.8,3.6,3.4)), labels=c("sown C4 grass biomass*", "sown C3 grass biomass", "sown forb biomass", "sown legume biomass"), col=c(c4_col, c3_col, forb_col, legume_col))
curve(exp(cbind(1, x, mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc3.subplot.sown)), add=T, col=c4_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), x,
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc3.subplot.sown)), 
      from=0, to=13, add=T, col=c3_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                x, mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mc3.subplot.sown)), 
      from=0, to=19, add=T, col=forb_col, lwd=3)								

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), x)%*%fixef(mc3.subplot.sown)), add=T, col=legume_col, lwd=3)
mtext("(b) C3 grass colonist biomass", line=0.5, side=3, adj=0, cex=1)


mforb.subplot.sown<-lmer(log(subplot.colonist.forb.Mass.g.m2+1) ~ sqrt(subplot.sown.c4.grass.Mass.g.m2) + sqrt(subplot.sown.c3.grass.Mass.g.m2) + sqrt(subplot.sown.forb.Mass.g.m2) + sqrt(subplot.sown.legume.Mass.g.m2) + (1|Plot), subplot_data)
summary(mforb.subplot.sown)
with(subplot_data, plot(subplot.colonist.forb.Mass.g.m2+1 ~ jitter(sqrt(subplot.sown.c4.grass.Mass.g.m2), amount=0.5), pch=19, col=c4_col, log="y",
                        ylab="Forb colonist biomass + 1", xlab="sqrt(sown biomass)", ylim=c(1,200)))
with(subplot_data, points(jitter(sqrt(subplot.sown.c3.grass.Mass.g.m2), amount=0.5), subplot.colonist.forb.Mass.g.m2+1, pch=19, col=c3_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.forb.Mass.g.m2), amount=0.5), subplot.colonist.forb.Mass.g.m2+1, pch=19, col=forb_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.legume.Mass.g.m2), amount=0.5), subplot.colonist.forb.Mass.g.m2+1, pch=19, col=legume_col))
text(c(17,17,17,17), exp(c(5.25,5,4.75, 4.5)), labels=c("sown C4 grass biomass***", "sown C3 grass biomass***", "sown forb biomass", "sown legume biomass"), col=c(c4_col, c3_col, forb_col, legume_col))
curve(exp(cbind(1, x, mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mforb.subplot.sown)), add=T, col=c4_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), x,
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mforb.subplot.sown)), 
      from=0, to=13, add=T, col=c3_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                x, mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mforb.subplot.sown)), 
      from=0, to=19, add=T, col=forb_col, lwd=3)								

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), x)%*%fixef(mforb.subplot.sown)), add=T, col=legume_col, lwd=3)
mtext("(c) Forb colonist biomass", line=0.5, side=3, adj=0, cex=1)


mlegume.subplot.sown<-lmer(log(subplot.colonist.legume.Mass.g.m2+1) ~ sqrt(subplot.sown.c4.grass.Mass.g.m2) + sqrt(subplot.sown.c3.grass.Mass.g.m2) + sqrt(subplot.sown.forb.Mass.g.m2) + sqrt(subplot.sown.legume.Mass.g.m2) + (1|Plot), subplot_data)
summary(mlegume.subplot.sown)
with(subplot_data, plot(subplot.colonist.legume.Mass.g.m2+1 ~ jitter(sqrt(subplot.sown.c4.grass.Mass.g.m2), amount=0.5), pch=19, col=c4_col, log="y",
                        ylab="Legume colonist biomass + 1", xlab="sqrt(sown biomass)", ylim=c(1,30)))
with(subplot_data, points(jitter(sqrt(subplot.sown.c3.grass.Mass.g.m2), amount=0.5), subplot.colonist.legume.Mass.g.m2+1, pch=19, col=c3_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.forb.Mass.g.m2), amount=0.5), subplot.colonist.legume.Mass.g.m2+1, pch=19, col=forb_col))
with(subplot_data, points(jitter(sqrt(subplot.sown.legume.Mass.g.m2), amount=0.5), subplot.colonist.legume.Mass.g.m2+1, pch=19, col=legume_col))
text(c(17,17,17,17), exp(c(3.4,3.2,3,2.8)), labels=c("sown C4 grass biomass*", "sown C3 grass biomass", "sown forb biomass", "sown legume biomass"), col=c(c4_col, c3_col, forb_col, legume_col))
curve(exp(cbind(1, x, mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mlegume.subplot.sown)), add=T, col=c4_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), x,
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mlegume.subplot.sown)), 
      from=0, to=13, add=T, col=c3_col, lwd=3)

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                x, mean(sqrt(subplot_data$subplot.sown.legume.Mass.g.m2), na.rm=T))%*%fixef(mlegume.subplot.sown)), 
      from=0, to=19, add=T, col=forb_col, lwd=3)								

curve(exp(cbind(1,mean(sqrt(subplot_data$subplot.sown.c4.grass.Mass.g.m2), na.rm=T), mean(sqrt(subplot_data$subplot.sown.c3.grass.Mass.g.m2), na.rm=T),
                mean(sqrt(subplot_data$subplot.sown.forb.Mass.g.m2)), x)%*%fixef(mlegume.subplot.sown)), add=T, col=legume_col, lwd=3)
mtext("(d) Legume colonist biomass", line=0.5, side=3, adj=0, cex=1)

dev.off()


