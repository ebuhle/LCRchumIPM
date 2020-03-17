#setwd("C:\\data\\BPA Projects\\Chum\\chum fecundity")
setwd("~/BentleyKT/Projects/Chum/Analysis/Fecundity")
library(R2jags)
library(rjags)

#plotfunction========================================================
JAGSplot<-function(JAGS.sims.list){
  dat<-JAGS.sims.list
  pdf(paste("JAGSplots",Sys.Date(),".pdf",sep="_"))
  for(i in 1:length(names(dat))){
    tdat<-dat[[i]]
    if(length(dim(tdat))<4){
      if(length(dim(tdat))==3){
        for(j in 1:dim(tdat)[3]){
          plotdat<-as.matrix(tdat[,,j])
          index<-j
          bp<-boxplot.matrix(plotdat,outline=F,plot=F)
          bp$stats[1,]<-apply(plotdat,2,function(x) quantile(x,0.025))
          bp$stats[2,]<-apply(plotdat,2,function(x) quantile(x,0.25))
          bp$stats[3,]<-apply(plotdat,2,function(x) quantile(x,0.5))
          bp$stats[4,]<-apply(plotdat,2,function(x) quantile(x,0.75))
          bp$stats[5,]<-apply(plotdat,2,function(x) quantile(x,0.975))
          bp$out<-NULL
          bxp(bp,outline=F)
          mtext(3,text=paste(names(dat)[[i]],index))
        }
      }else{plotdat<-tdat;index<-""
      bp<-boxplot.matrix(plotdat,outline=F,plot=F)
      bp$stats[1,]<-apply(plotdat,2,function(x) quantile(x,0.025))
      bp$stats[2,]<-apply(plotdat,2,function(x) quantile(x,0.25))
      bp$stats[3,]<-apply(plotdat,2,function(x) quantile(x,0.5))
      bp$stats[4,]<-apply(plotdat,2,function(x) quantile(x,0.75))
      bp$stats[5,]<-apply(plotdat,2,function(x) quantile(x,0.975))
      bp$out<-NULL
      bxp(bp,outline=F)
      mtext(3,text=paste(names(dat)[[i]],index))
      }
    }
  }
  dev.off()
}

#------------------------------------------------------------------------------------------------------------
# Bring in data and format
#------------------------------------------------------------------------------------------------------------
all.fec.dat<-read.csv("Data_ChumFecundity_fromHatcheryPrograms_2017-01-25.csv", header=TRUE, na.strings=c("NA", ""))
head(all.fec.dat)

colnames(all.fec.dat)[which(names(all.fec.dat) == "Reproductive.Effort")] <- "Perct.Reprod.Eff"
colnames(all.fec.dat)[which(names(all.fec.dat) == "Estimated.Fecundity")] <- "Fecundity"
colnames(all.fec.dat)[which(names(all.fec.dat) == "Green.egg.avg.weight")] <- "Avg.Egg.Wt"
is.numeric(all.fec.dat$Estimated.Fecundity)

#Remove fish with no fecundity estimate
dat<-all.fec.dat[!is.na(all.fec.dat$Fecundity),]

#Remove females with RE <16% from dataset (and also obvious outliers = females with RE >30%)
dat<-dat[dat$Perct.Reprod.Eff>=16 & dat$Perct.Reprod.Eff<30, ]

#remove fish with no age
dat<-dat[dat$Age%in%c(3,4,5),]

names(dat)
#-----------------------------------------------------------------------------------------------------------
#jags data
#-----------------------------------------------------------------------------------------------------------
JagsDat<-list(
  obs=dim(dat)[1],
  fecundity=dat$Fecundity,
  age=as.numeric(as.factor(dat$Age)),
  yr=as.numeric(as.factor(dat$BY)),
  basin=as.numeric(as.factor(dat$Stock)),
  ages=length(unique(dat$Age)),
  basins=length(unique(dat$Stock)),
  years=length(unique(dat$BY))
)
#parameters
pars<-c(
  #global mean
  "global_mu",

  #age fixed effects
  "age_eff",
  
  #year random effect
  "yr_eff",
  "sigma_yr",
  
  #basin effect
  "basin_eff",
  "sigma_basin",
  
  #year and basin-specific residuals
  "resid",
  "sigma_resid",
  
  #predictions
  "unknown_resid", 
  "unknown_basin_eff",
  "unknown_year_eff",
  "unk_basin_unknown_year_age",
  "unk_basin_specific_year_age",
  "known_basin_specific_year_unknown_resid_age",
  "known_basin_unknown_year_unknown_resid_age",
  
  #pointwise log likelihood
  "fecundity_lik"
  
)

#sometimes needed to stop random error message in JAGS
runif(1)

#run model
start.time<-Sys.time()
print(start.time)
jagsfit <- jags.parallel(JagsDat, 
                         inits=NULL, 
                         model.file="Chum_Fecundity_hierarchical_1.31.2017.txt",
                         n.chains=3, 
                         n.thin=10, 
                         n.burnin=5000, 
                         n.iter=30000,
                         parameters.to.save=pars)
#save model run R workspace
save(jagsfit,file=paste("Chum_Fecundity_hierarchical_",Sys.Date(),".rda",sep=""))
jagsfit.mcmc <- as.mcmc(jagsfit)
end.time<-Sys.time();print(end.time-start.time)

#convergence diagnostics
#codamenu()
#jagsfit.mcmc

#look at estimates
write.csv(jagsfit$BUGSoutput$summary,paste("Chum_Fecundity_hierarchical_summary_",Sys.Date(),".csv",sep=""))
JAGSplot(jagsfit$BUGSoutput$sims.list)

#==========================================
#waic and loo
#============
attach.jags(jagsfit)
PLL<-log(fecundity_lik)
if(!require(loo))install.packages("loo") 
library(loo) 
waic(PLL)
loo(PLL)

tau<-.1
sd<-1/sqrt(tau)
plot(density(exp(rnorm(100000,0,sd))),xlim=c(0,10000))

