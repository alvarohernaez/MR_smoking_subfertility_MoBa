rm(list=ls())

# Load required libraries #

# Shortcut format functions #

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.0001,signif(x,1),
                   ifelse(abs(x)<0.001,signif(x,1),
                          ifelse(abs(x)<0.1,round(x,3),
                                 ifelse(abs(x)<1,round(x,2),signif(x,3)))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",round(x,3)))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

pval_ic_guapa <- function(ratio, lo, up) 
{
  z<-qnorm(1-0.05/2)
  x<-log(ratio)/((log(up)-log(lo))/(2*z))
  pval<-pval_guapa(exp(-0.717*x-(0.416*x^2)))
  return(pval)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

setwd("N:/data/durable/Projects/Magnus_MR_smoking")


#####################
### MAIN ANALYSES ###
#####################

dir.create("./Outputs")
dir.create("./Outputs/descriptive")
dir.create("./Outputs/results")
dir.create("./Outputs/results/nlmr_source")

setwd("N:/data/durable/Projects/Magnus_MR_smoking/Outputs")


##############################
### POPULATION DESCRIPTION ###
##############################

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")

datx<-subset2(dat,"dat$smkinit_na_mom==0 & dat$preg_plan==1")
datx$parity<-with(datx,ifelse(parity==0,1,2))

xxx<-datx[,c("agedelivery_mom","eduyears_mom","bmi_mom","parity","subfertility_12plus",
             "smkinit_mom","agesmk_mom","cigday_mom","smkces_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subfertility_12plus,
                               xxx, method=c("smkinit_mom"=3,"agesmk_mom"=2,"cigday_mom"=2,"smkces_mom"=3,
                                             "bmi_mom"=2,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

comp<-NULL
comp<-createTable(compareGroups(subfertility_12plus~.
                                -sel,
                                xxx, method=c("smkinit_mom"=3,"agesmk_mom"=2,"cigday_mom"=2,"smkces_mom"=3,
                                              "bmi_mom"=2,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab1)<-c("Mothers-All","Mothers-Non-subfertile","Mothers-Subfertile","Mothers-P-value","Mothers-N")
write.table(tab1,file="./descriptive/descriptive_mothers.csv",sep=";",col.names=NA)


datx<-subset2(dat,"dat$smkinit_na_dad==0 & dat$preg_plan==1")
datx$cigday_dad<-with(datx,ifelse(cigday_na_dad==0,cigday_dad,NA))
datx$parity<-with(datx,ifelse(parity==0,1,2))

xxx<-datx[,c("agedelivery_dad","eduyears_dad","bmi_dad","parity","smkinit_dad","cigday_dad","subfertility_12plus")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subfertility_12plus,
                               xxx, method=c("smkinit_dad"=3,"cigday_dad"=2,"bmi_dad"=2,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

comp<-NULL
comp<-createTable(compareGroups(subfertility_12plus~.
                                -sel,
                                xxx, method=c("smkinit_dad"=3,"cigday_dad"=2,"bmi_dad"=2,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab2<-NULL
tab2<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab2)<-c("Fathers-All","Fathers-Non-subfertile","Fathers-Subfertile","Fathers-P-value","Fathers-N")
write.table(tab2,file="./descriptive/descriptive_fathers.csv",sep=";",col.names=NA)


### SELECTION BIAS - INCLUDED vs NON-INCLUDED ###

datx<-dat
datx$parity<-with(datx,ifelse(parity==0,1,2))
datx$smkinit_na_mom<-with(datx,ifelse(preg_plan==0,1,smkinit_na_mom))

xxx<-datx[,c("smkinit_na_mom","agedelivery_mom","eduyears_mom","bmi_mom",
             "parity","subfertility_12plus",
             "smkinit_mom","agesmk_mom","cigday_mom","smkces_mom")]

all<-NULL
all<-createTable(compareGroups(smkinit_na_mom~.,
                               xxx, method=c("smkinit_mom"=3,"agesmk_mom"=2,"cigday_mom"=2,"smkces_mom"=3,
                                             "bmi_mom"=2,"parity"=3,"subfertility_12plus"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_mothers.csv",sep=";",col.names=NA)

datx$smkinit_na_dad<-with(datx,ifelse(preg_plan==0,1,smkinit_na_dad))

xxx<-datx[,c("smkinit_na_dad","agedelivery_dad","eduyears_dad","bmi_dad",
             "parity","subfertility_12plus",
             "smkinit_dad","cigday_dad")]

all<-NULL
all<-createTable(compareGroups(smkinit_na_dad~.,
                               xxx, method=c("smkinit_dad"=3,"cigday_dad"=2,
                                             "bmi_dad"=2,"parity"=3,"subfertility_12plus"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_fathers.csv",sep=";",col.names=NA)


###################################################
### ASSOCIATION BETWEEN SMOKING TRAITS AND GRSs ###
###################################################

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")
library(pROC)

# ASSOCIATION BETWEEN GRSs AND BINARY VARIABLES #

vars01<-c("smkinit_mom","smkinit_dad","smkces_mom")
vars02<-c("smkinit_grs_mom","smkinit_grs_dad","smkces_grs_mom")
vars04<-c("m_id_2374","f_id_2374","m_id_2374")
vars05<-c("smkinit_na_mom","smkinit_na_dad","smkces_na_mom")
vars06<-c("355","355","23")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars05[i]]==0 & dat$preg_plan==1")
  mod01<-glm(formula=as.factor(datx[,vars01[i]])~datx[,vars02[i]],
             data=datx, family="binomial")
  estimate<-as.numeric(summary(mod01)$coefficients[2,1])
  se<-as.numeric(summary(mod01)$coefficients[2,2])
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(estimate))
  ic95a<-round(exp(estimate-(z*se)),3)
  ic95b<-round(exp(estimate+(z*se)),3)
  coef<-ic_guapa(hr,ic95a,ic95b)
  llrtest<-round(as.numeric(anova(mod01,test=c("F"))[2,5]),0) # log-likelihood ratio test (~F in binomial GLM)
  pseudor2<-paste("'",guapa(NagelkerkeR2(mod01)$R2*100),sep="") # pseudoR2 by the Nagelkerke method
  n_snps<-vars06[i]
  datx$prob<-predict.glm(mod01,type=c("response"))
  roc_auc<-guapa(roc(as.factor(datx[,vars01[i]])~prob,data=datx)$auc)
  mean_grs<-paste(guapa(mean(datx[,vars02[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps,mean_grs,coef,llrtest,roc_auc,pseudor2))
}

rownames(tab)<-vars01
write.table(tab,file="./descriptive/assoc_grs_binary.csv",sep=";",col.names=NA)


# ASSOCIATION BETWEEN GRS AND CONTINUOUS VARIABLES #

vars01<-c("cigday_mom","cigday_dad","agesmk_mom")
vars02<-c("cigday_grs_mom","cigday_grs_dad","agesmk_grs_mom")
vars03<-c("m_id_2374","f_id_2374","m_id_2374")
vars04<-c("50","50","10")
vars05<-c("cigday_na_mom","cigday_na_dad","agesmk_na_mom")
vars07<-c("Cigarettes per week-GRS, mothers","Cigarettes per week-GRS, fathers",
          "Age of initiation of smoking-GRS, mothers")
vars08<-c("Cigarettes per week","Cigarettes per week",
          "Age of initiation of smoking (years)")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars05[i]]==0 & dat$preg_plan==1")
  mod01<-lm_robust(datx[,vars01[i]]~datx[,vars02[i]], 
                   data=datx, clusters=datx[,vars03[i]], se_type="stata")
  coef<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
              guapa(summary(mod01)$coefficients[2,5]),"; ",
              guapa(summary(mod01)$coefficients[2,6]),"]",sep="")
  fstat<-round(mod01$fstatistic[1],0)
  r2<-paste("'",guapa(mod01$adj.r.squared*100),sep="")
  n_snps<-vars04[i]
  mean_grs<-paste(guapa(mean(datx[,vars02[i]]))," (",guapa(sd(datx[,vars02[i]])),")",sep="")
  tab<-rbind(tab,cbind(n_snps,mean_grs,coef,fstat,r2))
  
  aaa<-datx[,vars02[i]]
  model<-glm(datx[,vars01[i]]~bs(aaa),data=datx)
  mod01<-glm(datx[,vars01[i]]~aaa, family="gaussian", data=datx)
  p_lin<-guapa(summary(mod01)$coefficients[2,4])
  p_nonlin<-guapa(lrtest(model,mod01)[2,5])
  
  z<-qnorm(1-0.05/2)
  res<-termplot(model,term='bs(aaa)',rug=FALSE,se=TRUE,plot=FALSE)
  res<-res$aaa
  corr<-res$y[1]
  ajus<-mean(datx[,vars01[i]][which(datx[,vars02[i]]==res$x[1])],na.rm=TRUE)
  res$y<-res$y-corr+ajus
  ci<-cbind(res$y,res$y-z*res$se,res$y+z*res$se)
  alc_n<-paste("./descriptive/assoc_grs_",vars01[i],".jpg",sep="")
  xlabel<-vars07[i]
  ylabel<-vars08[i]
  
  jpeg(filename=alc_n,width=3000,height=3000,res=600,pointsize=11.5)
  par(las=1,cex=1,mar=c(5,5,2,2))
  matplot(res$x, ci, lty=c(1,0,0),lwd=2,type="l",col="black",xlab=xlabel,ylab=ylabel)
  polygon(c(res$x,rev(res$x)),c(ci[,2],rev(ci[,3])),col=grey(0.8),border="white")
  matpoints(res$x,ci,lty=c(1,0,0),lwd=2,type="l",col="black")
  abline(h=0, lty=3, col="black")
  abline(v=0, lty=3, col="black")
  leg<-paste("Linear component: P-value = ",p_lin,"\nNon-linear component: P-value = ",p_nonlin,sep="")
  legend("bottomleft",legend=leg, bty="n", cex=0.7, inset=c(0.025,0.025))
  dev.off()
}

rownames(tab)<-vars01
write.table(tab,file="./descriptive/assoc_grs_continuous.csv",sep=";",col.names=NA)


##########################################
### ASSOCIATION OF GRS WITH COVARIATES ###
##########################################

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")
dat$parity_cont<-with(dat,ifelse(parity_cont==4,3,parity_cont))

# Horizontal pleiotropy - Mothers #

vars01<-c("smkinit_grs_mom_z","smkces_grs_mom_z","cigday_grs_mom_z","agesmk_grs_mom_z")
vars02<-c("smkinit_na_mom","smkces_na_mom","cigday_na_mom","agesmk_na_mom")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars02[i]]==0 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  mod01<-lm_robust(agedelivery_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod02<-lm_robust(eduyears_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod03<-lm_robust(bmi_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod04<-lm_robust(parity_cont~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  coef01<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
                guapa(summary(mod01)$coefficients[2,5]),"; ",
                guapa(summary(mod01)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod01)$coefficients[2,4]),")",sep="")
  coef02<-paste(guapa(summary(mod02)$coefficients[2,1])," [",
                guapa(summary(mod02)$coefficients[2,5]),"; ",
                guapa(summary(mod02)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod02)$coefficients[2,4]),")",sep="")
  coef03<-paste(guapa(summary(mod03)$coefficients[2,1])," [",
                guapa(summary(mod03)$coefficients[2,5]),"; ",
                guapa(summary(mod03)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod03)$coefficients[2,4]),")",sep="")
  coef04<-paste(guapa(summary(mod04)$coefficients[2,1])," [",
                guapa(summary(mod04)$coefficients[2,5]),"; ",
                guapa(summary(mod04)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod04)$coefficients[2,4]),")",sep="")
  
  tab<-rbind(tab,cbind(coef01,coef02,coef03,coef04))
}

rownames(tab)<-vars01
colnames(tab)<-c("age","eduyears","bmi","parity")
write.table(tab,file="./descriptive/hp_mom.csv",sep=";",col.names=NA)


vars01<-c("smkces_grs_mom_z","cigday_grs_mom_z","agesmk_grs_mom_z")
vars02<-c("smkces_na_mom","cigday_na_mom","agesmk_na_mom")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat$smkinit_na_mom==0 & dat[,vars02[i]]==1 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  mod01<-lm_robust(agedelivery_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod02<-lm_robust(eduyears_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod03<-lm_robust(bmi_mom~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod04<-lm_robust(parity_cont~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  coef01<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
                guapa(summary(mod01)$coefficients[2,5]),"; ",
                guapa(summary(mod01)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod01)$coefficients[2,4]),")",sep="")
  coef02<-paste(guapa(summary(mod02)$coefficients[2,1])," [",
                guapa(summary(mod02)$coefficients[2,5]),"; ",
                guapa(summary(mod02)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod02)$coefficients[2,4]),")",sep="")
  coef03<-paste(guapa(summary(mod03)$coefficients[2,1])," [",
                guapa(summary(mod03)$coefficients[2,5]),"; ",
                guapa(summary(mod03)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod03)$coefficients[2,4]),")",sep="")
  coef04<-paste(guapa(summary(mod04)$coefficients[2,1])," [",
                guapa(summary(mod04)$coefficients[2,5]),"; ",
                guapa(summary(mod04)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod04)$coefficients[2,4]),")",sep="")
  
  tab<-rbind(tab,cbind(coef01,coef02,coef03,coef04))
}

rownames(tab)<-vars01
colnames(tab)<-c("age","eduyears","bmi","parity")
write.table(tab,file="./descriptive/hp_mom_norelevance.csv",sep=";",col.names=NA)


# Horizontal pleiotropy - Fathers #

vars01<-c("smkinit_grs_dad_z","cigday_grs_dad_z")
vars02<-c("smkinit_na_dad","cigday_na_dad")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars02[i]]==0 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  mod01<-lm_robust(agedelivery_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod02<-lm_robust(eduyears_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod03<-lm_robust(bmi_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod04<-lm_robust(parity_cont~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  coef01<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
                guapa(summary(mod01)$coefficients[2,5]),"; ",
                guapa(summary(mod01)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod01)$coefficients[2,4]),")",sep="")
  coef02<-paste(guapa(summary(mod02)$coefficients[2,1])," [",
                guapa(summary(mod02)$coefficients[2,5]),"; ",
                guapa(summary(mod02)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod02)$coefficients[2,4]),")",sep="")
  coef03<-paste(guapa(summary(mod03)$coefficients[2,1])," [",
                guapa(summary(mod03)$coefficients[2,5]),"; ",
                guapa(summary(mod03)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod03)$coefficients[2,4]),")",sep="")
  coef04<-paste(guapa(summary(mod04)$coefficients[2,1])," [",
                guapa(summary(mod04)$coefficients[2,5]),"; ",
                guapa(summary(mod04)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod04)$coefficients[2,4]),")",sep="")
  
  tab<-rbind(tab,cbind(coef01,coef02,coef03,coef04))
}

rownames(tab)<-vars01
colnames(tab)<-c("age","eduyears","bmi","parity")
write.table(tab,file="./descriptive/hp_dad.csv",sep=";",col.names=NA)


vars01<-c("cigday_grs_dad_z")
vars02<-c("cigday_na_dad")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat$smkinit_na_dad==0 & dat[,vars02[i]]==1 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  mod01<-lm_robust(agedelivery_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod02<-lm_robust(eduyears_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod03<-lm_robust(bmi_dad~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  mod04<-lm_robust(parity_cont~datx[,vars01[i]], 
                   data=datx, clusters=m_id_2374, se_type="stata")
  coef01<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
                guapa(summary(mod01)$coefficients[2,5]),"; ",
                guapa(summary(mod01)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod01)$coefficients[2,4]),")",sep="")
  coef02<-paste(guapa(summary(mod02)$coefficients[2,1])," [",
                guapa(summary(mod02)$coefficients[2,5]),"; ",
                guapa(summary(mod02)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod02)$coefficients[2,4]),")",sep="")
  coef03<-paste(guapa(summary(mod03)$coefficients[2,1])," [",
                guapa(summary(mod03)$coefficients[2,5]),"; ",
                guapa(summary(mod03)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod03)$coefficients[2,4]),")",sep="")
  coef04<-paste(guapa(summary(mod04)$coefficients[2,1])," [",
                guapa(summary(mod04)$coefficients[2,5]),"; ",
                guapa(summary(mod04)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod04)$coefficients[2,4]),")",sep="")
  
  tab<-rbind(tab,cbind(coef01,coef02,coef03,coef04))
}

rownames(tab)<-vars01
colnames(tab)<-c("age","eduyears","bmi","parity")
write.table(tab,file="./descriptive/hp_dad_norelevance.csv",sep=";",col.names=NA)


##########################################################################
### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##########################################################################

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")

dat$agesmk_mom_z<-as.numeric(scale(dat$agesmk_mom))
dat$parity_cont<-with(dat,ifelse(parity_cont==4,3,parity_cont))

vars00<-c("gen_pred_smkinit_mom","gen_pred_smkinit_dad",
          "gen_pred_agesmk_mom","gen_pred_smkces_mom",
          "gen_pred_cigday_mom","gen_pred_cigday_dad")
vars01<-c("gen_pred_smkinit_mom_z","gen_pred_smkinit_dad_z",
          "gen_pred_agesmk_mom_z","gen_pred_smkces_mom_z",
          "gen_pred_cigday_mom_z","gen_pred_cigday_dad_z")

for(i in 1:length(vars01))
{
  dat[,vars01[i]]<-as.numeric(with(dat,scale(dat[,vars00[i]])))
} 


vars02<-c("smkinit_mom","smkinit_dad","agesmk_mom_z","smkces_mom","cigday_mom_z","cigday_dad_z")
vars05<-c("agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_mom","agedelivery_mom","agedelivery_dad")
vars06<-c("eduyears_mom","eduyears_dad","eduyears_mom","eduyears_mom","eduyears_mom","eduyears_dad")
vars07<-c("bmi_mom","bmi_dad","bmi_mom","bmi_mom","bmi_mom","bmi_dad")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_mom","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_mom","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_mom","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_mom","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_mom","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_mom","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_mom","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_mom","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_mom","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_mom","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","m_id_2374","m_id_2374","f_id_2374")
vars19<-c("smkinit_na_mom","smkinit_na_dad","agesmk_na_mom","smkces_na_mom","cigday_na_mom","cigday_na_dad")
vars20<-c("SmkInit, mothers","SmkInit, fathers",
          "AgeSmk, mothers","SmkCes, mothers",
          "CigDay, mothers","CigDay, fathers")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  mod01<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars02[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod01)[2,1])
  se<-as.numeric(summary(mod01)[2,2])
  coef01<-risk_se_ic_guapa(estimate,se)
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars02[i]]
                               +datx[,vars05[i]]+datx[,vars06[i]]+datx[,vars07[i]]+parity_cont,
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod02)[2,1])
  se<-as.numeric(summary(mod02)[2,2])
  coef02<-risk_se_ic_guapa(estimate,se)
  forest_mv<-risk_se_ic_guapa2(estimate,se)
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  
  aaa<-datx[,vars02[i]]
  mod_base<-glm(formula=as.factor(subfertility_12plus)~datx[,vars05[i]]
                +datx[,vars06[i]]+datx[,vars07[i]]+parity_cont,
                data=datx, family="binomial")
  mod_lin<-glm(formula=as.factor(subfertility_12plus)~aaa
               +datx[,vars05[i]]+datx[,vars06[i]]+datx[,vars07[i]]+parity_cont,
               data=datx, family="binomial")
  mod_nlin<-glm(formula=as.factor(subfertility_12plus)~bs(aaa,degree=3)
                +datx[,vars05[i]]+datx[,vars06[i]]+datx[,vars07[i]]+parity_cont,
                data=datx, family="binomial")
  p_lin<-pval_guapa(lrtest(mod_base,mod_lin)[2,5])
  p_nonlin<-pval_guapa(lrtest(mod_base,mod_nlin)[2,5])
  p_lrtest<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  forest_mr<-risk_se_ic_guapa2(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,
                       p_lin,p_nonlin,p_lrtest,
                       coef03,pval03,coef04,pval04,forest_mv,forest_mr))
}

colnames(tab)<-c("Det-OR (raw)","Det-pval (raw)","Det-OR (adj.)","Det-pval (adj.)",
                 "P-linear (adj.)","P-nonlinear (adj.)","P-lrtest (adj.)",
                 "GRS-OR (raw)","GRS-pval (raw)","GRS-OR (adj.)","GRS-pval (adj.)",
                 "forest_mv","forest_mr")
rownames(tab)<-vars20
write.table(tab,file="./results/linear.csv",sep=";",col.names=NA)


### NO RELEVANCE ANALYSES ###

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")

vars00<-c("agesmk_grs_mom","smkces_grs_mom","cigday_grs_mom","cigday_grs_dad")
vars01<-c("agesmk_grs_mom_z","smkces_grs_mom_z","cigday_grs_mom_z","cigday_grs_dad_z")

for(i in 1:length(vars01))
{
  dat[,vars01[i]]<-as.numeric(with(dat,scale(dat[,vars00[i]])))
} 


vars08<-c("pc1_mom","pc1_mom","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_mom","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_mom","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_mom","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_mom","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_mom","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_mom","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_mom","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_mom","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_mom","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","m_id_2374","m_id_2374","f_id_2374")
vars19<-c("agesmk_na_mom","smkces_na_mom","cigday_na_mom","cigday_na_dad")
vars20<-c("smkinit_na_mom","smkinit_na_mom","smkinit_na_mom","smkinit_na_dad")
vars21<-c("AgeSmk, mothers","SmkCes, mothers","CigDay, mothers","CigDay, fathers")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==1 & dat[,vars20[i]]==0 & !is.na(dat[,vars01[i]]) & dat$preg_plan==1")
  
  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-rbind(tab,cbind(coef03,pval03,coef04,pval04))
}

colnames(tab)<-c("GRS-OR (raw)","GRS-pval (raw)","GRS-OR (adj.)","GRS-pval (adj.)")
rownames(tab)<-vars21
write.table(tab,file="./results/no_relevance_linear.csv",sep=";",col.names=NA)


#######################################################################
### LOGISTIC REGRESSION: NON-LINEAR ASSOCIATIONS (SMOOTHED SPLINES) ###
#######################################################################

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")
dat$agesmk_mom<-with(dat,ifelse(agesmk_mom<9,NA,
                                ifelse(agesmk_mom>38,NA,agesmk_mom)))

vars01<-c("agesmk_mom","cigday_mom","cigday_dad")
vars04<-c("agedelivery_mom","agedelivery_mom","agedelivery_dad")
vars05<-c("eduyears_mom","eduyears_mom","eduyears_dad")
vars06<-c("bmi_mom","bmi_mom","bmi_dad")
vars17<-c("Age of smoking initiation (years), mothers","Cigarettes per week, mothers","Cigarettes per week, fathers")
vars19<-c("agesmk_na_mom","cigday_na_mom","cigday_na_dad")
vars20<-c(16,0,0)
vars21<-c(9,0,0)
vars22<-c(38,150,180)
vars23<-c(0,0.9,0.9)
vars24<-c(2,1.5,1.5)

# EXECUTE MANUALLY REPEATING THE LOOPS i=1/2/3 #

tab<-NULL
for(i in 1:length(vars01))
  
{
  varstot<-c(vars01[i],vars04[i],vars05[i],vars06[i],vars19[i],
             "parity","subfertility_12plus","preg_plan")
  dat2<-na.omit(dat[,varstot])
  dat2<-subset2(dat2,"dat2[,vars19[i]]==0 & dat2$preg_plan==1")
  aaa<-dat2[,vars01[i]]
  
  mod01<-glm(formula=as.factor(subfertility_12plus)~bs(aaa,degree=3)
             +dat2[,vars04[i]]+dat2[,vars05[i]]+dat2[,vars06[i]]+parity,
             data=dat2, family="binomial")
  mod02<-glm(formula=as.factor(subfertility_12plus)~aaa
             +dat2[,vars04[i]]+dat2[,vars05[i]]+dat2[,vars06[i]]+parity,
             data=dat2, family="binomial")
  p_lin<-pval_guapa(summary(mod02)$coefficients[2,4])
  p_nonlin<-pval_guapa(lrtest(mod01,mod02)[2,5])
  
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,vars20[i])
  #min_val<-temp$x[which(temp$y==min(temp$y,na.rm=TRUE))]
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./results/spline_",vars01[i],".jpg",sep="")
  labely<-c("Subfertility (odds ratio, adjusted)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  riskmin<-plot.data[plot.data$y==min(plot.data$y,na.rm=TRUE),]
  riskmin<-paste(guapa(riskmin[1]),": ",ic_guapa(guapa(riskmin[2]),guapa(riskmin[3]),guapa(riskmin[4]))," (P=",
                 pval_ic_guapa(riskmin[2],riskmin[3],riskmin[4]),")",sep="")
  
  tab<-rbind(tab,cbind(p_lin,p_nonlin,riskmin))
  
  figure<-ggplot(plot.data, aes(x=x, ymax=5))
  figure<-figure + geom_hline(aes(yintercept=1), colour="black", linetype=2) + 
    geom_line(aes(y=yest), color="black") + 
    geom_line(aes(y=lci), color="grey") + 
    geom_line(aes(y=uci), color="grey") + 
    scale_x_continuous(limits = c(vars21[i],vars22[i])) +
    scale_y_continuous(limits = c(vars23[i],vars24[i])) +
    theme_bw() +
    labs(x=vars17[i],y=labely) +
    theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20),
          axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + 
    geom_point(aes(x=vars20[i], y=1),data=plot.data, colour="black", size=3) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  jpeg(filename=name,width = 8000, height = 8000, res=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
}

rownames(tab)<-vars01
write.table(tab,file="./results/spline_logreg_coefs.csv",sep=";")


##################################################
### MULTIVARIABLE MR / AGE-STRATIFIED ANALYSES ###
##################################################

########################
### MULTIVARIABLE MR ###
########################

### GENERATION OF ADJUSTED GENETICALLY-PREDICTED VALUES ###

setwd("N:/data/durable/Projects/Magnus_MR_smoking")

dir.create("./Outputs/sensitivity")
dir.create("./Outputs/sensitivity/mvmr")

setwd("N:/data/durable/Projects/Magnus_MR_smoking/Outputs/sensitivity/mvmr")

load("N:/data/durable/Projects/Magnus_MR_BMI/R/bmi_grs.RData")
bmi_grs<-dat
load("N:/data/durable/Projects/Magnus_MR_BMI/R/eduyears_grs.RData")
eduyears_grs<-dat
load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")

bmi_grs<-rename.vars(bmi_grs,from=c("id","bmi_grs"),to=c("sentrixid_mom","bmi_grs_mom"))
dat<-merge2(dat,bmi_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
bmi_grs<-rename.vars(bmi_grs,from=c("sentrixid_mom","bmi_grs_mom"),to=c("sentrixid_dad","bmi_grs_dad"))
dat<-merge2(dat,bmi_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

eduyears_grs<-rename.vars(eduyears_grs,from=c("id","eduyears_grs"),to=c("sentrixid_mom","eduyears_grs_mom"))
dat<-merge2(dat,eduyears_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
eduyears_grs<-rename.vars(eduyears_grs,from=c("sentrixid_mom","eduyears_grs_mom"),to=c("sentrixid_dad","eduyears_grs_dad"))
dat<-merge2(dat,eduyears_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)


### GENETICALLY-PREDICTED LIKELIHOOD OF BINARY VARIABLES ###

dat<-subset2(dat,"dat$preg_plan==1")

gp_smkinit_mom<-NULL
gp_smkinit_dad<-NULL
gp_smkces_mom<-NULL

vars01<-c("bmi_mom","bmi_dad","bmi_mom")
vars02<-c("bmi_grs_mom","bmi_grs_dad","bmi_grs_mom")
vars03<-c("bmi_grs_mom_z","bmi_grs_dad_z","bmi_grs_mom_z")
vars04<-c("eduyears_mom","eduyears_dad","eduyears_mom")
vars05<-c("eduyears_grs_mom","eduyears_grs_dad","eduyears_grs_mom")
vars06<-c("eduyears_grs_mom_z","eduyears_grs_dad_z","eduyears_grs_mom_z")
vars07<-c("smkinit_mom","smkinit_dad","smkces_mom")
vars08<-c("smkinit_grs_mom","smkinit_grs_dad","smkces_grs_mom")
vars09<-c("smkinit_grs_mom_z","smkinit_grs_dad_z","smkces_grs_mom_z")
vars10<-c("smkinit_na_mom","smkinit_na_dad","smkces_na_mom")
vars11<-c("gp_smkinit_na_mom","gp_smkinit_na_dad","gp_smkces_na_mom")
vars12<-c("m_id_2374","f_id_2374","m_id_2374")

vars13<-c("gp_smkinit_bmi_mom","gp_smkinit_bmi_dad","gp_smkces_bmi_mom")
vars14<-c("gp_smkinit_bmi_mom_z","gp_smkinit_bmi_dad_z","gp_smkces_bmi_mom_z")
vars15<-c("gp_smkinit_eduyears_mom","gp_smkinit_eduyears_dad","gp_smkces_eduyears_mom")
vars16<-c("gp_smkinit_eduyears_mom_z","gp_smkinit_eduyears_dad_z","gp_smkces_eduyears_mom_z")
vars17<-c("gp_smkinit_mom","gp_smkinit_dad","gp_smkces_mom")
vars18<-c("gp_smkinit_mom_z","gp_smkinit_dad_z","gp_smkces_mom_z")

varsaa<-list(gp_smkinit_mom,gp_smkinit_dad,gp_smkces_mom)

for(i in 1:length(vars01))
  
{
  dat[,vars11[i]]<-with(dat,ifelse(is.na(dat[,vars01[i]]) | is.na(dat[,vars02[i]]) | 
                                     is.na(dat[,vars04[i]]) | is.na(dat[,vars05[i]]) | 
                                     is.na(dat[,vars07[i]]) | is.na(dat[,vars08[i]]) | 
                                     dat[,vars10[i]]==1,1,0))
  
  datx<-subset2(dat,"dat[,vars11[i]]==0")
  a<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars02[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars02[i]])))
  a_z<-as.numeric(with(datx,scale(a)))
  b<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars05[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars05[i]])))
  b_z<-as.numeric(with(datx,scale(b)))
  c<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars08[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars08[i]])))
  c_z<-as.numeric(with(datx,scale(c)))
  
  mod01<-lm_robust(datx[,vars01[i]]~a_z+c_z, 
                   data=datx, clusters=datx[,vars12[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope01<-as.numeric(summary(mod01)$coefficients[2,1])
  slope02<-as.numeric(summary(mod01)$coefficients[3,1])
  dat[,vars13[i]]<-interc+
    (as.numeric(with(dat,scale(dat[,vars02[i]])))*slope01)+
    (as.numeric(with(dat,scale(dat[,vars08[i]])))*slope02)
  dat[,vars14[i]]<-as.numeric(with(dat,scale(dat[,vars13[i]])))
  
  mod01<-lm_robust(datx[,vars04[i]]~b_z+c_z, 
                   data=datx, clusters=datx[,vars12[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope01<-as.numeric(summary(mod01)$coefficients[2,1])
  slope02<-as.numeric(summary(mod01)$coefficients[3,1])
  dat[,vars15[i]]<-interc+
    (as.numeric(with(dat,scale(dat[,vars05[i]])))*slope01)+
    (as.numeric(with(dat,scale(dat[,vars08[i]])))*slope02)
  dat[,vars16[i]]<-as.numeric(with(dat,scale(dat[,vars15[i]])))
  
  mod02<-glm(formula=as.factor(datx[,vars07[i]])~c_z
             +a_z+b_z,
             data=datx, family="binomial")
  datx[,vars17[i]]<-predict.glm(mod02, type="response")
  datx[,vars18[i]]<-as.numeric(with(datx,scale(datx[,vars17[i]])))
  varsaa[[i]]<-as.data.frame(cbind(datx$preg_id_2374,datx[,vars17[i]],datx[,vars18[i]]))
  colnames(varsaa[[i]])<-c("preg_id_2374",vars17[i],vars18[i])
}

names(varsaa)<-c("gp_smkinit_mom","gp_smkinit_dad","gp_smkces_mom")
gp_smkinit_mom<-varsaa$gp_smkinit_mom
gp_smkinit_dad<-varsaa$gp_smkinit_dad
gp_smkces_mom<-varsaa$gp_smkces_mom
dat[,c("gen_pred_smkinit_mom","gen_pred_smkinit_mom_z",
       "gen_pred_smkinit_dad","gen_pred_smkinit_dad_z",
       "gen_pred_smkces_mom","gen_pred_smkces_mom_z")]<-NULL
dat<-merge2(dat,gp_smkinit_mom,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,gp_smkinit_dad,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,gp_smkces_mom,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)


### GENETICALLY-PREDICTED VALUES OF CONTINUOUS PARAMETERS ###

vars01<-c("bmi_mom","bmi_dad","bmi_mom")
vars02<-c("bmi_grs_mom","bmi_grs_dad","bmi_grs_mom")
vars03<-c("bmi_grs_mom_z","bmi_grs_dad_z","bmi_grs_mom_z")
vars04<-c("eduyears_mom","eduyears_dad","eduyears_mom")
vars05<-c("eduyears_grs_mom","eduyears_grs_dad","eduyears_grs_mom")
vars06<-c("eduyears_grs_mom_z","eduyears_grs_dad_z","eduyears_grs_mom_z")
vars07<-c("cigday_mom","cigday_dad","agesmk_mom")
vars08<-c("cigday_grs_mom","cigday_grs_dad","agesmk_grs_mom")
vars09<-c("cigday_grs_mom_z","cigday_grs_dad_z","agesmk_grs_mom_z")
vars10<-c("cigday_na_mom","cigday_na_dad","agesmk_na_mom")
vars11<-c("gp_cigday_na_mom","gp_cigday_na_dad","gp_agesmk_na_mom")
vars12<-c("m_id_2374","f_id_2374","m_id_2374")

vars13<-c("gp_cigday_bmi_mom","gp_cigday_bmi_dad","gp_agesmk_bmi_mom")
vars14<-c("gp_cigday_bmi_mom_z","gp_cigday_bmi_dad_z","gp_agesmk_bmi_mom_z")
vars15<-c("gp_cigday_eduyears_mom","gp_cigday_eduyears_dad","gp_agesmk_eduyears_mom")
vars16<-c("gp_cigday_eduyears_mom_z","gp_cigday_eduyears_dad_z","gp_agesmk_eduyears_mom_z")
vars17<-c("gp_cigday_mom","gp_cigday_dad","gp_agesmk_mom")
vars18<-c("gp_cigday_mom_z","gp_cigday_dad_z","gp_agesmk_mom_z")

for(i in 1:length(vars01))
  
{
  dat[,vars11[i]]<-with(dat,ifelse(is.na(dat[,vars01[i]]) | is.na(dat[,vars02[i]]) | 
                                     is.na(dat[,vars04[i]]) | is.na(dat[,vars05[i]]) | 
                                     is.na(dat[,vars07[i]]) | is.na(dat[,vars08[i]]) | 
                                     dat[,vars10[i]]==1,1,0))
  
  datx<-subset2(dat,"dat[,vars11[i]]==0")
  a<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars02[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars02[i]])))
  a_z<-as.numeric(with(datx,scale(a)))
  b<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars05[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars05[i]])))
  b_z<-as.numeric(with(datx,scale(b)))
  c<-with(datx,ifelse(datx[,vars11[i]]==0,datx[,vars08[i]],
                      ifelse(datx[,vars11[i]]==1,NA,datx[,vars08[i]])))
  c_z<-as.numeric(with(datx,scale(c)))
  
  mod01<-lm_robust(datx[,vars01[i]]~a_z+c_z, 
                   data=datx, clusters=datx[,vars12[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope01<-as.numeric(summary(mod01)$coefficients[2,1])
  slope02<-as.numeric(summary(mod01)$coefficients[3,1])
  dat[,vars13[i]]<-interc+
    (as.numeric(with(dat,scale(dat[,vars02[i]])))*slope01)+
    (as.numeric(with(dat,scale(dat[,vars08[i]])))*slope02)
  dat[,vars14[i]]<-as.numeric(with(dat,scale(dat[,vars13[i]])))
  
  mod01<-lm_robust(datx[,vars04[i]]~b_z+c_z, 
                   data=datx, clusters=datx[,vars12[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope01<-as.numeric(summary(mod01)$coefficients[2,1])
  slope02<-as.numeric(summary(mod01)$coefficients[3,1])
  dat[,vars15[i]]<-interc+
    (as.numeric(with(dat,scale(dat[,vars05[i]])))*slope01)+
    (as.numeric(with(dat,scale(dat[,vars08[i]])))*slope02)
  dat[,vars16[i]]<-as.numeric(with(dat,scale(dat[,vars15[i]])))
  
  mod02<-lm_robust(datx[,vars04[i]]~c_z+a_z+b_z, 
                   data=datx, clusters=datx[,vars12[i]], se_type="stata")
  interc<-as.numeric(summary(mod02)$coefficients[1,1])
  slope01<-as.numeric(summary(mod02)$coefficients[2,1])
  slope02<-as.numeric(summary(mod02)$coefficients[3,1])
  slope03<-as.numeric(summary(mod02)$coefficients[4,1])
  dat[,vars17[i]]<-interc+
    (as.numeric(with(dat,scale(dat[,vars08[i]])))*slope01)+
    (as.numeric(with(dat,scale(dat[,vars02[i]])))*slope02)+
    (as.numeric(with(dat,scale(dat[,vars05[i]])))*slope03)
  dat[,vars18[i]]<-as.numeric(with(dat,scale(dat[,vars11[i]])))
}


### ROBUSTNESS OF GENETICALLY PREDICTED BMI AND EDUYEARS ###

vars01<-c("bmi_mom","bmi_dad","eduyears_mom","eduyears_dad")
vars02<-c("bmi_grs_mom","bmi_grs_dad","eduyears_grs_mom","eduyears_grs_dad")
vars03<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars04<-c("896","896","1159","1159")
vars05<-c("gp_smkinit_na_mom","gp_smkinit_na_dad","gp_smkinit_na_mom","gp_smkinit_na_dad")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars05[i]]==0")
  mod01<-lm_robust(datx[,vars01[i]]~datx[,vars02[i]], 
                   data=datx, clusters=datx[,vars03[i]], se_type="stata")
  coef<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
              guapa(summary(mod01)$coefficients[2,5]),"; ",
              guapa(summary(mod01)$coefficients[2,6]),"]",sep="")
  fstat<-round(mod01$fstatistic[1],0)
  r2<-paste("'",guapa(mod01$adj.r.squared*100),sep="")
  
  datx<-subset2(datx,"!is.na(datx[,vars01[i]]) & !is.na(datx[,vars02[i]])")
  n_obs<-length(which(!is.na(datx[,vars02[i]])))
  n_snps<-vars04[i]
  mean_grs<-paste(guapa(mean(datx[,vars02[i]]))," (",guapa(sd(datx[,vars02[i]])),")",sep="")
  tab<-rbind(tab,cbind(n_obs,n_snps,mean_grs,coef,fstat,r2))
}

rownames(tab)<-vars01
write.table(tab,file="./assoc_grs_bmi_eduyears.csv",sep=";",col.names=NA)


### LINEAR MR ###

vars01<-c("gp_smkinit_mom_z","gp_smkinit_dad_z","gp_agesmk_mom_z","gp_smkces_mom_z",
          "gp_cigday_mom_z","gp_cigday_dad_z")
vars02<-c("gp_smkinit_eduyears_mom_z","gp_smkinit_eduyears_dad_z","gp_agesmk_eduyears_mom_z","gp_smkces_eduyears_mom_z",
          "gp_cigday_eduyears_mom_z","gp_cigday_eduyears_dad_z")
vars03<-c("gp_smkinit_bmi_mom_z","gp_smkinit_bmi_dad_z","gp_agesmk_bmi_mom_z","gp_smkces_bmi_mom_z",
          "gp_cigday_bmi_mom_z","gp_cigday_bmi_dad_z")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_mom","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_mom","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_mom","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_mom","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_mom","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_mom","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_mom","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_mom","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_mom","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_mom","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","m_id_2374","m_id_2374","f_id_2374")
vars19<-c("gp_smkinit_na_mom","gp_smkinit_na_dad","gp_agesmk_na_mom","gp_smkces_na_mom",
          "gp_cigday_na_mom","gp_cigday_na_dad")
vars20<-c("smkinit_mom","smkinit_dad","agesmk_mom","smkces_mom","cigday_mom","cigday_dad")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0")
  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars02[i]]+datx[,vars03[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  tab<-rbind(tab,cbind(coef03,pval03))
}

colnames(tab)<-c("OR","pval")
rownames(tab)<-vars20
write.table(tab,file="./mvmr_linear.csv",sep=";",col.names=NA)


### AGE-STRATIFIED ANALYSES ###
###############################

### BELOW MEDIAN ###

load("N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")
dat<-subset2(dat,"dat$preg_plan==1")
dat$agedelivery_mom_m<-as.numeric(ntile(dat$agedelivery_mom, 2))
dat$agedelivery_dad_m<-as.numeric(ntile(dat$agedelivery_dad, 2))

vars00<-c("gen_pred_smkinit_mom","gen_pred_smkinit_dad",
          "gen_pred_agesmk_mom","gen_pred_smkces_mom",
          "gen_pred_cigday_mom","gen_pred_cigday_dad")
vars01<-c("gen_pred_smkinit_mom_z","gen_pred_smkinit_dad_z",
          "gen_pred_agesmk_mom_z","gen_pred_smkces_mom_z",
          "gen_pred_cigday_mom_z","gen_pred_cigday_dad_z")

for(i in 1:length(vars01))
{
  dat[,vars01[i]]<-as.numeric(with(dat,scale(dat[,vars00[i]])))
} 

vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_mom","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_mom","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_mom","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_mom","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_mom","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_mom","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_mom","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_mom","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_mom","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_mom","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","m_id_2374","m_id_2374","f_id_2374")
vars19<-c("smkinit_na_mom","smkinit_na_dad","agesmk_na_mom","smkces_na_mom","cigday_na_mom","cigday_na_dad")
vars20<-c("SmkInit, mothers","SmkInit, fathers",
          "AgeSmk, mothers","SmkCes, mothers",
          "CigDay, mothers","CigDay, fathers")
vars21<-c("agedelivery_mom_m","agedelivery_dad_m","agedelivery_mom_m",
          "agedelivery_mom_m","agedelivery_mom_m","agedelivery_dad_m")


tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars01[i]]) & dat[,vars21[i]]==1")
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-rbind(tab,cbind(coef04,pval04))
}

colnames(tab)<-c("GRS-OR (adj.)","GRS-pval (adj.)")
rownames(tab)<-vars20
write.table(tab,file="./age_below_linear.csv",sep=";",col.names=NA)


### ABOVE MEDIAN ###

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars01[i]]) & dat[,vars21[i]]==2")
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-rbind(tab,cbind(coef04,pval04))
}

colnames(tab)<-c("GRS-OR (adj.)","GRS-pval (adj.)")
rownames(tab)<-vars20
write.table(tab,file="./age_above_linear.csv",sep=";",col.names=NA)
