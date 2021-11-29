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


######################################
### GENERATION OF WORKING DATABASE ###
######################################

### MERGING MOM/DAD GRSs ###

load("./source_files/smkinit_grs.RData")
smkinit_grs<-dat
load("./source_files/agesmk_grs.RData")
agesmk_grs<-dat
load("./source_files/cigday_grs.RData")
cigday_grs<-dat
load("./source_files/smkces_grs.RData")
smkces_grs<-dat
load("./source_files/lifesmk_grs.RData")
lifesmk_grs<-dat

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_raw.RData")
dat$sentrixid_mom<-with(dat,ifelse(sentrixid_mom=="",NA,sentrixid_mom))
dat$sentrixid_dad<-with(dat,ifelse(sentrixid_dad=="",NA,sentrixid_dad))

smkinit_grs<-rename.vars(smkinit_grs,from=c("id","smkinit_grs"),to=c("sentrixid_mom","smkinit_grs_mom"))
dat<-merge2(dat,smkinit_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
smkinit_grs<-rename.vars(smkinit_grs,from=c("sentrixid_mom","smkinit_grs_mom"),to=c("sentrixid_dad","smkinit_grs_dad"))
dat<-merge2(dat,smkinit_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

agesmk_grs<-rename.vars(agesmk_grs,from=c("id","agesmk_grs"),to=c("sentrixid_mom","agesmk_grs_mom"))
dat<-merge2(dat,agesmk_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
agesmk_grs<-rename.vars(agesmk_grs,from=c("sentrixid_mom","agesmk_grs_mom"),to=c("sentrixid_dad","agesmk_grs_dad"))
dat<-merge2(dat,agesmk_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

cigday_grs<-rename.vars(cigday_grs,from=c("id","cigday_grs"),to=c("sentrixid_mom","cigday_grs_mom"))
dat<-merge2(dat,cigday_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
cigday_grs<-rename.vars(cigday_grs,from=c("sentrixid_mom","cigday_grs_mom"),to=c("sentrixid_dad","cigday_grs_dad"))
dat<-merge2(dat,cigday_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

smkces_grs<-rename.vars(smkces_grs,from=c("id","smkces_grs"),to=c("sentrixid_mom","smkces_grs_mom"))
dat<-merge2(dat,smkces_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
smkces_grs<-rename.vars(smkces_grs,from=c("sentrixid_mom","smkces_grs_mom"),to=c("sentrixid_dad","smkces_grs_dad"))
dat<-merge2(dat,smkces_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

lifesmk_grs<-rename.vars(lifesmk_grs,from=c("id","lifesmk_grs"),to=c("sentrixid_mom","lifesmk_grs_mom"))
dat<-merge2(dat,lifesmk_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
lifesmk_grs<-rename.vars(lifesmk_grs,from=c("sentrixid_mom","lifesmk_grs_mom"),to=c("sentrixid_dad","lifesmk_grs_dad"))
dat<-merge2(dat,lifesmk_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

smkinit_grs<-NULL
agesmk_grs<-NULL
cigday_grs<-NULL
smkces_grs<-NULL
lifesmk_grs<-NULL

length(which(!is.na(dat$smkinit_grs_mom)))
length(which(!is.na(dat$agesmk_grs_mom)))
length(which(!is.na(dat$cigday_grs_mom)))
length(which(!is.na(dat$smkces_grs_mom)))
length(which(!is.na(dat$lifesmk_grs_mom)))
length(which(!is.na(dat$smkinit_grs_dad)))
length(which(!is.na(dat$agesmk_grs_dad)))
length(which(!is.na(dat$cigday_grs_dad)))
length(which(!is.na(dat$smkces_grs_dad)))
length(which(!is.na(dat$lifesmk_grs_dad)))


### MERGING QC INDICATORS ###

qc_indic<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/Stata/Marker_QC_Ind_Alex_group.txt",
                                   header=TRUE,sep=","))
qc_indic$V2<-NULL
qc_indic<-rename.vars(qc_indic,from=c("V1"),to=c("sentrixid_mom"))
qc_indic$qc_gen_mom<-1
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
qc_indic<-rename.vars(qc_indic,from=c("sentrixid_mom","qc_gen_mom"),to=c("sentrixid_dad","qc_gen_dad"))
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
dat$qc_gen_mom<-with(dat,ifelse(is.na(qc_gen_mom),0,qc_gen_mom))
dat$qc_gen_dad<-with(dat,ifelse(is.na(qc_gen_dad),0,qc_gen_dad))

moms<-as.data.frame(read_dta("N:/data/durable/RAW/MoBaGenomics/_key/2019_10_01_MoBaGenetics_mother_2374.dta"))
names(moms)<-tolower(names(moms))
moms<-rename.vars(moms,from=c("sentrixid_m","batch_m"),to=c("sentrixid_mom","batch_mom"))
moms<-moms[,c("sentrixid_mom","batch_mom")]
dat<-merge2(dat,moms,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
dads<-as.data.frame(read_dta("N:/data/durable/RAW/MoBaGenomics/_key/2019_10_01_MoBaGenetics_father_2374.dta"))
names(dads)<-tolower(names(dads))
dads<-rename.vars(dads,from=c("sentrixid_f","batch_f"),to=c("sentrixid_dad","batch_dad"))
dads<-dads[,c("sentrixid_dad","batch_dad")]
dat<-merge2(dat,dads,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(dat$batch_mom)))
length(which(!is.na(dat$batch_dad)))


### BMI AND HEIGHT MEASURED VALUES ###

dat$bmi_mom<-dat$aa85/((dat$aa87/100)^2)
dat$bmi_dad<-dat$aa89/((dat$aa88/100)^2)

dat$bmi_mom<-with(dat,ifelse(bmi_mom<13 | bmi_mom>60,NA,bmi_mom))
dat$bmi_dad<-with(dat,ifelse(bmi_dad<13 | bmi_dad>60,NA,bmi_dad))


### HIGHEST EDUCATIONAL LEVEL, COMPLETED OR ONGOING ###

# Mothers #

dat$aa1125<-with(dat,ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125>aa1124),aa1125,
                            ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125<aa1124),NA,aa1125)))
dat$aa1128x<-with(dat,ifelse(aa1128==1 | aa1129==1,1,0))
dat$aa1128x<-with(dat,ifelse(is.na(aa1128x),0,aa1128x))

dat$edu_mom<-with(dat,ifelse(aa1124==1,1,
                             ifelse(aa1124==2,2,
                                    ifelse(aa1124==3,3,
                                           ifelse(aa1124==4,4,
                                                  ifelse(aa1124==5,5,
                                                         ifelse(aa1124==6,6,NA)))))))
dat$edu_mom<-with(dat,ifelse(is.na(aa1125),edu_mom,
                             ifelse(aa1125==1,1,
                                    ifelse(aa1125==2,2,
                                           ifelse(aa1125==3,3,
                                                  ifelse(aa1125==4,4,
                                                         ifelse(aa1125==5,5,
                                                                ifelse(aa1125==6,6,NA))))))))
dat$edu_mom<-with(dat,ifelse(is.na(edu_mom),0,edu_mom))

dat$edu_mom<-with(dat,ifelse(edu_mom==0 & aa1128x==1,3,
                             ifelse((edu_mom>0 & edu_mom<3) & aa1128x==1,3,
                                    ifelse(edu_mom>=3 & aa1128x==1,edu_mom,
                                           ifelse(edu_mom==0 & aa1128x==0,NA,
                                                  ifelse((edu_mom>0 & edu_mom<3) & aa1128x==0,edu_mom,
                                                         ifelse(edu_mom>=3 & aa1128x==0,edu_mom,NA)))))))

dat$eduyears_mom<-with(dat,ifelse(edu_mom==1,10,
                                  ifelse(edu_mom==2,10,
                                         ifelse(edu_mom==3,13,
                                                ifelse(edu_mom==4,13,
                                                       ifelse(edu_mom==5,19,
                                                              ifelse(edu_mom==6,20,
                                                                     ifelse(is.na(edu_mom),NA,NA))))))))

dat$edu_mom<-with(dat,ifelse(edu_mom==1,1,
                             ifelse(edu_mom==2,1,
                                    ifelse(edu_mom==3,2,
                                           ifelse(edu_mom==4,2,
                                                  ifelse(edu_mom==5,3,
                                                         ifelse(edu_mom==6,4,
                                                                ifelse(is.na(edu_mom),NA,NA))))))))

# Fathers #

dat$aa1127<-with(dat,ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127>aa1126),aa1127,
                            ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127<aa1126),NA,aa1127)))
dat$aa1130x<-with(dat,ifelse(aa1130==1 | aa1131==1,1,0))
dat$aa1130x<-with(dat,ifelse(is.na(aa1130x),0,aa1130x))

dat$edu_dad<-with(dat,ifelse(aa1126==1,1,
                             ifelse(aa1126==2,2,
                                    ifelse(aa1126==3,3,
                                           ifelse(aa1126==4,4,
                                                  ifelse(aa1126==5,5,
                                                         ifelse(aa1126==6,6,NA)))))))
dat$edu_dad<-with(dat,ifelse(is.na(aa1127),edu_dad,
                             ifelse(aa1127==1,1,
                                    ifelse(aa1127==2,2,
                                           ifelse(aa1127==3,3,
                                                  ifelse(aa1127==4,4,
                                                         ifelse(aa1127==5,5,
                                                                ifelse(aa1127==6,6,NA))))))))
dat$edu_dad<-with(dat,ifelse(is.na(edu_dad),0,edu_dad))

dat$edu_dad<-with(dat,ifelse(edu_dad==0 & aa1130x==1,3,
                             ifelse((edu_dad>0 & edu_dad<3) & aa1130x==1,3,
                                    ifelse(edu_dad>=3 & aa1130x==1,edu_dad,
                                           ifelse(edu_dad==0 & aa1130x==0,NA,
                                                  ifelse((edu_dad>0 & edu_dad<3) & aa1130x==0,edu_dad,
                                                         ifelse(edu_dad>=3 & aa1130x==0,edu_dad,NA)))))))

dat$eduyears_dad<-with(dat,ifelse(edu_dad==1,10,
                                  ifelse(edu_dad==2,10,
                                         ifelse(edu_dad==3,13,
                                                ifelse(edu_dad==4,13,
                                                       ifelse(edu_dad==5,19,
                                                              ifelse(edu_dad==6,20,
                                                                     ifelse(is.na(edu_dad),NA,NA))))))))
dat$edu_dad<-with(dat,ifelse(edu_dad==1,1,
                             ifelse(edu_dad==2,1,
                                    ifelse(edu_dad==3,2,
                                           ifelse(edu_dad==4,2,
                                                  ifelse(edu_dad==5,3,
                                                         ifelse(edu_dad==6,4,
                                                                ifelse(is.na(edu_dad),NA,NA))))))))


### DEFINITION OF PARITY (number of previous deliveries) ###

dat$parity<-with(dat,ifelse(paritet_5==0,0,
                            ifelse(paritet_5==1,1,
                                   ifelse(paritet_5==2,1,
                                          ifelse(paritet_5==3,1,
                                                 ifelse(paritet_5==4,1,NA))))))
dat$parity_cont<-dat$paritet_5

length(which(is.na(dat$parity)))
# 486 NAs in "paritet_5" 


### DEFINITION OF EVER SMOKERS (smkinit) ###

# Mothers #

dat$smkinit_mom<-with(dat,ifelse(aa1355==2 | aa1356>1 | aa1357!=0 | aa1358!=0 | aa1359>1 | aa1360!=0 | aa1361!=0 | 
                                   !is.na(aa1362) | aa1363==2 | !is.na(aa1364) | !is.na(aa1365),1,0))
dat$smkinit_mom<-with(dat,ifelse(is.na(smkinit_mom),0,smkinit_mom))
dat$smkinit_mom<-with(dat,ifelse(is.na(aa1355) & is.na(aa1356) & is.na(aa1357) & is.na(aa1358) & is.na(aa1359) & is.na(aa1360) & is.na(aa1361) & 
                                   is.na(aa1362) & is.na(aa1363) & is.na(aa1364) & is.na(aa1365),NA,smkinit_mom))

# Fathers #

dat$smkinit_dad<-with(dat,ifelse(aa1353==2 | aa1354==2 | ff214==2 | ff215>1 | ff216!=0 | ff217!=0 | ff218>1 | ff219!=0 | ff220!=0,1,0))
dat$smkinit_dad<-with(dat,ifelse(is.na(smkinit_dad),0,smkinit_dad))
dat$smkinit_dad<-with(dat,ifelse(is.na(aa1353) & is.na(aa1354) & is.na(ff214) & is.na(ff215) & is.na(ff216) & is.na(ff217) & is.na(ff218) & is.na(ff219) & is.na(ff220),NA,smkinit_dad))

dat$smkinit_na_mom<-with(dat,ifelse(is.na(smkinit_mom),1,0))
dat$smkinit_na_dad<-with(dat,ifelse(is.na(smkinit_dad),1,0))

dat$dadnumber<-with(dat,ifelse(!is.na(ff214) | !is.na(ff215) | !is.na(ff216) | !is.na(ff217) | !is.na(ff218) | !is.na(ff219) | !is.na(ff220),1,0))


### DEFINITION OF CigDay ###

# Mothers #

dat$xxx<-with(dat,ifelse(!is.na(aa1360) & is.na(aa1361),aa1360,
                         ifelse(is.na(aa1360) & !is.na(aa1361),aa1361*7,
                                ifelse(!is.na(aa1360) & !is.na(aa1361) & (aa1360/7)>aa1361,aa1360,
                                       ifelse(!is.na(aa1360) & !is.na(aa1361) & (aa1360/7)<aa1361,aa1361*7,
                                              ifelse(!is.na(aa1360) & !is.na(aa1361) & (aa1360/7)==aa1361,aa1360,NA))))))
dat$yyy<-with(dat,ifelse(!is.na(aa1357) & is.na(aa1358),aa1357,
                         ifelse(is.na(aa1357) & !is.na(aa1358),aa1358*7,
                                ifelse(!is.na(aa1357) & !is.na(aa1358) & (aa1357/7)>aa1358,aa1357,
                                       ifelse(!is.na(aa1357) & !is.na(aa1358) & (aa1357/7)<aa1358,aa1358*7,
                                              ifelse(!is.na(aa1357) & !is.na(aa1358) & (aa1357/7)==aa1358,aa1357,NA))))))
dat$cigday_mom<-with(dat,ifelse(!is.na(xxx) & is.na(yyy),xxx,
                                ifelse(is.na(xxx) & !is.na(yyy),yyy,
                                       ifelse(!is.na(xxx) & !is.na(yyy) & xxx>yyy,xxx,
                                              ifelse(!is.na(xxx) & !is.na(yyy) & xxx<yyy,yyy,
                                                     ifelse(!is.na(xxx) & !is.na(yyy) & xxx==yyy,xxx,NA))))))
dat$cigday_mom<-with(dat,ifelse(cigday_mom==0,NA,cigday_mom))

dat$cigday_cat_mom<-with(dat, ifelse(cigday_mom<7,1,
                                     ifelse(cigday_mom<=35,2,
                                            ifelse(cigday_mom<=70,3,
                                                   ifelse(cigday_mom<=105,4,
                                                          ifelse(cigday_mom<=140,5,
                                                                 ifelse(cigday_mom>140,6,NA)))))))


# Fathers #

dat$xxx<-with(dat,ifelse(!is.na(ff219) & is.na(ff220),ff219,
                         ifelse(is.na(ff219) & !is.na(ff220),ff220*7,
                                ifelse(!is.na(ff219) & !is.na(ff220) & (ff219/7)>ff220,ff219,
                                       ifelse(!is.na(ff219) & !is.na(ff220) & (ff219/7)<ff220,ff220*7,
                                              ifelse(!is.na(ff219) & !is.na(ff220) & (ff219/7)==ff220,ff219,NA))))))
dat$yyy<-with(dat,ifelse(!is.na(ff216) & is.na(ff217),ff216,
                         ifelse(is.na(ff216) & !is.na(ff217),ff217*7,
                                ifelse(!is.na(ff216) & !is.na(ff217) & (ff216/7)>ff217,ff216,
                                       ifelse(!is.na(ff216) & !is.na(ff217) & (ff216/7)<ff217,ff217*7,
                                              ifelse(!is.na(ff216) & !is.na(ff217) & (ff216/7)==ff217,ff216,NA))))))
dat$cigday_dad<-with(dat,ifelse(!is.na(xxx) & is.na(yyy),xxx,
                                ifelse(is.na(xxx) & !is.na(yyy),yyy,
                                       ifelse(!is.na(xxx) & !is.na(yyy) & xxx>yyy,xxx,
                                              ifelse(!is.na(xxx) & !is.na(yyy) & xxx<yyy,yyy,
                                                     ifelse(!is.na(xxx) & !is.na(yyy) & xxx==yyy,xxx,NA))))))
dat$cigday_dad<-with(dat,ifelse(cigday_dad==0,NA,cigday_dad))

dat$cigday_cat_dad<-with(dat, ifelse(cigday_dad<7,1,
                                     ifelse(cigday_dad<=35,2,
                                            ifelse(cigday_dad<=70,3,
                                                   ifelse(cigday_dad<=105,4,
                                                          ifelse(cigday_dad<=140,5,
                                                                 ifelse(cigday_dad>140,6,NA)))))))


### DEFINITION OF CURRENT SMOKERS (CURRENT + QUITTERS <=2 YEARS BEFORE) ###

dat$currentsmk_mom<-with(dat,ifelse(smoking_mom>1,1,0))
dat$xxx<-dat$mors_alder-dat$aa1364
dat$xxx<-with(dat,ifelse(xxx<(-5),NA,
                         ifelse(xxx<(0),0,xxx)))
dat$currentsmk_mom<-with(dat,ifelse(!is.na(currentsmk_mom) & !is.na(xxx) & currentsmk_mom==0 & xxx<=2,1,currentsmk_mom))

dat$currentsmk_dad<-with(dat,ifelse(smoking_dad>1,1,0))

dat$cigday_na_mom<-with(dat,ifelse(currentsmk_mom==1,0,1))
dat$cigday_na_mom<-with(dat,ifelse(is.na(currentsmk_mom),1,cigday_na_mom))
dat$cigday_na_mom<-with(dat,ifelse(is.na(cigday_mom),1,cigday_na_mom))

dat$cigday_na_dad<-with(dat,ifelse(currentsmk_dad==1,0,1))
dat$cigday_na_dad<-with(dat,ifelse(is.na(currentsmk_dad),1,cigday_na_dad))
dat$cigday_na_dad<-with(dat,ifelse(is.na(cigday_dad),1,cigday_na_dad))


### DEFINITION OF AgeSmkInit and AgeSmkCes IN MOTHERS (NO DATA IN FATHERS) ###

dat$agesmkinit_mom<-dat$aa1362
dat$agesmk_mom<-dat$aa1362
dat$agesmkces_mom<-dat$aa1364


### DEFINITION OF SmkCes IN MOTHERS (NO DATA IN FATHERS) ###

dat$smkces_mom<-with(dat,ifelse(aa1363==2 & is.na(aa1365),1,0))
dat$smkces_mom<-with(dat,ifelse(is.na(smkces_mom),0,smkces_mom))
dat$smkces_mom<-with(dat,ifelse(is.na(smoking_mom),NA,smkces_mom))


### MERGING PRINCIPAL COMPONENTS ###

pcs<-as.data.frame(read.csv("N:/data/durable/Projects/Magnus_MR_BMI/dat/chr_tmp/combined_filtered_indep_indep.csv",
                            header=TRUE,sep=",",dec=".")[,3:13])
names(pcs)<-c("sentrixid_mom","pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom",
              "pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom")
pcs$sentrixid_mom<-as.character(pcs$sentrixid_mom)
dat<-merge2(dat,pcs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)

pcs<-rename.vars(pcs,
                 from=c("sentrixid_mom","pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom",
                        "pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom"),
                 to=c("sentrixid_dad","pc1_dad","pc2_dad","pc3_dad","pc4_dad","pc5_dad",
                      "pc6_dad","pc7_dad","pc8_dad","pc9_dad","pc10_dad"))
pcs$sentrixid_dad<-as.character(pcs$sentrixid_dad)
dat<-merge2(dat,pcs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(dat$pc1_mom)))
length(which(!is.na(dat$pc2_mom)))

length(which(!is.na(dat$pc1_dad)))
length(which(!is.na(dat$pc2_dad)))


### DEFINITION OF SUBFERTILITY ###

dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$art_nona<-with(dat,ifelse(!is.na(art) & art>0,1,0))

dat$subfertility_12plus<-with(dat,ifelse(aa48<12,0,
                                         ifelse(aa48>=12,1,NA)))
dat$subfertility_12plus<-with(dat,ifelse(art_nona==0,subfertility_12plus,
                                         ifelse(art_nona==1,1,NA)))
dat$subfertility_12plus<-with(dat,ifelse(exclude==0 & is.na(subfertility_12plus),0,
                                         ifelse(exclude==0 & subfertility_12plus==0,0,
                                                ifelse(exclude==0 & subfertility_12plus==1,1,
                                                       ifelse(exclude==1 & is.na(subfertility_12plus),NA,
                                                              ifelse(exclude==1 & subfertility_12plus==0,NA,
                                                                     ifelse(exclude==1 & subfertility_12plus==1,NA,NA)))))))

length(which(!is.na(dat$subfertility)))
table(dat$subfertility)

dat$preg_plan<-with(dat,ifelse(is.na(aa46),9,
                               ifelse(aa46==0,9,
                                      ifelse(aa46==1,0,
                                             ifelse(aa46==2,1,9)))))
attributes(dat$preg_plan)$value.label<-c("0=No","1=Yes","9=NA")
attributes(dat$preg_plan)$label<-c("Planned pregnancy")


### FINAL FORMAT OF DATA ###

dat<-rename.vars(dat,
                 from=c("mors_alder","fars_alder"),
                 to=c("agedelivery_mom","agedelivery_dad"))

dat$exclude<-with(dat,ifelse(is.na(subfertility),1,
                             ifelse(flerfodsel==1,1,0)))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
dat$flerfodsel<-with(dat,ifelse(is.na(flerfodsel),0,flerfodsel))

dat$exclude_batch_mom<-with(dat,ifelse(batch_mom=="TED",1,
                                       ifelse(batch_mom=="NORMENT-JAN15",1,
                                              ifelse(batch_mom=="NORMENT-JUN15",1,0))))
dat$exclude_batch_mom<-with(dat,ifelse(is.na(exclude_batch_mom),0,exclude_batch_mom))

dat$exclude_batch_dad<-with(dat,ifelse(batch_dad=="TED",1,
                                       ifelse(batch_dad=="NORMENT-JAN15",1,
                                              ifelse(batch_dad=="NORMENT-JUN15",1,0))))
dat$exclude_batch_dad<-with(dat,ifelse(is.na(exclude_batch_dad),0,exclude_batch_dad))

dat$exclude_skjema_tbl1<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))

dat$exclude_subfertility_na<-with(dat,ifelse(is.na(subfertility),1,0))

dat$exclude_flerfodsel<-with(dat,ifelse(flerfodsel==1,1,0))

nam_ok<-c("sentrixid_mom","sentrixid_dad","m_id_2374","f_id_2374","preg_id_2374","flerfodsel","preg_plan","art_nona",
          "smoking_mom","smoking_dad","smkinit_mom","smkinit_dad","smkinit_na_mom","smkinit_na_dad",
          "cigday_mom","cigday_dad","cigday_cat_mom","cigday_cat_dad","cigday_na_mom","cigday_na_dad",
          "agesmk_mom","agesmkinit_mom","smkces_mom","agesmkces_mom",
          "agedelivery_mom","agedelivery_dad","bmi_mom","bmi_dad","eduyears_mom","eduyears_dad",
          "parity","parity_cont","subfertility","subfertility_12plus","batch_mom","batch_dad",
          "exclude","exclude_batch_mom","exclude_batch_dad","exclude_skjema_tbl1","exclude_flerfodsel","exclude_subfertility_na",
          "smkinit_grs_mom","smkinit_grs_dad","agesmk_grs_mom","agesmk_grs_dad",
          "cigday_grs_mom","cigday_grs_dad","smkces_grs_mom","smkces_grs_dad",
          "pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom","pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom",
          "pc1_dad","pc2_dad","pc3_dad","pc4_dad","pc5_dad","pc6_dad","pc7_dad","pc8_dad","pc9_dad","pc10_dad")
dat<-dat[,nam_ok]


######################################
### DEFINITION OF STUDY FLOW CHART ###
######################################

# Total registers
length(which(!is.na(dat$preg_id_2374)))

# Singleton pregnancies (flerfodsel==0)
dat<-subset2(dat,"dat$exclude_flerfodsel==0")
dim(dat)[1]

# Exclusion by version_skjema_tbl1
table(dat$exclude_skjema_tbl1)
table(dat$exclude_skjema_tbl1,by=dat$subfertility_12plus)
dat<-subset2(dat,"dat$exclude_skjema_tbl1==0")
dim(dat)[1]

# Exclusion of mothers/fathers without any information on smoking
dat$smoking_na_mom<-with(dat,ifelse(is.na(smoking_mom) & is.na(smkinit_mom) & is.na(cigday_mom),1,0))
table(dat$smoking_na_mom)

dat$smoking_na_dad<-with(dat,ifelse(is.na(smoking_dad) & is.na(smkinit_dad) & is.na(cigday_dad),1,0))
table(dat$smoking_na_dad)

table(dat$dadnumber)

# Unique IDs of mothers/fathers with smoking information
length(unique(dat[which(dat$smoking_na_mom==0),c("m_id_2374")]))
length(unique(dat[which(dat$smoking_na_dad==0),c("f_id_2374")]))


# IDs of mothers/fathers with genotype data and smoking information
length(dat[which((!is.na(dat$smkinit_grs_mom) | !is.na(dat$cigday_grs_mom) | !is.na(dat$agesmk_grs_mom) | !is.na(dat$smkces_grs_mom)) 
                 & dat$smoking_na_mom==0 & !is.na(dat$sentrixid_mom)),c("sentrixid_mom")])
length(dat[which((!is.na(dat$smkinit_grs_dad) | !is.na(dat$cigday_dad)) 
                 & dat$smoking_na_dad==0 & !is.na(dat$sentrixid_dad)),c("sentrixid_dad")])


# Exclusion by batch_m / batch_d or lack of 
table(dat$exclude_batch_mom)
table(dat$exclude_batch_dad)
dat$exclude_genotype_mom<-with(dat,ifelse(exclude_batch_mom==1 | is.na(dat$sentrixid_mom),1,0))
dat$exclude_genotype_dad<-with(dat,ifelse(exclude_batch_dad==1 | is.na(dat$sentrixid_dad),1,0))

length(dat[which((!is.na(dat$smkinit_grs_mom) | !is.na(dat$cigday_grs_mom) | !is.na(dat$agesmk_grs_mom) | !is.na(dat$smkces_grs_mom)) 
                 & dat$smoking_na_mom==0 & dat$exclude_genotype_mom==0 & !is.na(dat$sentrixid_mom)),c("sentrixid_mom")])
length(dat[which((!is.na(dat$smkinit_grs_dad) | !is.na(dat$cigday_dad)) 
                 & dat$smoking_na_dad==0 & dat$exclude_genotype_dad==0 & !is.na(dat$sentrixid_dad)),c("sentrixid_dad")])


# Genotyped mothers and fathers with anthropometric measurements & quality OK
length(unique(dat[which((!is.na(dat$smkinit_grs_mom) | !is.na(dat$cigday_grs_mom) | !is.na(dat$agesmk_grs_mom) | !is.na(dat$smkces_grs_mom)) 
                        & dat$smoking_na_mom==0 & !is.na(dat$sentrixid_mom) & dat$exclude_genotype_mom==0),c("sentrixid_mom")]))
length(unique(dat[which((!is.na(dat$smkinit_grs_dad) | !is.na(dat$cigday_dad)) 
                        & dat$smoking_na_dad==0 & !is.na(dat$sentrixid_dad) & dat$exclude_genotype_dad==0),c("sentrixid_dad")]))

table(dat$exclude_subfertility_na)


# Mothers/fathers without covariates
dat$edu_na_mom<-with(dat,ifelse(is.na(eduyears_mom),1,0))
table(dat$edu_na_mom)
dat$edu_na_dad<-with(dat,ifelse(is.na(eduyears_dad),1,0))
table(dat$edu_na_dad)

dat$bmi_na_mom<-with(dat,ifelse(is.na(bmi_mom),1,0))
table(dat$bmi_na_mom)
dat$bmi_na_dad<-with(dat,ifelse(is.na(bmi_dad),1,0))
table(dat$bmi_na_dad)


### EXCLUSION CIRTERIA ACCORDING TO SMOKING INFORMATION AND GRSs ###

dat$smkinit_na_mom<-with(dat,ifelse(smkinit_na_mom==1 | is.na(smkinit_grs_mom) | exclude_genotype_mom==1,1,smkinit_na_mom))
dat$smkinit_na_dad<-with(dat,ifelse(smkinit_na_dad==1 | is.na(smkinit_grs_dad) | exclude_genotype_dad==1,1,smkinit_na_dad))
dat$cigday_na_mom<-with(dat,ifelse(cigday_na_mom==1 | is.na(cigday_grs_mom) | exclude_genotype_mom==1,1,cigday_na_mom))
dat$cigday_na_dad<-with(dat,ifelse(cigday_na_dad==1 | is.na(cigday_grs_dad) | exclude_genotype_dad==1,1,cigday_na_dad))
dat$agesmk_na_mom<-with(dat,ifelse(is.na(agesmk_mom) | smoking_na_mom==1 | smkinit_mom==0 | is.na(agesmk_grs_mom) | exclude_genotype_mom==1,1,0))
dat$smkces_na_mom<-with(dat,ifelse(is.na(smkces_mom) | smoking_na_mom==1 | smkinit_mom==0 | is.na(smkces_grs_mom) | exclude_genotype_mom==1,1,0))


#######################################################
### CALCULATION OF GENETICALLY-DETERMINED VARIABLES ###
#######################################################

# Calculation of genetically-predetermined continuous variables

vars01<-c("cigday_mom","cigday_dad","agesmk_mom")
vars02<-c("cigday_grs_mom","cigday_grs_dad","agesmk_grs_mom")
vars03<-c("cigday_mom_z","cigday_dad_z","agesmk_mom_z")
vars04<-c("cigday_grs_mom_z","cigday_grs_dad_z","agesmk_grs_mom_z")
vars05<-c("cigday_na_mom","cigday_na_dad","agesmk_na_mom")
vars09<-c("m_id_2374","f_id_2374","m_id_2374")
vars10<-c("gen_pred_cigday_mom","gen_pred_cigday_dad","gen_pred_agesmk_mom")
vars11<-c("environ_cigday_mom","environ_cigday_dad","environ_agesmk_mom")
vars12<-c("gen_pred_cigday_mom_z","gen_pred_cigday_dad_z","gen_pred_agesmk_mom_z")
vars13<-c("environ_cigday_mom_z","environ_cigday_dad_z","environ_agesmk_mom_z")

for(i in 1:length(vars01))
  
{
  dat[,vars03[i]]<-as.numeric(with(dat,scale(dat[,vars01[i]])))
  dat[,vars04[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  
  datx<-subset2(dat,"!is.na(dat[,vars02[i]]) & dat[,vars05[i]]==0")
  mod01<-lm_robust(datx[,vars01[i]]~datx[,vars04[i]], 
                   data=datx, clusters=datx[,vars09[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope01<-as.numeric(summary(mod01)$coefficients[2,1])
  dat[,vars10[i]]<-interc+(dat[,vars04[i]]*slope01)
  dat[,vars10[i]]<-with(dat,ifelse(dat[,vars05[i]]==0,dat[,vars10[i]],
                                   ifelse(dat[,vars05[i]]==1,NA,dat[,vars10[i]])))
  dat[,vars11[i]]<-dat[,vars01[i]]-dat[,vars10[i]]
  dat[,vars11[i]]<-with(dat,ifelse(dat[,vars05[i]]==0,dat[,vars11[i]],
                                   ifelse(dat[,vars05[i]]==1,NA,dat[,vars11[i]])))
  dat[,vars12[i]]<-as.numeric(with(dat,scale(dat[,vars10[i]])))
  dat[,vars13[i]]<-as.numeric(with(dat,scale(dat[,vars11[i]])))
}


# Calculation of genetically-predicted likelihood of binary outcomes

smkinit_mom<-NULL
smkinit_dad<-NULL
smkces_mom<-NULL

vars01<-c("smkinit_na_mom","smkinit_na_dad","smkces_na_mom")
vars05<-c("smkinit_mom","smkinit_dad","smkces_mom")
vars06<-c("smkinit_grs_mom","smkinit_grs_dad","smkces_grs_mom")
vars08<-c("smkinit_grs_mom_z","smkinit_grs_dad_z","smkces_grs_mom_z")
vars09<-c("m_id_2374","f_id_2374","m_id_2374")
vars14<-c("gen_pred_smkinit_mom","gen_pred_smkinit_dad","gen_pred_smkces_mom")
vars15<-c("environ_smkinit_mom","environ_smkinit_dad","environ_smkces_mom")
vars16<-c("gen_pred_smkinit_mom_z","gen_pred_smkinit_dad_z","gen_pred_smkces_mom_z")
vars17<-c("environ_smkinit_mom_z","environ_smkinit_dad_z","environ_smkces_mom_z")
varsaa<-list(smkinit_mom,smkinit_dad,smkces_mom)

for(i in 1:length(vars01))
  
{
  dat[,vars08[i]]<-as.numeric(with(dat,scale(dat[,vars06[i]])))
  
  datx<-subset2(dat,"dat[,vars01[i]]==0 & !is.na(dat[,vars06[i]])")
  mod02<-glm(formula=as.factor(datx[,vars05[i]])~datx[,vars08[i]],
             data=datx, family="binomial")
  datx[,vars14[i]]<-predict.glm(mod02, type="response")
  datx[,vars14[i]]<-with(datx,ifelse(datx[,vars01[i]]==0,datx[,vars14[i]],
                                     ifelse(datx[,vars01[i]]==1,NA,datx[,vars14[i]])))
  datx[,vars15[i]]<-datx[,vars05[i]]-datx[,vars14[i]]
  datx[,vars15[i]]<-with(datx,ifelse(datx[,vars01[i]]==0,datx[,vars15[i]],
                                     ifelse(datx[,vars01[i]]==1,NA,datx[,vars15[i]])))
  datx[,vars16[i]]<-as.numeric(with(datx,scale(datx[,vars14[i]])))
  datx[,vars17[i]]<-as.numeric(with(datx,scale(datx[,vars15[i]])))
  varsaa[[i]]<-as.data.frame(cbind(datx$preg_id_2374,datx[,vars14[i]],datx[,vars15[i]],datx[,vars16[i]],datx[,vars17[i]]))
  colnames(varsaa[[i]])<-c("preg_id_2374",vars14[i],vars15[i],vars16[i],vars17[i])
}

names(varsaa)<-c("smkinit_mom","smkinit_dad","smkces_mom")
smkinit_mom<-varsaa$smkinit_mom
smkinit_dad<-varsaa$smkinit_dad
smkces_mom<-varsaa$smkces_mom
dat<-merge2(dat,smkinit_mom,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,smkinit_dad,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,smkces_mom,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)


attributes(dat$agedelivery_mom)$label<-c("agedelivery_mom")
attributes(dat$agedelivery_dad)$label<-c("agedelivery_dad")
attributes(dat$bmi_mom)$label<-c("bmi_mom")
attributes(dat$bmi_dad)$label<-c("bmi_dad")
attributes(dat$agesmk_mom)$label<-c("agesmk_mom")

save(dat,file="N:/data/durable/Projects/Magnus_MR_smoking/MoBa_smoking.RData")