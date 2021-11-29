rm(list=ls())

# Load required libraries #
library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


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



##############################################
### GENERATE WORKING FILES FOR 2-SAMPLE MR ###
##############################################

setwd("N:/data/durable/Projects/Magnus_MR_smoking/2smr")

# EXPOSURE : SmkInit #

smkinit<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/2smr/gwas_smkinit_liu2019.csv",
                   header=TRUE,sep=";",dec=".")
smkinit$Phenotype<-c("SmkInit")
smkinit$units<-c("units")
smkinit$id<-c("Liu 2019")
smkinit<-rename.vars(smkinit,
                     from=c("chr","rsid","effect_allele_freq","n"),
                     to=c("chrom","SNP","eaf","samplesize"))
smkinit<-smkinit[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
                    "units","samplesize","id")]
write.table(smkinit,"./smkinit.txt",sep="\t",row.names=FALSE, quote=FALSE)

smkinit_harm<-read_exposure_data(filename = "./smkinit.txt",
                                 clump = FALSE,
                                 sep="\t",
                                 phenotype_col = "Phenotype",
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 eaf_col = "eaf",
                                 effect_allele_col = "effect_allele",
                                 other_allele_col = "other_allele",
                                 pval_col = "pval",
                                 units_col = "units",
                                 samplesize_col = "samplesize",
                                 id_col = "id",
                                 min_pval = 1e-200,
                                 log_pval = FALSE)
save(smkinit_harm,file="./smkinit_harm.RData")


# EXPOSURE : CigDay #

cigday<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/2smr/gwas_cigday_liu2019.csv",
                  header=TRUE,sep=";",dec=".")
cigday$Phenotype<-c("cigday")
cigday$units<-c("units")
cigday$id<-c("Liu 2019")
cigday<-rename.vars(cigday,
                    from=c("chr","rsid","effect_allele_freq","n"),
                    to=c("chrom","SNP","eaf","samplesize"))
cigday<-cigday[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
                  "units","samplesize","id")]
write.table(cigday,"./cigday.txt",sep="\t",row.names=FALSE, quote=FALSE)

cigday_harm<-read_exposure_data(filename = "./cigday.txt",
                                clump = FALSE,
                                sep="\t",
                                phenotype_col = "Phenotype",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                eaf_col = "eaf",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                units_col = "units",
                                samplesize_col = "samplesize",
                                id_col = "id",
                                min_pval = 1e-200,
                                log_pval = FALSE)
save(cigday_harm,file="./cigday_harm.RData")


# EXPOSURE : SmkCes #

smkces<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/2smr/gwas_smkces_liu2019.csv",
                  header=TRUE,sep=";",dec=".")
smkces$Phenotype<-c("smkces")
smkces$units<-c("units")
smkces$id<-c("Liu 2019")
smkces<-rename.vars(smkces,
                    from=c("chr","rsid","effect_allele_freq","n"),
                    to=c("chrom","SNP","eaf","samplesize"))
smkces<-smkces[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
                  "units","samplesize","id")]
write.table(smkces,"./smkces.txt",sep="\t",row.names=FALSE, quote=FALSE)

smkces_harm<-read_exposure_data(filename = "./smkces.txt",
                                clump = FALSE,
                                sep="\t",
                                phenotype_col = "Phenotype",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                eaf_col = "eaf",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                units_col = "units",
                                samplesize_col = "samplesize",
                                id_col = "id",
                                min_pval = 1e-200,
                                log_pval = FALSE)
save(smkces_harm,file="./smkces_harm.RData")


# EXPOSURE : agesmk #

agesmk<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/2smr/gwas_agesmk_liu2019.csv",
                  header=TRUE,sep=";",dec=".")
agesmk$Phenotype<-c("agesmk")
agesmk$units<-c("units")
agesmk$id<-c("Liu 2019")
agesmk<-rename.vars(agesmk,
                    from=c("chr","rsid","effect_allele_freq","n"),
                    to=c("chrom","SNP","eaf","samplesize"))
agesmk<-agesmk[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
                  "units","samplesize","id")]
write.table(agesmk,"./agesmk.txt",sep="\t",row.names=FALSE, quote=FALSE)

agesmk_harm<-read_exposure_data(filename = "./agesmk.txt",
                                clump = FALSE,
                                sep="\t",
                                phenotype_col = "Phenotype",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                eaf_col = "eaf",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",
                                units_col = "units",
                                samplesize_col = "samplesize",
                                id_col = "id",
                                min_pval = 1e-200,
                                log_pval = FALSE)
save(agesmk_harm,file="./agesmk_harm.RData")


# OUTCOME : Subfertility #

# Women #

subf_mom<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_CVRF/2smr/meta_analysis/gwama_mom.txt.out",
                                   header=TRUE,sep="\t"))
subf_mom<-rename.vars(subf_mom,
                      from=c("rs_number","reference_allele","OR","OR_se","p.value","n_samples"),
                      to=c("SNP","effect_allele","or","or_se","pval","samplesize"))
subf_mom<-subf_mom[,c("SNP","effect_allele","other_allele","eaf","or","or_se","pval","samplesize")]
subf_mom$beta<-log(subf_mom$or)
subf_mom$se<-subf_mom$or_se/subf_mom$or
subf_mom$Units<-c("Units")
subf_mom$Phenotype<-c("Subfertility")
write.table(subf_mom,"./all_smkinit/meta_analysis/subf_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)

subf_mom_harm <- read_outcome_data(
  snps = smkinit_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(smkinit_harm, subf_mom_harm, action = 2)
save(dat,file="./clean_databases/smkinit_subf_mom_2smr.RData")


subf_mom_harm <- read_outcome_data(
  snps = cigday_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cigday_harm, subf_mom_harm, action = 2)
dat<-dat[-c(which(dat$SNP=="rs10519203")),] # Outlier in scatterplot
save(dat,file="./clean_databases/cigday_subf_mom_2smr.RData")


subf_mom_harm <- read_outcome_data(
  snps = smkces_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(smkces_harm, subf_mom_harm, action = 2)
save(dat,file="./clean_databases/smkces_subf_mom_2smr.RData")


subf_mom_harm <- read_outcome_data(
  snps = agesmk_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(agesmk_harm, subf_mom_harm, action = 2)
save(dat,file="./clean_databases/agesmk_subf_mom_2smr.RData")


# Men #

subf_dad<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_CVRF/2smr/meta_analysis/gwama_dad.txt.out",
                                   header=TRUE,sep="\t"))
subf_dad<-rename.vars(subf_dad,
                      from=c("rs_number","reference_allele","OR","OR_se","p.value","n_samples"),
                      to=c("SNP","effect_allele","or","or_se","pval","samplesize"))
subf_dad<-subf_dad[,c("SNP","effect_allele","other_allele","eaf","or","or_se","pval","samplesize")]
subf_dad$beta<-log(subf_dad$or)
subf_dad$se<-subf_dad$or_se/subf_dad$or
subf_dad$Units<-c("Units")
subf_dad$Phenotype<-c("Subfertility")
write.table(subf_dad,"./all_smkinit/meta_analysis/subf_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)

subf_dad_harm <- read_outcome_data(
  snps = smkinit_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(smkinit_harm, subf_dad_harm, action = 2)
save(dat,file="./clean_databases/smkinit_subf_dad_2smr.RData")


subf_dad_harm <- read_outcome_data(
  snps = cigday_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cigday_harm, subf_dad_harm, action = 2)
save(dat,file="./clean_databases/cigday_subf_dad_2smr.RData")


subf_dad_harm <- read_outcome_data(
  snps = smkces_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(smkces_harm, subf_dad_harm, action = 2)
save(dat,file="./clean_databases/smkces_subf_dad_2smr.RData")


subf_dad_harm <- read_outcome_data(
  snps = agesmk_harm$SNP,
  filename = "./all_smkinit/meta_analysis/subf_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(agesmk_harm, subf_dad_harm, action = 2)
save(dat,file="./clean_databases/agesmk_subf_dad_2smr.RData")


######################
### TWO SAMPLE MRs ###
######################

setwd("N:/data/durable/Projects/Magnus_MR_smoking/2smr")

vars01<-c("smkinit","smkinit","agesmk","agesmk","smkces","smkces","cigday","cigday")
vars02<-c("mom","dad","mom","dad","mom","dad","mom","dad")

tab<-NULL
for(i in 1:length(vars01))
  
{
  namedat<-paste("./clean_databases/",vars01[i],"_subf_",vars02[i],"_all_2smr.RData",sep="")
  load(namedat)
  
  mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_results$pval<-pval_guapa(mr_results$pval)
  mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
  mr_results$or<-risk_se_ic_guapa(mr_results$b,mr_results$se)
  forest<-risk_se_ic_guapa2(mr_results$b[1],mr_results$se[1])
  mr_raps<-mr.raps(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome,diagnosis=FALSE)
  mr_raps_beta<-paste(guapa(mr_raps$beta.hat)," (",guapa(mr_raps$beta.se),")",sep="")
  mr_raps_or<-risk_se_ic_guapa(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_pval<-pval_guapa(mr_raps$beta.p.value)
  q1<-paste(round(mr_heterogeneity(dat)$Q[2],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[2]),")",sep="")
  q2<-paste(round(mr_heterogeneity(dat)$Q[1],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[1]),")",sep="")
  plei<-pval_guapa(mr_pleiotropy_test(dat)$pval)
  
  tab<-rbind(tab,cbind(mr_results$nsnp[1],
                       mr_results$beta[1],mr_results$or[1],mr_results$pval[1],
                       mr_results$beta[2],mr_results$or[2],mr_results$pval[2],
                       mr_results$beta[3],mr_results$or[3],mr_results$pval[3],
                       mr_results$beta[4],mr_results$or[4],mr_results$pval[4],
                       mr_raps_beta,mr_raps_or,mr_raps_pval,plei,q1,q2,forest))
}

colnames(tab)<-c("SNPs","IVW_b","IVW_or","IVW_p","Egger_b","Egger_or","Egger_p",
                 "WMe_b","WMe_or","Wme_p","WMo_b","WMo_or","WMo_p",
                 "RAPS_b","RAPS_or","RAPS_p","Egger_pleio","Cochran_Q","Rucker_Q","forestplot")
rownames(tab)<-paste(vars01,"_",vars02,sep="")
write.table(tab,file="./results/2smr.csv",sep=";",col.names=NA)


### PLOTS ###

# Women, smkinit #

load("./clean_databases/smkinit_subf_mom_all_2smr.RData")
mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/smkinit_subf_mom_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/smkinit_subf_mom_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/smkinit_subf_mom_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Men, smkinit #

load("./clean_databases/smkinit_subf_dad_all_2smr.RData")
mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/smkinit_subf_dad_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/smkinit_subf_dad_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/smkinit_subf_dad_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Women, cigday #

load("./clean_databases/cigday_subf_mom_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/cigday_subf_mom_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/cigday_subf_mom_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/cigday_subf_mom_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Men, cigday #

load("./clean_databases/cigday_subf_dad_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/cigday_subf_dad_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/cigday_subf_dad_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/cigday_subf_dad_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Women, smkces #

load("./clean_databases/smkces_subf_mom_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/smkces_subf_mom_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/smkces_subf_mom_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/smkces_subf_mom_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Men, smkces #

load("./clean_databases/smkces_subf_dad_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/smkces_subf_dad_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/smkces_subf_dad_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/smkces_subf_dad_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Women, agesmk #

load("./clean_databases/agesmk_subf_mom_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/agesmk_subf_mom_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/agesmk_subf_mom_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/agesmk_subf_mom_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()


# Men, agesmk #

load("./clean_databases/agesmk_subf_dad_all_2smr.RData")

mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))

jpeg("./results/plots/agesmk_subf_dad_all_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/plots/agesmk_subf_dad_all_forestplot.jpg", width = 4000, height = 20000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/plots/agesmk_subf_dad_all_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()
