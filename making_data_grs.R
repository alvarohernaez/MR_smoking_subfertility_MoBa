rm(list=ls())

# Load required libraries #

setwd("N:/data/durable/Projects/Magnus_MR_smoking")


#################################
### GENERATION OF SMKINIT GRS ###
#################################

# MoBa SELECTED SNPs ON SmkInit #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr4.raw",header=TRUE,sep="\t"))
chr05<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr5.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr9.raw",header=TRUE,sep="\t"))
chr10<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr10.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr11.raw",header=TRUE,sep="\t"))
chr12<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr12.raw",header=TRUE,sep="\t"))
chr13<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr13.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr16.raw",header=TRUE,sep="\t"))
chr17<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr17.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr21.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/chr22.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr05<-chr05[,c(1,7:dim(chr05)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr10<-chr10[,c(1,7:dim(chr10)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr12<-chr12[,c(1,7:dim(chr12)[2])]
chr13<-chr13[,c(1,7:dim(chr13)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr17<-chr17[,c(1,7:dim(chr17)[2])]
chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]
chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr05,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr10,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr12,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr13,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr17,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO SmkInit, META-ANALYSIS LIU ET AL, NAT GENET, 2019 #

# Load summary results of the meta-analysis #

smkinit_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkinit/gwas_smkinit_liu2019.csv",header=TRUE,sep=",",dec=".")

names(smkinit_ma)<-tolower(names(smkinit_ma))
smkinit_ma<-rename.vars(smkinit_ma,
                        from=c("tested.allele","tested.allele.freq","pvalue"),
                        to=c("tested_allele_ma","maf_ma","pval"))
smkinit_ma$tested_allele_ma<-tolower(smkinit_ma$tested_allele_ma)


### CHECK THAT THE 355 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(smkinit_ma$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND smkinit_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
smkinit_ma<-merge2(smkinit_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkinit_ma<-na.omit(smkinit_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

smkinit_ma$same_coding<-with(smkinit_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                               ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

smkinit_ma$flip<-with(smkinit_ma,ifelse(beta>=0 & same_coding==1,0,
                                        ifelse(beta<0 & same_coding==0,0,
                                               ifelse(beta>=0 & same_coding==0,1,
                                                      ifelse(beta<0 & same_coding==1,1,NA)))))

smkinit_ma$maf<-with(smkinit_ma,ifelse(same_coding==1,maf_ma,
                                       ifelse(same_coding==0,1-maf_ma,NA)))
smkinit_ma$coef<-with(smkinit_ma,ifelse(beta>0,beta,
                                        ifelse(beta<0,-beta,NA)))


### CALCULATION OF smkinit GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-smkinit_ma[smkinit_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-smkinit_ma[smkinit_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$smkinit_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","smkinit_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_smoking/source_files/smkinit_grs.RData")


#################################
### GENERATION OF AGESMK GRS ###
#################################

# MoBa SELECTED SNPs ON AgeSmk #

chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/chr4.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/chr8.raw",header=TRUE,sep="\t"))

chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]

dat<-merge2(chr02,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO agesmk, META-ANALYSIS LIU ET AL, NAT GENET, 2019 #

# Load summary results of the meta-analysis #

agesmk_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_agesmk/gwas_agesmk_liu2019.csv",header=TRUE,sep=",",dec=".")
names(agesmk_ma)<-tolower(names(agesmk_ma))
agesmk_ma<-rename.vars(agesmk_ma,
                       from=c("tested.allele","tested.allele.freq","pvalue"),
                       to=c("tested_allele_ma","maf_ma","pval"))
agesmk_ma$tested_allele_ma<-tolower(agesmk_ma$tested_allele_ma)


### CHECK THAT THE 10 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(agesmk_ma$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND agesmk_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
agesmk_ma<-merge2(agesmk_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
agesmk_ma<-na.omit(agesmk_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

agesmk_ma$same_coding<-with(agesmk_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                             ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

agesmk_ma$flip<-with(agesmk_ma,ifelse(beta>=0 & same_coding==1,0,
                                      ifelse(beta<0 & same_coding==0,0,
                                             ifelse(beta>=0 & same_coding==0,1,
                                                    ifelse(beta<0 & same_coding==1,1,NA)))))

agesmk_ma$maf<-with(agesmk_ma,ifelse(same_coding==1,maf_ma,
                                     ifelse(same_coding==0,1-maf_ma,NA)))
agesmk_ma$coef<-with(agesmk_ma,ifelse(beta>0,beta,
                                      ifelse(beta<0,-beta,NA)))


### CALCULATION OF agesmk GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-agesmk_ma[agesmk_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-agesmk_ma[agesmk_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$agesmk_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","agesmk_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_smoking/source_files/agesmk_grs.RData")


#################################
### GENERATION OF CIGDAY GRS ###
#################################

# MoBa SELECTED SNPs ON CigDay #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr4.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr9.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr11.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr16.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/chr21.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO cigday, META-ANALYSIS LIU ET AL, NAT GENET, 2019 #

# Load summary results of the meta-analysis #

cigday_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_cigday/gwas_cigday_liu2019.csv",header=TRUE,sep=",",dec=".")
names(cigday_ma)<-tolower(names(cigday_ma))
cigday_ma<-rename.vars(cigday_ma,
                       from=c("tested.allele","tested.allele.freq","pvalue"),
                       to=c("tested_allele_ma","maf_ma","pval"))
cigday_ma$tested_allele_ma<-tolower(cigday_ma$tested_allele_ma)


### CHECK THAT THE 50 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(cigday_ma$rsid,decreasing=FALSE))

dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND cigday_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
cigday_ma<-merge2(cigday_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
cigday_ma<-na.omit(cigday_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

cigday_ma$same_coding<-with(cigday_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                             ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

cigday_ma$flip<-with(cigday_ma,ifelse(beta>=0 & same_coding==1,0,
                                      ifelse(beta<0 & same_coding==0,0,
                                             ifelse(beta>=0 & same_coding==0,1,
                                                    ifelse(beta<0 & same_coding==1,1,NA)))))

cigday_ma$maf<-with(cigday_ma,ifelse(same_coding==1,maf_ma,
                                     ifelse(same_coding==0,1-maf_ma,NA)))
cigday_ma$coef<-with(cigday_ma,ifelse(beta>0,beta,
                                      ifelse(beta<0,-beta,NA)))


### CALCULATION OF cigday GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-cigday_ma[cigday_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-cigday_ma[cigday_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$cigday_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","cigday_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_smoking/source_files/cigday_grs.RData")


#################################
### GENERATION OF SMKCES GRS ###
#################################

# MoBa SELECTED SNPs ON SmkCes #

chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr3.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr9.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr11.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr15.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr20.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/chr22.raw",header=TRUE,sep="\t"))

chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr02,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO smkces, META-ANALYSIS LIU ET AL, NAT GENET, 2019 #

# Load summary results of the meta-analysis #

smkces_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_smoking/gwas/gwas_smkces/gwas_smkces_liu2019.csv",header=TRUE,sep=",",dec=".")
names(smkces_ma)<-tolower(names(smkces_ma))
smkces_ma<-rename.vars(smkces_ma,
                       from=c("tested.allele","tested.allele.freq","pvalue"),
                       to=c("tested_allele_ma","maf_ma","pval"))
smkces_ma$tested_allele_ma<-tolower(smkces_ma$tested_allele_ma)


### CHECK THAT THE 23 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(smkces_ma$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND smkces_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
smkces_ma<-merge2(smkces_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkces_ma<-na.omit(smkces_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

smkces_ma$same_coding<-with(smkces_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                             ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

smkces_ma$flip<-with(smkces_ma,ifelse(beta>=0 & same_coding==1,1,
                                      ifelse(beta<0 & same_coding==0,1,
                                             ifelse(beta>=0 & same_coding==0,0,
                                                    ifelse(beta<0 & same_coding==1,0,NA)))))

smkces_ma$maf<-with(smkces_ma,ifelse(same_coding==1,maf_ma,
                                     ifelse(same_coding==0,1-maf_ma,NA)))
smkces_ma$coef<-with(smkces_ma,ifelse(beta>0,beta,
                                      ifelse(beta<0,-beta,NA)))


### CALCULATION OF smkces GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-smkces_ma[smkces_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-smkces_ma[smkces_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$smkces_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","smkces_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_smoking/source_files/smkces_grs.RData")