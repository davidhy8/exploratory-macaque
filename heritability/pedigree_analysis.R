library(INLA)
#library(AnimalINLA)
library(dplyr)
library(MCMCglmm)
library(kinship2)
library(MasterBayes)
library(ggplot2)
library(gridExtra)
library(pedigree)
library(htmlTable)
library(naniar)
library(stringr)

# Phenotypes
phenoRAW <- read.csv("ped_jan2021_withIOP.csv")
pheno<-phenoRAW %>%
  dplyr::select(ID=animalid,year=dob, sex, behaviormom,sire,dam, iop_re_first, iop_le_first)%>%
  mutate(animal=ID)
#pheno$year<-as.factor((pheno$year))
pheno$year = str_sub(pheno$year, -4,-1)
pheno$year = as.numeric(pheno$year)

phenotypes = subset(phenoRAW, select= c(animalid, iop_re_first, iop_le_first))


# Pedigree
pedRAW <- read.csv("ped_jan2021.csv")
pedRAW = merge(pedRAW, phenotypes, by.x = "AnimalId", by.y = "animalid", all.x=TRUE)

#for calculating x-chromosome matrix, father (heterozygous parent) need be in column 2
ped<-pedRAW %>%
  dplyr::select(ID=AnimalId,sire=Sire,dam=Dam,sex=Sex, behaviormom = BehaviorMom, iop_re_first, iop_le_first)

# ped<-na_if(ped,y="?")
# ped<-na_if(ped,y="")
na_strings <- c("?", "")
ped = ped %>% replace_with_na_all(~.x %in% na_strings)

ped$ID<-gsub("#", "", ped$ID)
ped$sire<-gsub("#", "", ped$sire)
ped$dam<-gsub("#", "", ped$dam)

#females (homozygous) are sex 1, males are sex 2
ped<-ped %>%
  filter(sex!="u")%>%
  filter(sex!="U")%>%
  mutate(sex=ifelse(sex=="F"|sex=="f",1,2))

#clean pedigree and put in correct order
ped<-as.data.frame(ped)
ped1<-prunePed(ped,keep=as.factor(pheno$ID))
ped<-MasterBayes::orderPed(ped1)

# prior1 <- list(R = list(V = 1, nu = 0.002),
#               G = list(G1 = list(V = 1, nu = 0.002)))

prior1.1 <- list(R = list(V = 1, nu = 1),
                       G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior1.4.SEX <- list(
  G = list(G1 = list(V = diag(2), nu = 1.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)),
  R = list(V = diag(2), nu = 1.002)
)
# #pedigree analysis re
# model.simple.re <- MCMCglmm(iop_re_first ~ 1,
#                   random = ~ animal,
#                   family   = "gaussian",
#                   prior    = prior1,
#                   pedigree = ped[,1:3],
#                   data     = pheno,
#                   nitt     = 1500000,
#                   burnin   = 10000,
#                   thin     = 1000,
#                   verbose=TRUE)
# 
# saveRDS(model.simple.re, "model_simple_re.rds")
# plot(model.simple.re[["Sol"]])
# plot(model.simple.re[["VCV"]])
# 
# effectiveSize(model.simple.re[["Sol"]])
# effectiveSize(model.simple.re[["VCV"]])
# 
# herit.simple.re <-
#   model.simple.re[["VCV"]][ , "animal"] /
#   (model.simple.re[["VCV"]][ , "animal"] + model.simple.re[["VCV"]][ , "units"])
# 
# effectiveSize(herit.simple.re)
# mean(herit.simple.re)
# HPDinterval(herit.simple.re)
# 
# # pedigree analysis simple le
# model.simple.le <- MCMCglmm(iop_le_first ~ 1,
#                             random = ~ animal,
#                             family   = "gaussian",
#                             prior    = prior1,
#                             pedigree = ped[,1:3],
#                             data     = pheno,
#                             nitt     = 1500000,
#                             burnin   = 10000,
#                             thin     = 1000,
#                             verbose=TRUE)
# 
# saveRDS(model.simple.le, "model_simple_le.rds")
# plot(model.simple.le[["Sol"]])
# plot(model.simple.le[["VCV"]])
# 
# effectiveSize(model.simple.le[["Sol"]])
# effectiveSize(model.simple.le[["VCV"]])
# 
# herit.simple.le <-
#   model.simple.le[["VCV"]][ , "animal"] /
#   (model.simple.le[["VCV"]][ , "animal"] + model.simple.le[["VCV"]][ , "units"])
# 
# effectiveSize(herit.simple.le)
# mean(herit.simple.le)
# HPDinterval(herit.simple.le)

# Adding fixed effects age & sex
model.fixed.re <- MCMCglmm(iop_re_first ~ sex + year,
                            random = ~ animal,
                            family   = "gaussian",
                            prior    = prior1.1,
                            pedigree = ped[,1:3],
                            data     = pheno,
                            nitt     = 12500000,
                            burnin   = 10000,
                            thin     = 1000,
                            verbose=TRUE)

saveRDS(model.fixed.re, "model_fixed_re.rds")

#graphics.off()
#par("mar")
#par(mar=c(1,1,1,1))

pdf("model_fixed_re_Sol.pdf")
print(plot(model.fixed.re[["Sol"]]))
dev.off()

pdf("model_fixed_re_VCV.pdf")
plot(model.fixed.re[["VCV"]])
dev.off()

effectiveSize(model.fixed.re[["Sol"]])
effectiveSize(model.fixed.re[["VCV"]])

herit.fixed.re <-
  model.fixed.re[["VCV"]][ , "animal"] /
  (model.fixed.re[["VCV"]][ , "animal"] + model.fixed.re[["VCV"]][ , "units"])

autocorr.diag(model.fixed.re[["Sol"]])

effectiveSize(herit.fixed.re)
mean(herit.fixed.re)
HPDinterval(herit.fixed.re)

model.fixed.le <- MCMCglmm(iop_le_first ~ sex + year,
                           random = ~ animal,
                        
                           family   = "gaussian",
                           prior    = prior1.1,
                           pedigree = ped[,1:3],
                           data     = pheno,
                           nitt     = 12500000,
                           burnin   = 10000,
                           thin     = 1000,
                           verbose=TRUE)

saveRDS(model.fixed.le, "model_fixed_le.rds")

pdf("model_fixed_le_Sol.pdf")
plot(model.fixed.le[["Sol"]])
dev.off()

pdf("model_fixed_le_VCV.pdf")
plot(model.fixed.le[["VCV"]])
dev.off()

effectiveSize(model.fixed.le[["Sol"]])
effectiveSize(model.fixed.le[["VCV"]])

herit.fixed.le <-
  model.fixed.le[["VCV"]][ , "animal"] /
  (model.fixed.le[["VCV"]][ , "animal"] + model.fixed.le[["VCV"]][ , "units"])

effectiveSize(herit.fixed.le)
mean(herit.fixed.le)
HPDinterval(herit.fixed.le)

