#cone ratios MASTER clean analyses

#setup####################################################################################################
library(INLA)
library(AnimalINLA)
library(dplyr)
library(MCMCglmm)
library(kinship2)
library(MasterBayes)
library(ggplot2)
library(gridExtra)
library(pedigree)
library(htmlTable)

#import and clean data###################################################################################################################################################################################
phenoRAW <- read.csv("~/Dropbox/My Mac (ecooperimac27.anthro.nyu.edu)/Documents/rhesus_ocular/ddPCR.csv")
#remove infants and clean up column names for pheno data
pheno<-phenoRAW %>%
  dplyr::select(ID=Cayo,year=Year, sex=Sex,age=Age,Lratio,Sratio)%>%
  filter(age>0) %>%
  mutate(animal=ID)
pheno$year<-as.factor((pheno$year))

#clean pedigree#####
pedRAW <- read.csv("C:/Users/u6354548/Dropbox/My Mac (ecooperimac27.anthro.nyu.edu)/Documents/rhesus_ocular/ped_jan2021.csv")
#for calculating x-chromosome matrix, father (heterozygous parent) need be in column 2
ped<-pedRAW %>%
  dplyr::select(ID=AnimalId,sire=Sire,dam=Dam,sex=Sex)
ped<-na_if(ped,y="?")
ped<-na_if(ped,y="")
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


#import models##############################################################################################
Sratio_base<-readRDS("Sratio_base.rds")
Sratio_animal<-readRDS("Sratio_animal.rds")

Lratio_base<-readRDS("Lratio_base.rds")
Lratio_animal<-readRDS("Lratio_animal.rds")

#L-ratio base model#########################################################################################
prior1b<-list(R=list(V=1, nu=0.002))

Lratio_base<- MCMCglmm(Lratio ~ sex+age+as.factor(year),data = pheno,pr=TRUE,
                   nitt=2500000,thin=1000,burnin= 10000,prior=prior1b,verbose = TRUE)
summary(Lratio_base)
plot(Lratio_base$VCV)
Lratio_base$DIC

mean(subset(pheno$Lratio,pheno$sex=="F"))

saveRDS(Lratio_base,"Lratio_base.rds")

#L-ratio animal model########################################################################################
#this is Lratio_Aonly in coneRatio_testMods_Lratio

Lratio_animal<- MCMCglmm(scale(Lratio) ~ sex+age+year, 
                        random = ~ animal,
                        pedigree = ped[,1:3],data = pheno,
                        nitt=2500000,thin=1000,burnin= 100000,prior = prior1a,verbose = TRUE)
saveRDS(Lratio_animal,"Lratio_animal.rds")

plot(Lratio_animal$VCV)
Lratio_animal$DIC
summary(Lratio_animal)


herit_L<-(Lratio_animal$VCV[,'animal']/
            ((Lratio_animal$VCV[,'animal'])+(Lratio_animal$VCV[,'units'])))
HPDinterval(herit_L)
mean(herit_L)

#S-ratio base model########################################################################################
#this is 'Sratio2' in coneRatio_Sratio.R
prior1b<-list(R=list(V=1, nu=0.002))

Sratio_base<- MCMCglmm(Sratio ~ sex+age+year,data = pheno,
                   nitt=1500000,thin=1000,burnin= 10000,prior = prior1b,verbose = TRUE)
summary(Sratio_base)
plot(Sratio_base$VCV)
Sratio_base$DIC


#S-ratio animal model###############################################################################################
#this is 'Sratio1' in coneRatio_Sratio.R
prior1a<-list(G=list(G1=list(V = 1,nu=0.002)),
              R=list(V=1, nu=0.002))

Sratio_animal<- MCMCglmm(Sratio ~ sex+age+as.factor(year), 
                   random = ~ animal,
                   pedigree = ped[,1:3],data = pheno,
                   nitt=1500000,thin=1000,burnin= 10000,prior = prior1a,verbose = TRUE)
summary(Sratio_animal)
plot(Sratio_animal$VCV)

Sratio_animal$DIC

herit_S<-(Sratio_animal$VCV[,'animal']/
            ((Sratio_animal$VCV[,'animal'])+(Sratio_animal$VCV[,'units'])))
HPDinterval(herit_S)
mean(herit_S)


#table: fixed effects######################################################################################
#fixed effects from the two best models
#Lratio_base & Sratio_animal

fixedEffects_S<- list(Sratio_animal$Sol[,'(Intercept)'],Sratio_animal$Sol[,'sexM'],Sratio_animal$Sol[,'age'],Sratio_animal$Sol[,'year2019'])

fixedEffects_L<- list(Lratio_base$Sol[,'(Intercept)'],Lratio_base$Sol[,'sexM'],Lratio_base$Sol[,'age'],Lratio_base$Sol[,'as.factor(year)2019'])

posterior.mode(fixedEffects_S)
mean(Sratio_animal$Sol[,"year2019"])

tab_FE<-data.frame("Fixed Effects"=c("Intercept","Sex (Male)","Age","Year (2019)"),
                    "S-Ratio"=unlist(lapply(fixedEffects_S,function(x){
                      paste0(round(posterior.mode(x),4)," [",round(HPDinterval(x)[1], 3)," - ",
                             round(HPDinterval(x)[2], 3),"]")})),
                      "L-Ratio"=unlist(lapply(fixedEffects_L,function(x){
                        paste0(round(posterior.mode(x),4)," [",round(HPDinterval(x)[1], 3)," - ",
                               round(HPDinterval(x)[2], 3),"]")
                    })))
tab_FE
htmlTable(tab_FE)
#bar graph: variance components###########################################################################
ranef_S<- list(Sratio_animal$VCV[,'animal'],
                Sratio_animal$VCV[,'units'],
                Sratio_animal$VCV[,'animal'],
                Sratio_animal$VCV[,'units'])

#summary(ranef_S)

#put in a table
tabranef <- data.frame(term=c("Animal","Residual","Animal","Residual"),
                       estimate=unlist(lapply(allranef, mode)),
                       CIlow=unlist(lapply(allranef, function(x){
                         round(HPDinterval(x)[1], 3)
                       })),
                       CIhigh=unlist(lapply(allranef, function(x){
                         round(HPDinterval(x)[2], 3)
                       })))
totPhen<-tabranef%>%
  group_by(age)%>%
  summarise(totVar=sum(estimate))
totPhen


#unstacked bar graph with error bars
term.labs <- c("Residual","Individual Environment","Additive Genetic","Year")
names(term.labs) <- c("res", "aID","aani","ayr")


ggplot(data=ranef_S,aes(x=as.factor(age),y=estimate))+
  facet_grid(.~term,labeller=labeller (term=term.labs))+
  geom_col(aes(fill=term))+
  ylab("Variance")+
  xlab("Age Group")+
  geom_errorbar(aes(ymin=CIlow,ymax=CIhigh,x=as.factor(age)))+
  theme_minimal()+
  scale_fill_brewer(palette = 18)+
  theme(legend.position = "none")

#forest plot: fixed effects###########################################################################
fixedEffects_S<- list(Sratio_animal$Sol[,'(Intercept)'],Sratio_animal$Sol[,'sexM'],Sratio_animal$Sol[,'age'],Sratio_animal$Sol[,'year2019'])

fixedEffects_L<- list(Lratio_base$Sol[,'(Intercept)'],Lratio_base$Sol[,'sexM'],Lratio_base$Sol[,'age'],Lratio_base$Sol[,'as.factor(year)2019'])

fixedEffects_all<- list(Sratio_animal$Sol[,'(Intercept)'],Sratio_animal$Sol[,'sexM'],Sratio_animal$Sol[,'age'],Sratio_animal$Sol[,'year2019']
                        ,Lratio_base$Sol[,'(Intercept)'],Lratio_base$Sol[,'sexM'],Lratio_base$Sol[,'age'],Lratio_base$Sol[,'as.factor(year)2019'])

FE_sum<-data.frame(effect=c("Intercept", "Sex (male)", "Age", "Year (2019)","Intercept", "Sex (male)", "Age", "Year (2019)"),
                   model=c(rep("S",4),rep("L",4)),
                   value=unlist(lapply(fixedEffects_all, function(x){
                     round(mean(x), 5)})),
                   min=unlist(lapply(fixedEffects_all, function(x){
                     round(HPDinterval(x)[1], 5)
                   })),
                   max=unlist(lapply(fixedEffects_all, function(x){
                     round(HPDinterval(x)[2], 5)
                   })))
FE_sum<-filter(FE_sum,effect!="Intercept")
FE_sumS<-filter(FE_sum,model=="S")
FE_sumL<-filter(FE_sum,model=="L")


ggplot(data=FE_sum,aes(x=value,y=effect,group=model))+
  geom_errorbarh(aes(xmin=min,xmax=max))+
  theme_minimal()

library(cowplot)

L<-ggplot(data=FE_sumL,aes(x=value,y=effect))+
  geom_point()+
  geom_errorbarh(aes(xmin=min,xmax=max))+
  geom_vline(xintercept=0,color="red")+
  coord_cartesian(xlim=c(-0.2,0.2))+
  labs(title = 'L:M Ratio')+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())


S<-ggplot(data=FE_sumS,aes(x=value,y=effect))+
  geom_point()+
  geom_errorbarh(aes(xmin=min,xmax=max))+
  geom_vline(xintercept=0,color="red")+
  coord_cartesian(xlim=c(-0.05,0.05))+
  labs(title = 'S:L+M Ratio')+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())


plot_grid(L,S,labels=c('A','B'))




