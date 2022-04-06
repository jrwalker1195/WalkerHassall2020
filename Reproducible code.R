install.packages("lintr")
library(lintr)
fp<-"C:/Users/jrwal/MSc Summer Project/Stats and R code/Project R"
tidy_source(fp,comment=TRUE,file="new.R")



#Loading files
wc<-read.csv("wc.csv")
inverts<-read.csv("inverts.csv")
distance<-read.csv("SiteLockDistance.csv")
survival<-read.csv("daphniasurvival.csv")
distance2<-read.csv("SiteLockDistance2.CSV")

#Calculating means of pollutants for each site
wcmean<-aggregate(wc[,2:12],list(wc$Site),mean)

#Testing normality of pollutant residuals
NH4mod<-lm(wcmean$NH4~wcmean$Group.1)
NH4res<-NH4mod$residuals
shapiro.test(NH4res)

NO2mod<-lm(wcmean$NO2~wcmean$Group.1)
NO2res<-NO2mod$residuals
shapiro.test(NO2res)

NO2.NO3mod<-lm(wcmean$NO2.NO3~wcmean$Group.1)
NO2.NO3res<-NO2.NO3mod$residuals
shapiro.test(NO2.NO3res)

PO4mod<-lm(wcmean$PO4~wcmean$Group.1)
PO4res<-PO4mod$residuals
shapiro.test(PO4res)

Camod<-lm(wcmean$Ca~wcmean$Group.1)
Cares<-Camod$residuals
shapiro.test(Cares)

Namod<-lm(wcmean$Na~wcmean$Group.1)
Nares<-Namod$residuals
shapiro.test(Nares)

Znmod<-lm(wcmean$Zn~wcmean$Group.1)
Znres<-Znmod$residuals
shapiro.test(Znres)

Mnmod<-lm(wcmean$Mn~wcmean$Group.1)
Mnres<-Mnmod$residuals
shapiro.test(Mnres)

Femod<-lm(wcmean$Fe~wcmean$Group.1)
Feres<-Femod$residuals
shapiro.test(Feres)

OCmod<-lm(wcmean$OC~wcmean$Group.1)
OCres<-OCmod$residuals
shapiro.test(OCres)

ICmod<-lm(wcmean$IC~wcmean$Group.1)
ICres<-ICmod$residuals
shapiro.test(ICres)

#all normally distributed except Zn

#correlation matrix for colinearity between pollutants
install.packages("Hmisc")
library("Hmisc")
wcmean_cor<-cor(wcmean[,c(2:12)])
wcmean_cor_hmisc<-rcorr(as.matrix(wcmean[,c(2:12)]),type="pearson")
wcmean_cor_r<-wcmean_cor_hmisc$r #Correlation coefficients
wcmean_cor_p<-wcmean_cor_hmisc$P #P-values


#sqrt transforming Zn and retesting normality of residuals
Zn3<-sqrt(wc$Zn)
Zn3mod<-lm(Zn3~wc$Site)
Zn3res<-Zn3mod$residuals
shapiro.test(Zn3res) #still not normal

#rank matrix for colinearity between pollutants because Zn not normally distributed
wcmean_rank_hmisc<-rcorr(as.matrix(wcmean[,c(2:12)]),type="spearman")
wcmean_rank_r<-wcmean_rank_hmisc$r #coefficients
wcmean_rank_p<-wcmean_rank_hmisc$P #P-values


#Calculating diversity indices
install.packages("vegan")
library("vegan")
shannon_d<-diversity(inverts, index = "shannon")
write.csv(shannon_d,file='shannon_d.csv')

#calculated means in XL "invertdiversity"

invertdiversity<-read.csv("invertdiversity.csv")

#Testing normality of diversity residuals
Shannonmod<-lm(invertdiversity$Shannon~invertdiversity$Site)
Shannonres<-Shannonmod$residuals
shapiro.test(Shannonres) #normal

#Testing diversity against site
plot(invertdiversity$Site,invertdiversity$Shannon,xlab="Site",ylab="Diversity")
cor.test(invertdiversity$Site,invertdiversity$Shannon)
abline(lm(Shannon~Site,data=invertdiversity)) 

#Inverts mean
invertmean<-aggregate(inverts[,2:36],list(inverts$Site),mean)

# Calculate principle components that explain the correlated chemistry variables
TestData<-wcmean[,c(2:12)]
# Scale the chemistry so that magnitude of variables doesn't skew PCA
TestData<-scale(TestData)

# Carry out PCA
PCATest<-princomp(TestData)

# Plot PCA
biplot(PCATest,xlab="PC1",ylab="PC2")

# We can then extract the points from the plot
PCAPoints<-PCATest$scores[,c(1:2)]
plot(PCAPoints)
# ...and add the site numbers to see where they fall
plot(c(1:15),PCAPoints[,1],xlab="Site",ylab="PC1")
abline(lm(PCAPoints[,1]~c(1:15)))
plot(c(1:15),PCAPoints[,2],xlab="Site",ylab="PC2")
abline(lm(PCAPoints[,2]~c(1:15)))

# We can then plot the PCs against the location on the canal and test where there is a
# directional change along the stretch in PC1
plot(c(1:15),PCAPoints[,1])
TestLM<-lm(PCAPoints[,1]~c(1:15))
summary(TestLM)
abline(lm(PCAPoints[,1]~c(1:15))) #Significant change in PC1 out from city centre

# and PC2
plot(c(1:15),PCAPoints[,2])
TestLM<-lm(PCAPoints[,2]~c(1:15))
summary(TestLM)
abline(lm(PCAPoints[,2]~c(1:15))) #no significant change in PC2 out from city centre

# and PC3
PCAPoints2<-PCATest$scores[,c(1:3)]
plot(c(1:15),PCAPoints2[,3])
TestLM<-lm(PCAPoints2[,3]~c(1:15))
summary(TestLM)
abline(lm(PCAPoints2[,3]~c(1:15)))

#Testing KMO to validate PCA test
install.packages("psych")
library(psych)
KMO(wcmean) #KMO>0.5 for NH4, NO2 and Na


# Try plotting chemistry PCA against diversity
plot(PCAPoints[,1],invertdiversity$Shannon,xlab="PC1",ylab="Diversity")
abline(lm(invertdiversity$Shannon~PCAPoints[,1]))
plot(PCAPoints[,2],invertdiversity$Shannon,xlab="PC2",ylab="Diversity")
abline(lm(invertdiversity$Shannon~PCAPoints[,2]))
plot(PCAPoints2[,3],invertdiversity$Shannon,xlab="PC3",ylab="Diversity")
abline(lm(invertdiversity$Shannon~PCAPoints2[,3]))

#Testing correlation between PC and diversity
cor.test((PCAPoints[,1]),invertdiversity$Shannon) #significant
cor.test((PCAPoints[,2]),invertdiversity$Shannon) #not significant
cor.test((PCAPoints2[,3]),invertdiversity$Shannon) #not significant

# Hellinger transformation helps to down weight the cells in the community matrix where no
# or small numbers of animals were found
sptrans<-decostand(invertmean,"hellinger")

# CCA/RDA is a form of constrained ordination which allows for hypothesis testing using
# matrices as response variables and (optionally) as predictors. Alternatively, you can
# just use standard model format for the predictors.
PC1<-PCATest$scores[,c(1)]
PC2<-PCATest$scores[,c(2)]
PC3<-PCATest$scores[,c(3)]
TestCCA<-rda(sptrans~PC1+PC2+PC3)
#To get the significance of terms in the RDA/CCA, use anova() but be careful to use the
# "by" argument to specify that you want the values for each term in the model.
anova(TestCCA,by="term")
permutest(TestCCA,by="term")
plot(TestCCA, scaling=3.5)


SiteNames<-paste0("Site",c(1:15))
# Replot the RDA without the site numbers ("sp"=species, "cn"=predictors). Note that I have used "ylim"
# to specify the limits of the y-axis because otherwise the plot is too small to contain all the
# site names.
plot(TestCCA,display=c("sp", "cn"),ylim=c(-1.3,0.6),scaling=3)
text(scores(TestCCA)$sites,rownames(sptrans))





#Testing for relationships between variables in the RDA

#Lumbriculidae looks positively associated with both PC1 and PC2
(plot(PCATest$scores[,c(1)],sptrans$Lumbriculidae,xlab="PC1",ylab="Lumbriculidae"))
spLM<-lm(sptrans$Lumbriculidae~PCATest$scores[,c(1)])
abline(spLM)
summary(spLM)
#only just significant

#Asellidae looks negatively associated with PC2
(plot(PCATest$scores[,c(2)],sptrans$Asellidae))
spLM<-lm(sptrans$Asellidae~PCATest$scores[,c(2)])
abline(spLM)
summary(spLM)
#significant


#Testing normality of SiteLockDistance2 residuals
distancemod<-lm(distance2$Distance2~distance2$Site)
distanceres<-distancemod$residuals
shapiro.test(distanceres) #normal

#parametric test for relationship between diversity and logdistance
diversitydistance<-read.csv("diversitydistance.csv")
cor.test(diversitydistance$Shannon,distance2$Distance2)
plot(distance2$Distance2,diversitydistance$Shannon,xlab="Distance (m)",ylab="Diversity")
#not significant


#Testing effects of locks on community structure
#Group sites within locks together
#Only one sample site between locks 1&2 so not included - nothing to compare to

Lock<-c(rep("A",6),rep("B",2),rep("C",2),rep("D",4))
#where A=Sites 2-7, B=Sites 8-9. C=Sites 10-11, D=Sites 12-15

#Hellinger transformation
newDat<-decostand(invertmean[c(2:15),-1],"hellinger")

newMod<-rda(newDat~Lock)
newMod
anova(newMod) #only one p-value? using "by term" changes value
permutest(newMod)
plot(newMod,scaling=3.5)

#Asellidae looks positively associated with LockA
AsellidaeANOVA<-aov(newDat$Asellidae~Lock)
summary(AsellidaeANOVA)
TukeyHSD(AsellidaeANOVA)
barplot(newDat$Asellidae,names.arg=Lock)
#significant

#Chydoridae looks positively associated with LockD
chydANOVA<-aov(newDat$Chydoridae~Lock)
summary(chydANOVA)
TukeyHSD(chydANOVA)
barplot(newDat$Chydoridae,names.arg=Lock)
#significant

#Gammaridae looks positively associated with LockC
gamANOVA<-aov(newDat$Gammaridae~Lock)
summary(gamANOVA)
TukeyHSD(gamANOVA)
barplot(newDat$Chydoridae,names.arg=Lock)
#significant


#Cox proportional hazards model (survival)
install.packages(c("survival", "survminer"))
library("survival")

CoxData<-read.csv("daphnia_final.csv",header=TRUE)

time<-CoxData$Time
event<-CoxData$Event
site<-CoxData$Site
ind<-CoxData$IndividualID
summary(time)
summary(event)
site1=factor(site,levels=c("C","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","FO","WSA"))
site1
plot(site1,time,xlab="Site",ylab="Survival Time (Days)")

survivalmean<-aggregate(CoxData[,2],list(CoxData$Site),mean)
write.csv(survivalmean,file="survivalmean.csv")
site2=factor(survivalmean$Group.1,levels=c("FO","1","2","3","4","5","6","7","8","9","10","11","WSA","12","13","14","15"))
plot(site2,survivalmean$x,xlab="Site",ylab="Survival Time")
survivallm<-(lm(survivalmean$x~site2))
abline(lm(survivalmean$x~site2))

site2numeric<-as.numeric(site2) # convert site to continuous variable
plot(site2numeric,survivalmean$x)
survivallm<-(lm(survivalmean$x~site2numeric))
summary(survivallm)
abline(lm(survivalmean$x~site2numeric))

survivalmean2<-read.csv("survivalmean2.csv")
site3=factor(survivalmean2$Site,levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))
plot(site3,survivalmean2$Survival)
survivallm<-(lm(survivalmean2$Site~site3))
summary(survivallm)
abline(lm(site3~survivalmean2$Site))

site3numeric<-as.numeric(site3)
plot(site3numeric,survivalmean2$Survival,xlab="Site",ylab="Survival Time (Days)")
survivallm<-(lm(survivalmean2$Survival~site3numeric))
summary(survivallm)
abline(lm(survivalmean2$Survival~site3numeric))

survivalmean3<-read.csv("survivalmean3.csv")
site4=factor(survivalmean3$Site,levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))
site4numeric<-as.numeric(site4)
plot(site4numeric,survivalmean3$Survival)
survivallm<-(lm(survivalmean3$Site~site4numeric))
summary(survivallm)
abline(lm(site4numeric~survivalmean3$Site))

# Kaplan-Meier non-parametric analysis
kmsurvival<-survfit(Surv(time,event)~1)
summary(kmsurvival)
plot(kmsurvival,xlab="Time (Days)",ylab="Survival Probability")

# Kaplan-Meier non-parametric analysis by group
kmsurvival_group<-survfit(Surv(time,event)~site)
summary(kmsurvival_group)
plot(kmsurvival_group,xlab="Time",ylab="Survival Probability")
survdiff(Surv(time,event)~site)

#Cox proportional hazard model
coxph(formula=Surv(time,event)~site1) # gives risk of death compared with control


#Testing survival against wc PC's
survivalmean3<-read.csv("survivalmean3.csv")
survivalmean3mod<-lm(survivalmean3$Survival~survivalmean3$Site)
survivalmean3res<-survivalmean3mod$residuals
shapiro.test(survivalmean3res)

plot(PCAPoints[1:12,1],survivalmean3$Survival,xlab="PC1",ylab="Survival")
abline(lm(survivalmean3$Survival~PCAPoints[1:12,1]))
plot(PCAPoints[1:12,2],survivalmean3$Survival,xlab="PC2",ylab="Survival")
abline(lm(survivalmean3$Survival~PCAPoints[1:12,2]))
plot(PCAPoints2[1:12,3],survivalmean3$Survival,xlab="PC3",ylab="Survival")
abline(lm(survivalmean3$Survival~PCAPoints2[1:12,3]))

cor.test((PCAPoints[1:12,1]),survivalmean3$Survival)
cor.test((PCAPoints[1:12,2]),survivalmean3$Survival)
cor.test((PCAPoints2[1:12,3]),survivalmean3$Survival)

#Testing survival against diversity
diversity2<-read.csv("invertdiversity2.csv")
plot(diversity2$Shannon,survivalmean3$Survival,xlab="Diversity",ylab="Survival")
abline(lm(survivalmean3$Survival~diversity2$Shannon))

cor.test(diversity2$Shannon,survivalmean3$Survival)
