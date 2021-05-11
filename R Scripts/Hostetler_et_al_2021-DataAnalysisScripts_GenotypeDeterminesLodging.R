#Genotype determines root lodging susceptibility
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/")

#Library Load###
library(dplyr)
library(arsenal)
library(lme4)
library(data.table)
library(dplyr)

#Figure 2A-B####
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv(file = "LodgingData_2020_long_02162021.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
lodging = subset(lodging, Lodging.Index<5)
accessions <- read.csv(file = "Unique_55Genotypes_Plots.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
accessions = accessions[!duplicated(accessions$Accession),]
accessions$Genotype = accessions$Accession
accessions$Accession = NULL
lodging$Plot = as.character(lodging$Plot)
summary = summary(comparedf(lodging, accessions, by="Genotype")) #compare lodging data and accessions list and only keep genotypes used in this study 
notShared = summary$obs.table
lodging = lodging[!lodging$Genotype %in% notShared$Genotype,] #remove any genotypes from the lodging file that are not shared in accessions list 
hist(lodging$Lodging.Index)
df1 = subset(lodging, Lodged=="Y")
df2 = subset(lodging, Lodged=="N")
df1$Lodged2 = "1" 
df2$Lodged2 = "0"
lodging = rbind(df1, df2)
lodging %>% count(Lodged2)
str(lodging)
lodging$Lodged2 = as.numeric(lodging$Lodged2)
write.csv(lodging, file = "Figure2A-B_04122021.csv", quote = F, row.names = T)


#Lodging heritability#### 
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "LodgingData_2020.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
head(data)
colnames(data)
data = data[,c(4,10,11)]
str(data)
#prepping BLUPS file
BLUPS<-data.frame(data$Genotype)
BLUPS<-unique(BLUPS[,1])
BLUPS<-data.frame(BLUPS)
colnames(BLUPS)<-"Genotype" #rename column to Genotype so it matches data file
#Creating a data frame for variances and heritability
statis <- matrix(NA,nrow=7,ncol=1)
rownames(statis) <- c("varErr","varG","h2","LSD","pval","var","stddev")
stat<-c("varErr","varG","h2","LSD","pval","var","stddev")
statistics<-data.frame(stat)
#creating loop
name<-colnames(data)
name
name<-name[c(-1)] #Remove all columns that are not traits
for(i in 1:length(name)){
  data5<-data[c(1,i+1)]
  names(data5)[2]<-paste("trait")  
  
  fm <- lmer(trait ~ (1|Genotype), data=data5)
  
  fixEffect <- as.matrix(fixef(fm)) 
  randEffect <- as.matrix(ranef(fm))
  interc <- fixEffect[1,1]
  EntryEffect <- randEffect[["Genotype",1]]
  BLUP <- interc+EntryEffect 
  names(BLUP)[1]<-paste(name[i]) 
  BLUP$Genotype <-rownames(BLUP)
  BLUPS <- merge(BLUPS,BLUP,by="Genotype",all.x = TRUE)
  
  varcorr <- VarCorr(fm) 
  
  varErr <- attr(varcorr,'sc')^2 
  
  varG <- as.vector(varcorr$'Genotype') 
  
  LSD <- 1.96*sqrt(varErr) #calculate least significant difference; Z score for 95% confidence=1.96
  coefs <- data.frame(coef(summary(fm))) 
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))) 
  pval<-coefs[,4]
  var<- var(data5$trait,na.rm=TRUE) 
  stddev<- sqrt(var)
  
  h2 <- varG/(varG + varErr) #calculate heritability w/ broad sense calculation
  statis[,1] <- round(c(varErr,varG,h2,LSD,pval,var,stddev),9) 
  statis2 <- data.frame(statis,stat) 
  colnames(statis2)[1]<-paste(name[i])
  statistics <- merge(statistics,statis2,by="stat")
}
statistics
write.table(statistics, file="TableS4_LodgingData.csv", sep=",", row.names = FALSE)