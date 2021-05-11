###Brace root phenotypes vary between genotypes
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/")

#Library Load
library(arsenal)
library(dplyr)
library(dplyr)
library(agricolae)
library(rcompanion)
library(corrplot)
library(lme4)
library(factoextra)
library(ggfortify)

###Phenotype Profile Prep, Merge with height, Convert to cm ####
cat("\014")
rm(list=ls()) 
ls() 
data1 <- read.csv(file = "2019_root_photo_data_filtered.csv", header = TRUE, na.strings = "NA")
data2 <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio.csv", header = TRUE, na.strings = "NA")
data2_2018 = subset(data2, Year=="2018")
data2_2019 = subset(data2, Year=="2019")
summary = summary(comparedf(data2_2018, data2_2019, by="Accession")) #Compare 2018 and 2019, determine which genotypes are present in both years and have a none/all ratio > 0
notShared = summary$obs.table
data3 = data2[!data2$Accession %in% notShared$Accession,]
data3 = data3[data3$Ratio..None.All....IMU > 0,]
data3$PLOT_ID = data3$Plot.ID
data3$Plot.ID = NULL
data1$ID = paste(data1$PLOT_ID, data1$Plant_Number, sep="_")
data3$ID = paste(data3$PLOT_ID, data3$Plant.Number, sep="_")
summary2 = summary(comparedf(data3, data1, by="ID")) #compare phenotype data and ratio data to determine which genotypes we have biomechanics and phenotype data for 
notShared2 = summary2$obs.table
data4 = data1[!data1$ID %in% notShared2$ID,]
#Subset A and B measurements so we can make a wide format and have one entry per plant
data4_A = subset(data4, AorB == "A")
data4_B = subset(data4, AorB == "B")
data5 = merge(data4_A, data4_B, by = "ID") #X = A, Y = B
colnames(data5)
data6 = data5[,c(1,3)] #Create new data frame with A and B measurements merged appropriately 
data6$Plant_Number = data6$Plant_Number.x
data6$Plant_Number.x = NULL
#Combine matching phenotypes 
data6$stalk_width = ((data5$stalk_width.x + data5$stalk_width.y) / 2)
data6$leftmost_single_root_width_A = (data5$leftmost_single_root_width.x)
data6$leftmost_single_root_width_B = (data5$leftmost_single_root_width.y)
data6$num_whorls = ((data5$num_whorls.x + data5$num_whorls.y) / 2)
data6$w1 = ((data5$w1.x + data5$w1.y))
data6$w2 = ((data5$w2.x + data5$w2.y))
data6$w3 = ((data5$w3.x + data5$w3.y))
data6$w4 = ((data5$w4.x + data5$w4.y))
data6$root_heightonstalk_left = ((data5$root_heightonstalk_left.x + data5$root_heightonstalk_right.y) / 2)
data6$root_heightonstalk_right = ((data5$root_heightonstalk_right.x + data5$root_heightonstalk_left.y) / 2)
data6$stalk_to_leftmost_rootgrounding = ((data5$stalk_to_leftmost_rootgrounding.x + data5$stalk_to_rightmost_rootgrounding.y) / 2)
data6$stalk_to_rightmost_rootgrounding = ((data5$stalk_to_rightmost_rootgrounding.x + data5$stalk_to_leftmost_rootgrounding.y) / 2)
data6$root_angledeg_left = ((data5$root_angledeg_left.x + data5$root_angledeg_right.y) / 2)
data6$root_angledeg_right = ((data5$root_angledeg_right.x + data5$root_angledeg_left.y) / 2)
data6$spread_width = ((data5$spread_width.x + data5$spread_width.y) / 2)
#Create a data frame that has both phenotype data and ratio data 
colnames(data3)
data7 = select(data3, Year, Accession, Plant.Age, Load.Cell.Height, Operator, Ratio..None.All....IMU, 
               Brace.Root.Whorls.in.the.Soil, PLOT_ID, ID)
data8 = merge(data7, data6, by = "ID")
#Add EI data for status A (all brace roots intact) to data8 file 
data9 <- read.csv(file = "Inbred_Subpop_2Years_EI_All.csv", header = TRUE, na.strings = "NA")
data9_2019 = subset(data9, Year=="2019") #only work with 2019 data
data9_A = subset(data9_2019, Brace.Root.Status=="A") #only work with raw A data (all brace roots intact)
data9_A$ID = paste(data9_A$Plot.ID, data9_A$Plant.Number, sep="_") 
summary3 = summary(comparedf(data9_A, data8, by="ID"))  
notShared3 = summary3$obs.table #remove individuals that are not present in both data frames
data10 = data9_A[!data9_A$ID %in% notShared3$ID,]
data10 = data10[data10$Flexural.Rigidity..EI....IMU > 0,] #remove any measurements that have an EI < 0
data11 = data10[,c(14,28)] #only select EI value and ID
data12 = merge(data8, data11, by = "ID") #merge two data frames 
colnames(data12)
data13 = data12[,c(1:3,9:10,26,7,8,11,14:18,25)]
#Create phenotypic profiles
data13$single_root_width = (data12$leftmost_single_root_width_A + data12$leftmost_single_root_width_B)/2
data13$root_heightonstalk = (data12$root_heightonstalk_left + data12$root_heightonstalk_right)/2
data13$stalk_to_rootgrounding = (data12$stalk_to_leftmost_rootgrounding + data12$stalk_to_rightmost_rootgrounding)/2
data13$root_angle = (data12$root_angledeg_left + data12$root_angledeg_right)/2
#Merge with height data
data14 <- read.csv(file = "2019PC_Height.csv", header = TRUE, na.strings = "NA")
head(data14)
data14$ID = paste(data14$plot, data14$plant, sep="_")
data14 = data14[,c(8,4)]
data13 = merge(data13, data14, by = "ID")
data13[data13$ID=="1318_2", "height.cm."] = 146 #fix data entry error; was 146 not 46
#Convert all inches to cm
data13$stalk_width.cm = (data13$stalk_width)*2.54
data13$spread_width.cm = (data13$spread_width)*2.54
data13$single_root_width.cm = (data13$single_root_width)*2.54
data13$root_heightonstalk.cm = (data13$root_heightonstalk)*2.54
data13$stalk_to_rootgrounding.cm = (data13$stalk_to_rootgrounding)*2.54
data13$stalk_to_rootgrounding.cm = (data13$stalk_to_rootgrounding)*2.54
colnames(data13)
data15 = data13[,c(2,1,4:5,3,6:8,10:13,19:25)]
colnames(data15)
write.csv(data15, file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", row.names = FALSE, quote = F)

##Determine if phenotype significantly differs among genotypes#### 
data16 = data15
colnames(data16)
data17 = data16[,c(3,5:19)]
#Note repeat these steps for each trait
colnames(data17)
data17$x = data17$w3
shapiro.test(data17$x)
hist(data17$x)
data17$tukey <- transformTukey(
  data17$x,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data17$tukey)
res.aov = aov(data17$tukey ~ Accession, data = data17)
summary(res.aov)
tukey.test = HSD.test(res.aov, trt = 'Accession')
tukey.test
output1 = as.data.frame(tukey.test$means)
output2 = as.data.frame(tukey.test$groups)
#Summary for not-transformed data
res.aov = aov(data17$x ~ Accession, data = data17)
summary(res.aov)
tukey.test = HSD.test(res.aov, trt = 'Accession')
output1 = as.data.frame(tukey.test$means)
output1

##PCA ####
df = data15
df = df[,c(8:11,13:19)]
pca <- prcomp(df, scale. = TRUE)
var <- get_pca_var(pca)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
pdf("Figure3.pdf", width = 8.5, height = 11)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
dev.off()
#Table S8 and S9 Information
summary(pca) #S8 Info
var <- get_pca_var(pca)
var$contrib #S9 Info

#Correlation of 2019 phenotypes (all reps inluded)#####
cor(df)
corrplot(cor(df), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")

#Heritability Info#### 
colnames(data15)
data15 = data15[,c(5,8:19)]
str(data15)
BLUPS<-data.frame(data15$Accession)
BLUPS<-unique(BLUPS[,1])
BLUPS<-data.frame(BLUPS)
colnames(BLUPS)<-"Accession" 

statis <- matrix(NA,nrow=7,ncol=1)
rownames(statis) <- c("varErr","varG","h2","LSD","pval","var","stddev")
stat<-c("varErr","varG","h2","LSD","pval","var","stddev")
statistics<-data.frame(stat)

#creating loop
name<-colnames(data15)
name
name<-name[c(-1)]
for(i in 1:length(name)){
  
  data18<-data15[c(1,i+1)] 
  names(data18)[2]<-paste("trait")
  
  fm <- lmer(trait ~ (1|Accession), data=data18)
  
  fixEffect <- as.matrix(fixef(fm)) 
  randEffect <- as.matrix(ranef(fm)) 
  interc <- fixEffect[1,1] 
  EntryEffect <- randEffect[["Accession",1]]
  BLUP <- interc+EntryEffect 
  names(BLUP)[1]<-paste(name[i]) 
  BLUP$Accession <-rownames(BLUP)
  BLUPS <- merge(BLUPS,BLUP,by="Accession",all.x = TRUE)
  
  varcorr <- VarCorr(fm) 
  
  varErr <- attr(varcorr,'sc')^2 
  
  varG <- as.vector(varcorr$'Accession') 
  
  LSD <- 1.96*sqrt(varErr) #calculate least significant difference; Z score for 95% confidence=1.96
  coefs <- data.frame(coef(summary(fm))) 
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))) 
  pval<-coefs[,4]
  var<- var(data18$trait,na.rm=TRUE)
  stddev<- sqrt(var) 
  
  h2 <- varG/(varG + varErr) #calculate heritability w/ broad sense calculation
  statis[,1] <- round(c(varErr,varG,h2,LSD,pval,var,stddev),9)
  statis2 <- data.frame(statis,stat) 
  colnames(statis2)[1]<-paste(name[i])
  statistics <- merge(statistics,statis2,by="stat")
}
write.table(statistics, file="TableS4_PhenotypeData.csv", sep=",", row.names = FALSE)


