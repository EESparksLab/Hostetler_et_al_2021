#Updated script for Hostetler_et_al_2021 PCE-Version 2-2021/12
#Scripts are broken down into section headers associated with the manuscript  
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/")
#Library Load
library(arsenal)
library(dplyr)
library(agricolae)
library(rcompanion)
library(corrplot)
library(lme4)
library(factoextra)
library(ggfortify)
library(car)
library(tidyverse)
library(tidymodels)
library(ranger)
library(kknn)
library(janitor)
library(data.table)
#RESULTS: The brace root contribution to anchorage varies among genotypes####
##ANOVA on the brace root contribution to anchorage - Table S3-S4####
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
data$Plot.ID = as.factor(data$Plot.ID)
data$Year = as.factor(data$Year)
attach(data)
par(mfrow=c(1,1))
lm_ratio <- lm(Ratio..None.All....IMU ~ Accession * Year)
lm_ratio_aov=anova(lm_ratio) 
lm_ratio_aov
lm_ratio_aov=aov(lm_ratio) 
Ratio_Tukey1=TukeyHSD(lm_ratio_aov, which = "Accession:Year") 
TukeyOutput1=Ratio_Tukey1$`Accession:Year`
TukeyOutput1 = as.data.frame(TukeyOutput1)
#write.csv(TukeyOutput1, file = "TukeyOutput1_12152021.csv", row.names = TRUE, quote = F)
HSD.test(lm_ratio_aov, trt = c("Accession", "Year"), console = TRUE)
par(mfrow=c(2,2))
plot(lm_ratio)
par(mfrow=c(1,1))
plot(Ratio..None.All....IMU)
hist(Ratio..None.All....IMU)
hist(lm_ratio$residuals)
detach(data)

##Heritability Info - Table S5 #####
#Ratio
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
head(data)
colnames(data)
data = data[,c(1,2,7,8)]
str(data)
i <- c(1:2)
data[,i] <- apply(data[ ,i], 2,
                  function(x) as.character(x))
BLUPS<-data.frame(data$Accession)
BLUPS<-unique(BLUPS[,1]) 
BLUPS<-data.frame(BLUPS)
colnames(BLUPS)<-"Accession"
statis <- matrix(NA,nrow=9,ncol=1)
rownames(statis) <- c("n Reps","Error Var","Genotypic Var", "Year Var", "Heritability","LSD","pval","var","stddev")
stat<-c("n Reps","Error Var","Genotypic Var", "Year Var", "Heritability","LSD","pval","var","stddev")
statistics<-data.frame(stat)

name<-colnames(data) 
name
name<-name[c(-1,-2)]
for(i in 1:length(name)){
  data5<-data[c(1,2,i+2)] 
  names(data5)[3]<-paste("trait")
  
  fm <- lmer(trait ~ (1|Accession) + (1|Year) + (1|Accession:Year), data=data5)
  
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
  
  varYear <- as.vector(varcorr$'Year')
  reps <- unique(data5$Year)
  nRep <- length(reps)
  LSD <- 1.96*sqrt(varErr) #calculate least significant difference; Z score for 95% confidence=1.96
  coefs <- data.frame(coef(summary(fm))) 
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))) 
  pval<-coefs[,4]
  var<- var(data5$trait,na.rm=TRUE)
  stddev<- sqrt(var)
  
  h2 <- varG/(varG + varErr/nRep) #broad sense heritability
  statis[,1] <- round(c(nRep,varErr,varG,varYear,h2,LSD,pval,var,stddev),9) 
  statis2 <- data.frame(statis,stat)
  colnames(statis2)[1]<-paste(name[i])
  statistics <- merge(statistics,statis2,by="stat")
}
write.table(statistics, file="TableS5_RatioData.csv", sep=",", row.names = FALSE)

#RESULTS: Plant phenotypes vary among genotypes####
cat("\014")
rm(list=ls()) 
ls() 
data16 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
##ANOVAs: Phenotypes differ between genotypes - Table S7-S8#### 
colnames(data16)
data17 = data16[,c(3,5:8,10:19)]
#Note repeat these steps for each trait
colnames(data17)
data17$x = data17$stalk_to_rootgrounding.cm
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
res.aov = aov(data17$tukey ~ Accession*PLOT_ID, data = data17)
summary(res.aov)
tukey.test = HSD.test(res.aov, trt = c("Accession","PLOT_ID"))
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
df = data16
df = df[,c(8,10:11,13:19)]
pca <- prcomp(df, scale. = TRUE)
var <- get_pca_var(pca)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)

###Table S9 and S10 Information####
summary(pca) #S9 Info
var <- get_pca_var(pca)
var$contrib #S10 Info

##Correlation of 2019 phenotypes (all reps inluded)-FigureS6#####
cor(df)
corrplot(cor(df), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")

##Heritability Info - Table S5#### 
data15 = data16
colnames(data15)
data15 = data15[,c(5,8,10:19)]
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
#write.table(statistics, file="TableS4_PhenotypeData.csv", sep=",", row.names = FALSE)




#RESULTS: Multiple phenotypes predict the contribution of brace roots to anchorage ####
##Phenotype correlations with the brace root contribution to anchorage & lodging - Table 1 ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("./Figure Files/Table1_10212021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(2:4,6:16)]
df=cor(df_cor)
df = as.data.frame(df)
write.csv(df, "Table1_Restults.csv", quote = F, row.names = T)




##Multiple Regression Modeling - The brace root contribution to anchorage#####
##Outcome = BRC, within year (2019) paired predictors + BRC
cat("\014")
rm(list=ls()) 
ls() 
getwd()
data <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data)
plant_data = data[,c(7,8,10:19)]
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
tibble::glimpse(plant_test)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(ratio_none_all_imu, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
  recipes::update_role(height_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_width_cm, new_role = "predictor") %>%
  recipes::update_role(spread_width_cm, new_role = "predictor") %>%
  recipes::update_role(single_root_width_cm, new_role = "predictor") %>%
  recipes::update_role(root_heightonstalk_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_to_rootgrounding_cm, new_role = "predictor") 
plant_recipe = recipes::prep(plant_recipe)
plant_head = recipes::juice(plant_recipe)
plant_train = recipes::bake(plant_recipe, plant_train)
plant_test = recipes::bake(plant_recipe, plant_test)
linear_lm = linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression") 
plant_workflow = workflows::workflow() %>% 
  workflows::add_recipe(plant_recipe) %>% 
  workflows::add_model(linear_lm) 
plant_fit = parsnip::fit(plant_workflow, plant_train) 
plant_pred = predict(plant_fit, plant_test)
plant_test = dplyr::bind_cols(plant_test, plant_pred) 
plant_metrics = yardstick::metrics(plant_test, truth = ratio_none_all_imu, estimate = .pred)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
plant_control = tune::control_grid(
  verbose = FALSE, 
  allow_par = TRUE, 
  extract = function(x){extract_model(x)},
  save_pred = TRUE 
)
plant_cv = tune::fit_resamples(plant_workflow, plant_folds, control = plant_control)
plant_metrics = tune::collect_metrics(plant_cv)  
plant_metrics
plant_predictions = tune::collect_predictions(plant_cv)
plant_cv$.extracts[[1]][[1]] 

#Outcome = BRC, BRC averaged within year (phenos averaged, BRC averaged)
cat("\014")
rm(list=ls()) 
ls() 
getwd()
data <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data)
data = data[,c(5,7,8,10:19)]
colnames(data)
MEANS = data[,c(1)]
MEANS = as.data.frame(MEANS)
colnames(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data) 
name 
name<-name[c(-1)]
for(i in 1:length(name)){
  data1 = data[c(1,i+1)]
  names(data1)[2] = paste("trait")
  data3 = data1 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data3 = as.data.frame(data3) 
  test = data3[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
data4 = MEANS
colnames(data4)
plant_data = data4[,c(2:13)]
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9)
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
tibble::glimpse(plant_test)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(ratio_none_all_imu, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
  recipes::update_role(height_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_width_cm, new_role = "predictor") %>%
  recipes::update_role(spread_width_cm, new_role = "predictor") %>%
  recipes::update_role(single_root_width_cm, new_role = "predictor") %>%
  recipes::update_role(root_heightonstalk_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_to_rootgrounding_cm, new_role = "predictor") 
plant_recipe = recipes::prep(plant_recipe)
plant_head = recipes::juice(plant_recipe)
plant_train = recipes::bake(plant_recipe, plant_train)
plant_test = recipes::bake(plant_recipe, plant_test)
linear_lm = linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression") 
plant_workflow = workflows::workflow() %>% 
  workflows::add_recipe(plant_recipe) %>% 
  workflows::add_model(linear_lm) 
plant_fit = parsnip::fit(plant_workflow, plant_train) 
plant_pred = predict(plant_fit, plant_test)
plant_test = dplyr::bind_cols(plant_test, plant_pred) 
plant_metrics = yardstick::metrics(plant_test, truth = ratio_none_all_imu, estimate = .pred)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
plant_control = tune::control_grid(
  verbose = FALSE, 
  allow_par = TRUE, 
  extract = function(x){extract_model(x)},
  save_pred = TRUE 
)
plant_cv = tune::fit_resamples(plant_workflow, plant_folds, control = plant_control)
plant_metrics = tune::collect_metrics(plant_cv)  
plant_metrics
plant_predictions = tune::collect_predictions(plant_cv)
plant_cv$.extracts[[1]][[1]] 

#Outcome = BRC, BRC averaged across years for each genotype (all predictors included)
cat("\014")
rm(list=ls()) 
ls() 
getwd()
data <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, na.strings = "NA")
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data)
data = data[,c(2,7)]
colnames(data)
MEANS = data[,c(1)]
MEANS = as.data.frame(MEANS)
colnames(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data) 
name 
name<-name[c(-1)]
for(i in 1:length(name)){
  data1 = data[c(1,i+1)]
  names(data1)[2] = paste("trait")
  data3 = data1 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data3 = as.data.frame(data3) 
  test = data3[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
data4 = MEANS
colnames(data2)
data2 = data2[,c(5,8,10:19)]
data5 = merge(data4, data2, by = "Accession")
colnames(data5)
plant_data = data5[,c(2:13)]
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
tibble::glimpse(plant_test)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(ratio_none_all_imu, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
  recipes::update_role(height_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_width_cm, new_role = "predictor") %>%
  recipes::update_role(spread_width_cm, new_role = "predictor") %>%
  recipes::update_role(single_root_width_cm, new_role = "predictor") %>%
  recipes::update_role(root_heightonstalk_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_to_rootgrounding_cm, new_role = "predictor") 
plant_recipe = recipes::prep(plant_recipe)
plant_head = recipes::juice(plant_recipe)
plant_train = recipes::bake(plant_recipe, plant_train)
plant_test = recipes::bake(plant_recipe, plant_test)
linear_lm = linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression") 
plant_workflow = workflows::workflow() %>% 
  workflows::add_recipe(plant_recipe) %>% 
  workflows::add_model(linear_lm) 
plant_fit = parsnip::fit(plant_workflow, plant_train) 
plant_pred = predict(plant_fit, plant_test)
plant_test = dplyr::bind_cols(plant_test, plant_pred) 
plant_metrics = yardstick::metrics(plant_test, truth = ratio_none_all_imu, estimate = .pred)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
plant_control = tune::control_grid(
  verbose = FALSE, 
  allow_par = TRUE, 
  extract = function(x){extract_model(x)},
  save_pred = TRUE 
)
plant_cv = tune::fit_resamples(plant_workflow, plant_folds, control = plant_control)
plant_metrics = tune::collect_metrics(plant_cv)  
plant_metrics
plant_predictions = tune::collect_predictions(plant_cv)
plant_cv$.extracts[[1]][[1]] 

#Outcome = BRC, BRC averaged across years for each genotype (predictors averaged across years)
cat("\014")
rm(list=ls()) 
ls() 
getwd()
data <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, na.strings = "NA")
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data)
data = data[,c(2,7)]
colnames(data)
MEANS = data[,c(1)]
MEANS = as.data.frame(MEANS)
colnames(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data) 
name 
name<-name[c(-1)]
for(i in 1:length(name)){
  data1 = data[c(1,i+1)]
  names(data1)[2] = paste("trait")
  data3 = data1 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data3 = as.data.frame(data3) 
  test = data3[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
data4 = MEANS
colnames(data2)
data2 = data2[,c(5,8,10:19)]

MEANS = data2[,c(1)]
MEANS = as.data.frame(MEANS)
colnames(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data2) 
name 
name<-name[c(-1)]
for(i in 1:length(name)){
  data1 = data2[c(1,i+1)]
  names(data1)[2] = paste("trait")
  data3 = data1 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data3 = as.data.frame(data3) 
  test = data3[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
data5 = MEANS
data5 = merge(data4, data5, by = "Accession")
colnames(data5)
plant_data = data5[,c(2:13)]
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9)
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
tibble::glimpse(plant_test)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(ratio_none_all_imu, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
  recipes::update_role(height_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_width_cm, new_role = "predictor") %>%
  recipes::update_role(spread_width_cm, new_role = "predictor") %>%
  recipes::update_role(single_root_width_cm, new_role = "predictor") %>%
  recipes::update_role(root_heightonstalk_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_to_rootgrounding_cm, new_role = "predictor") 
plant_recipe = recipes::prep(plant_recipe)
plant_head = recipes::juice(plant_recipe)
plant_train = recipes::bake(plant_recipe, plant_train)
plant_test = recipes::bake(plant_recipe, plant_test)
linear_lm = linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression") 
plant_workflow = workflows::workflow() %>% 
  workflows::add_recipe(plant_recipe) %>% 
  workflows::add_model(linear_lm) 
plant_fit = parsnip::fit(plant_workflow, plant_train) 
plant_pred = predict(plant_fit, plant_test)
plant_test = dplyr::bind_cols(plant_test, plant_pred) 
plant_metrics = yardstick::metrics(plant_test, truth = ratio_none_all_imu, estimate = .pred)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
plant_control = tune::control_grid(
  verbose = FALSE, 
  allow_par = TRUE, 
  extract = function(x){extract_model(x)},
  save_pred = TRUE 
)
plant_cv = tune::fit_resamples(plant_workflow, plant_folds, control = plant_control)
plant_metrics = tune::collect_metrics(plant_cv)  
plant_metrics
plant_predictions = tune::collect_predictions(plant_cv)
plant_cv$.extracts[[1]][[1]] 

##Random Forest Models####
##Outcome = BRC, within year (2019) paired predictors + BRC
cat("\014")
rm(list=ls()) 
ls() 
data12 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
str(data12)
i <- c(1:5)                         
data12[,i] <- apply(data12[ ,i], 2,          
                    function(x) as.character(x))
data12$Ratio..None.All....IMU = scale(data12$Ratio..None.All....IMU, center=TRUE, scale=TRUE)
data12$Ratio..None.All....IMU = as.numeric(data12$Ratio..None.All....IMU)
mean = mean(data12$Ratio..None.All....IMU)
sd = sd(data12$Ratio..None.All....IMU)
high = mean+sd
low = mean-sd
data12$Ratio_cat = data12$Ratio..None.All....IMU
data12$Ratio_cat = as.numeric(data12$Ratio_cat)
for (i in 1:length(data12$Ratio_cat)){
  if (data12$Ratio..None.All....IMU[i] <= low) {
    data12$Ratio_cat[i] = "low"
  } else if (data12$Ratio..None.All....IMU[i] < high) {
    data12$Ratio_cat[i] = "average"
  } else if (data12$Ratio..None.All....IMU[i] >= high) {
    data12$Ratio_cat[i] = "high"
  }
}
colnames(data12)
plant_data = data12[,c(20,8,10:19)]
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(ratio_cat ~ .) %>%
  recipes::step_string2factor(ratio_cat)
rf_tune = parsnip::rand_forest(trees = 100, mtry = tune(), min_n = tune()) %>%
  parsnip::set_engine("randomForest") %>%
  parsnip::set_mode("classification")
tune_wf = workflows::workflow() %>%
  workflows::add_recipe(rf_recipe) %>%
  workflows::add_model(rf_tune)
control = control_grid(extract = function(x) extract_model(x))
tune_res = tune::tune_grid(tune_wf, resample=plant_folds, grid=5, control=control)
tunes_metrics = tune::collect_metrics(tune_res)
best_params = tune::select_best(tune_res, "accuracy")
final_rf = tune::finalize_model(rf_tune, best_params)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_wf
final_res = final_wf %>%
  tune::last_fit(plant_split)
collect_metrics(final_res)
model = tune_res$.extracts[[1]]$.extracts[[1]]
model
randomForest::importance(model)

##Between year genotype; None/All averaged across years, all plants for pheno included
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, na.strings = "NA")
data1 = data %>%
  group_by(Accession) %>%
  summarise(Ratio_avg = mean(Ratio..None.All....IMU))
head(data1)
hist(data1$Ratio_avg)
data1$Ratio_avg2 = scale(data1$Ratio_avg, center=TRUE, scale=TRUE)
hist(data1$Ratio_avg2)
mean = mean(data1$Ratio_avg2)
sd = sd(data1$Ratio_avg2)
high = mean+sd
low = mean-sd
data1$Ratio_Category = data1$Ratio_avg2
for (i in 1:length(data1$Ratio_Category)){
  if (data1$Ratio_avg2[i] <= low) {
    data1$Ratio_Category[i] = "low"
  } else if (data1$Ratio_avg2[i] < high) {
    data1$Ratio_Category[i] = "average"
  } else if (data1$Ratio_avg2[i] >= high) {
    data1$Ratio_Category[i] = "high"
  }
}
data1 = data1[,c(1,2,4)]
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data2)
data2 = data2[,c(5,8:19)]
data2 = merge(data1, data2, by = "Accession")
str(data2)
colnames(data2)
plant_data = data2[,c(3:4,6:15)]
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(ratio_category ~ .) %>%
  recipes::step_string2factor(ratio_category)
rf_tune = parsnip::rand_forest(trees = 500, mtry = tune(), min_n = tune()) %>%
  parsnip::set_engine("randomForest") %>%
  parsnip::set_mode("classification")
tune_wf = workflows::workflow() %>%
  workflows::add_recipe(rf_recipe) %>%
  workflows::add_model(rf_tune)
control = control_grid(extract = function(x) extract_model(x))
tune_res = tune::tune_grid(tune_wf, resample=plant_folds, grid=5, control=control) 
tunes_metrics = tune::collect_metrics(tune_res)
best_params = tune::select_best(tune_res, "accuracy")
final_rf = tune::finalize_model(rf_tune, best_params)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_wf
final_res = final_wf %>%
  tune::last_fit(plant_split)
collect_metrics(final_res)
model = tune_res$.extracts[[1]]$.extracts[[1]]
model
randomForest::importance(model)

#RESULTS: Genotype determines root lodging susceptibility#####
##Figure S7 Prep####
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

##Heritability Info - Table S5#### 
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
write.table(statistics, file="TableS5_LodgingData.csv", sep=",", row.names = FALSE)
#RESULTS: The contribution of brace roots to anchorage is a good proxy for root lodging susceptibility####
##Multiple Regression Modeling - Lodging susceptibility #####
##Outcome = lodging susceptibility, lodging susceptibility averaged for each genotype (all 2019 plants included as predictors)
cat("\014")
rm(list=ls()) 
ls() 
getwd()
data <- read.csv(file = "LodgingData_2020.csv", header = TRUE, na.strings = "NA")
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
data$Accession = data$Genotype
colnames(data)
data = data[,c(13,10,11)]
colnames(data)
MEANS = data[,c(1)]
MEANS = as.data.frame(MEANS)
colnames(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data) 
name 
name<-name[c(-1)]
for(i in 1:length(name)){
  data1 = data[c(1,i+1)]
  names(data1)[2] = paste("trait")
  data3 = data1 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data3 = as.data.frame(data3) 
  test = data3[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
data4 = MEANS
data5 = merge(data4, data2, by = "Accession")
colnames(data5)
plant_data = data5[,c(2,10,12:21)]
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
tibble::glimpse(plant_test)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(percent_lodging, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
  recipes::update_role(height_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_width_cm, new_role = "predictor") %>%
  recipes::update_role(spread_width_cm, new_role = "predictor") %>%
  recipes::update_role(single_root_width_cm, new_role = "predictor") %>%
  recipes::update_role(root_heightonstalk_cm, new_role = "predictor") %>%
  recipes::update_role(stalk_to_rootgrounding_cm, new_role = "predictor") 
plant_recipe = recipes::prep(plant_recipe)
plant_head = recipes::juice(plant_recipe)
plant_train = recipes::bake(plant_recipe, plant_train)
plant_test = recipes::bake(plant_recipe, plant_test)
linear_lm = linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression") 
plant_workflow = workflows::workflow() %>% 
  workflows::add_recipe(plant_recipe) %>% 
  workflows::add_model(linear_lm) 
plant_fit = parsnip::fit(plant_workflow, plant_train) 
plant_pred = predict(plant_fit, plant_test)
plant_test = dplyr::bind_cols(plant_test, plant_pred) 
plant_metrics = yardstick::metrics(plant_test, truth = percent_lodging, estimate = .pred)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
plant_control = tune::control_grid(
  verbose = FALSE, 
  allow_par = TRUE, 
  extract = function(x){extract_model(x)},
  save_pred = TRUE 
)
plant_cv = tune::fit_resamples(plant_workflow, plant_folds, control = plant_control)
plant_metrics = tune::collect_metrics(plant_cv)  
plant_metrics
plant_predictions = tune::collect_predictions(plant_cv)
plant_cv$.extracts[[1]][[1]] 

##Random Forest - Lodging susceptibility ####
##Outcome = lodging susceptibility, lodging susceptibility averaged for each genotype (all 2019 plants included as predictors)
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv(file = "LodgingData_2020.csv", header = TRUE, na.strings = "NA")
hist(data$Percent.Lodging)
head(data)
colnames(data)
data1 = data %>%
  group_by(Genotype) %>%
  summarise(susceptibility_avg = mean(Percent.Lodging))
data1$susceptibility_avg2 = scale(data1$susceptibility_avg, center=TRUE, scale=TRUE)
data1$susceptibility_avg2 = as.numeric(data1$susceptibility_avg2)
mean = mean(data1$susceptibility_avg2)
sd = sd(data1$susceptibility_avg2)
high = mean+sd
low = mean-sd
data1$susceptibility_avg_category2 = data1$susceptibility_avg2
data1$susceptibility_avg_category2 = as.numeric(data1$susceptibility_avg_category2)
for (i in 1:length(data1$susceptibility_avg_category2)){
  if (data1$susceptibility_avg_category2[i] <= low) {
    data1$susceptibility_avg_category2[i] = "low"
  } else if (data1$susceptibility_avg_category2[i] < high) {
    data1$susceptibility_avg_category2[i] = "average"
  } else if (data1$susceptibility_avg_category2[i] >= high) {
    data1$susceptibility_avg_category2[i] = "high"
  }
}
data1[data1$Genotype=="2369", "susceptibility_avg_category2"] = "low"
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data2)
data1$Accession = data1$Genotype
colnames(data1)
data3 = data1[,c(5,2,4)]
data4 = merge(data3, data2, by = "Accession")
colnames(data4)
plant_data = data4[,c(3,10,12:21)]
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(susceptibility_avg_category2 ~ .) %>%
  recipes::step_string2factor(susceptibility_avg_category2)
rf_tune = parsnip::rand_forest(trees = 500, mtry = tune(), min_n = tune()) %>%
  parsnip::set_engine("randomForest") %>%
  parsnip::set_mode("classification")
tune_wf = workflows::workflow() %>%
  workflows::add_recipe(rf_recipe) %>%
  workflows::add_model(rf_tune)
control = control_grid(extract = function(x) extract_model(x))
tune_res = tune::tune_grid(tune_wf, resample=plant_folds, grid=5, control=control)
tunes_metrics = tune::collect_metrics(tune_res) 
best_params = tune::select_best(tune_res, "accuracy")
final_rf = tune::finalize_model(rf_tune, best_params)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_wf 
final_res = final_wf %>%
  tune::last_fit(plant_split) 
collect_metrics(final_res)
model = tune_res$.extracts[[1]]$.extracts[[1]]
model
randomForest::importance(model)


##Random Forest - Lodging susceptibility w/ the brace root contribution to anchorage as a predictor ####
##Outcome = lodging susceptibility, lodging susceptibility averaged for each genotype (all 2019 plants included as predictors + associated 2019 brace root contribution to anchorage)
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv(file = "LodgingData_2020.csv", header = TRUE, na.strings = "NA")
hist(data$Percent.Lodging)
head(data)
colnames(data)
data1 = data %>%
  group_by(Genotype) %>%
  summarise(susceptibility_avg = mean(Percent.Lodging))
data1$susceptibility_avg2 = scale(data1$susceptibility_avg, center=TRUE, scale=TRUE)
data1$susceptibility_avg2 = as.numeric(data1$susceptibility_avg2)
mean = mean(data1$susceptibility_avg2)
sd = sd(data1$susceptibility_avg2)
high = mean+sd
low = mean-sd
data1$susceptibility_avg_category2 = data1$susceptibility_avg2
data1$susceptibility_avg_category2 = as.numeric(data1$susceptibility_avg_category2)
for (i in 1:length(data1$susceptibility_avg_category2)){
  if (data1$susceptibility_avg_category2[i] <= low) {
    data1$susceptibility_avg_category2[i] = "low"
  } else if (data1$susceptibility_avg_category2[i] < high) {
    data1$susceptibility_avg_category2[i] = "average"
  } else if (data1$susceptibility_avg_category2[i] >= high) {
    data1$susceptibility_avg_category2[i] = "high"
  }
}
data1[data1$Genotype=="2369", "susceptibility_avg_category2"] = "low"
#Merge Files 
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data2)
data1$Accession = data1$Genotype
colnames(data1)
data3 = data1[,c(5,2,4)]
data4 = merge(data3, data2, by = "Accession")
colnames(data4)
plant_data = data4[,c(3,9:10,12:21)]
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
tibble::glimpse(plant_test)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(susceptibility_avg_category2 ~ .) %>%
  recipes::step_string2factor(susceptibility_avg_category2)
rf_tune = parsnip::rand_forest(trees = 500, mtry = tune(), min_n = tune()) %>%
  parsnip::set_engine("randomForest") %>%
  parsnip::set_mode("classification")
tune_wf = workflows::workflow() %>%
  workflows::add_recipe(rf_recipe) %>%
  workflows::add_model(rf_tune)
control = control_grid(extract = function(x) extract_model(x))
tune_res = tune::tune_grid(tune_wf, resample=plant_folds, grid=5, control=control)
tunes_metrics = tune::collect_metrics(tune_res)
best_params = tune::select_best(tune_res, "accuracy")
final_rf = tune::finalize_model(rf_tune, best_params)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_wf 
final_res = final_wf %>%
  tune::last_fit(plant_split)
collect_metrics(final_res)
model = tune_res$.extracts[[1]]$.extracts[[1]]
model
randomForest::importance(model)
