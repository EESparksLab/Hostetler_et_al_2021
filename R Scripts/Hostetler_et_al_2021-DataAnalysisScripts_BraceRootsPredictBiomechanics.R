#Brace root phenotypes predict plant biomechanics
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/")

#Library Load
library(car)
library(lme4)
library(tidyverse)
library(tidymodels)
library(ranger)
library(kknn)
library(janitor)
library(agricolae)

###ANOVA on biomechanics ####
#EI
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "Inbred_Subpop_2Years_EI_Processed.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
colnames(data)
data = data[,c(1:3,8)]
attach(data)
lm_flex <- lm(Flexural.Rigidity..EI....IMU ~ Accession*Year)
summary(lm_flex)
anova(lm_flex)
lm_flex_aov=aov(lm_flex) 
HSD.test(lm_flex_aov, trt = c("Accession", "Year"), console = TRUE)
par(mfrow=c(2,2))
plot(lm_flex)
par(mfrow=c(1,1))
#transformation data to meet assumptions
par(mfrow=c(2,1))
plot(data$Flexural.Rigidity..EI....IMU)
hist(data$Flexural.Rigidity..EI....IMU)
par(mfrow=c(2,2))
hist(data$Flexural.Rigidity..EI....IMU)
data$trans_Flexural.Rigidity..EI....IMU_2 = sqrt(data$Flexural.Rigidity..EI....IMU)
hist(data$trans_Flexural.Rigidity..EI....IMU_2) #best
data$trans_Flexural.Rigidity..EI....IMU_3 = log10(data$Flexural.Rigidity..EI....IMU)
hist(data$trans_Flexural.Rigidity..EI....IMU_3)
data$trans_Flexural.Rigidity..EI....IMU_4 = (1/(data$Flexural.Rigidity..EI....IMU))
hist(data$trans_Flexural.Rigidity..EI....IMU_4)
data$trans_Flexural.Rigidity..EI....IMU_3 = NULL
data$trans_Flexural.Rigidity..EI....IMU_4 = NULL
lm_flex2 <- lm(data$trans_Flexural.Rigidity..EI....IMU_2 ~ Accession * Year) #create a linear model and look at how factors predict phenotype
anova(lm_flex2)
lm_flex_aov2=aov(lm_flex2) 
HSD.test(lm_flex_aov2, trt = c("Accession", "Year"), console = TRUE)
detach(data)

#Ratio
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
data$Plot.ID = as.factor(data$Plot.ID)
data$Year = as.factor(data$Year)
attach(data)
par(mfrow=c(1,1))
lm_ratio <- lm(Ratio..None.All....IMU ~ Accession + Year)
lm_ratio_aov=anova(lm_ratio) 
lm_ratio_aov
lm_ratio_aov=aov(lm_ratio) 
Ratio_Tukey1=TukeyHSD(lm_ratio_aov, which = "Accession:Year") 
TukeyOutput1=Ratio_Tukey1$`Accession:Year`
TukeyOutput1 = as.data.frame(TukeyOutput1)
HSD.test(lm_ratio_aov, trt = c("Accession", "Year"), console = TRUE)
par(mfrow=c(2,2))
plot(lm_ratio)
par(mfrow=c(1,1))
plot(Ratio..None.All....IMU)
hist(Ratio..None.All....IMU)
hist(lm_ratio$residuals)
outlierTest(lm_ratio)
detach(data)

#Heritability Calculations #####
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
write.table(statistics, file="TableS4_RatioData.csv", sep=",", row.names = FALSE)

#EI
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv(file = "Inbred_Subpop_2Years_EI_Processed.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
data = subset(data, Brace.Root.Status == "A")
head(data)
colnames(data)
data = data[,c(1,2,8,9)]
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
  
  #Fit the model here
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
  LSD <- 1.96*sqrt(varErr) 
  coefs <- data.frame(coef(summary(fm))) 
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))) 
  pval<-coefs[,4]
  var<- var(data5$trait,na.rm=TRUE) 
  stddev<- sqrt(var) 

  h2 <- varG/(varG + varErr/nRep) 
  statis[,1] <- round(c(nRep,varErr,varG,varYear,h2,LSD,pval,var,stddev),9) 
  statis2 <- data.frame(statis,stat)
  colnames(statis2)[1]<-paste(name[i])
  statistics <- merge(statistics,statis2,by="stat")
}
write.table(statistics, file="TableS4_EIData.csv", sep=",", row.names = FALSE)


#Random Forest Models ####
#Flexural Stiffness - Genotypic averages
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv(file = "Inbred_Subpop_2Years_EI_Processed.csv", header = TRUE, na.strings = "NA")
str(data)
colnames(data)
data3 = data %>%
  group_by(Accession) %>%
  summarise(EI_avg = mean(Flexural.Rigidity..EI....IMU))
hist(data3$EI_avg)
data3$EI_avg = scale(data3$EI_avg, center=TRUE, scale=TRUE)
mean = mean(data3$EI_avg)
sd = sd(data3$EI_avg)
high = mean+sd
low = mean-sd
data3$EI_Category = data3$EI_avg
for (i in 1:length(data3$EI_Category)){
  if (data3$EI_avg[i] <= low) {
    data3$EI_Category[i] = "low"
  } else if (data3$EI_avg[i] < high) {
    data3$EI_Category[i] = "average"
  } else if (data3$EI_avg[i] >= high) {
    data3$EI_Category[i] = "high"
  }
}
data4 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data4)
data4 = data4[,c(5,8:19)]
data4 = merge(data3, data4, by = "Accession")
plant_data = data4[,c(3:15)]
plant_data = data4[,c(3:9,11:15)] #Run this line if you are excluding height as a predictor
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(ei_category ~ .) %>%
  recipes::step_string2factor(ei_category)
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

#Ratio None/All - Genotypic averages
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
plant_data = data2[,c(3:15)]
plant_data = data2[,c(3:9,11:15)] #Run this when excluding height as a predictor from the model
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(ratio_category ~ .) %>%
  recipes::step_string2factor(ratio_category)
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

#Flexural Stiffness - Paired Data
cat("\014")
rm(list=ls()) 
ls() 
data12 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
str(data12)
i <- c(1:5)                          
data12[,i] <- apply(data12[ ,i], 2,        
                    function(x) as.character(x))

colnames(data12)
data12$Flexural.Rigidity..EI....IMU = scale(data12$Flexural.Rigidity..EI....IMU, center=TRUE, scale=TRUE)
data12$Flexural.Rigidity..EI....IMU = as.numeric(data12$Flexural.Rigidity..EI....IMU)
mean = mean(data12$Flexural.Rigidity..EI....IMU)
sd = sd(data12$Flexural.Rigidity..EI....IMU)
high = mean+sd
low = mean-sd
data12$EI_cat = data12$Flexural.Rigidity..EI....IMU
data12$EI_cat = as.numeric(data12$EI_cat)
for (i in 1:length(data12$EI_cat)){
  if (data12$Flexural.Rigidity..EI....IMU[i] <= low) {
    data12$EI_cat[i] = "low"
  } else if (data12$Flexural.Rigidity..EI....IMU[i] < high) {
    data12$EI_cat[i] = "average"
  } else if (data12$Flexural.Rigidity..EI....IMU[i] >= high) {
    data12$EI_cat[i] = "high"
  }
}

plant_data = data12[,c(20,8:19)]
plant_data = data12[,c(20,8:13,15:19)] #Run this when excluding height as a predictor 
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(ei_cat ~ .) %>%
  recipes::step_string2factor(ei_cat)
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


#Ratio None/All - Paired Data
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
plant_data = data12[,c(20,8:19)]
plant_data = data12[,c(20,8:13,15:19)] #Run these parameters when excluding height
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



