#Inclusion of plant biomechanics improves the prediction accuracy of root lodging
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/")

#Load Library
library(dplyr)
library(corrplot)
library(tidyverse)
library(tidymodels)
library(ranger)
library(kknn)
library(janitor)

#Correlation - Genotypic level biomechanics and lodging (Prep for Figure 5A File) ####
cat("\014")
rm(list=ls()) 
ls() 
#Merge 2018 and 2019 biomechancis data
data <- read.csv(file = "Inbred_Subpop_2Years_EI_Processed.csv", header = TRUE, na.strings = "NA")
data2 <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", header = TRUE, na.strings = "NA")
colnames(data)
data$ID_Year_Accession_PlotID_PlantNumber = paste(data$Year, data$Accession, data$Plot.ID, data$Plant.Number, sep="_")
data2$ID_Year_Accession_PlotID_PlantNumber = paste(data2$Year, data2$Accession, data2$Plot.ID, data2$Plant.Number, sep="_")
data = data[,c(11,8)]
data2 = data2[,c(10,1,2,3,4,7,8)]
data3 = merge(data, data2, by = "ID_Year_Accession_PlotID_PlantNumber")
colnames(data3)
data3 = data3[,c(4,2,7,8)]
#Calculate genotype average across both years and  plots
MEANS = data3[,c(1)]
MEANS = as.data.frame(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data3)
name
name<-name[c(-1)]
for(i in 1:length(name)){
  data4 = data3[c(1,i+1)]
  names(data4)[2] = paste("trait")
  data5 = data4 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data5 = as.data.frame(data5) 
  test = data5[,c(1,2)]
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
head(MEANS)
biomech = MEANS
#Merge with lodging data
data6 <- read.csv(file = "LodgingData_2020_long_02162021.csv", header = TRUE, na.strings = "NA")
data6 = na.omit(data6)
#Calculate means of lodging data across all plots and merge with mechanics means
data6$Accession = data6$Genotype
colnames(data6)
data6 = data6[,c(9,6,7)]
MEANS = data6[,c(1)]
MEANS = as.data.frame(MEANS)
MEANS<-unique(MEANS[,c(1)])
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data6)
name
name<-name[c(-1)]
for(i in 1:length(name)){
  data7 = data6[c(1,i+1)]
  names(data7)[2] = paste("trait")
  data8 = data7 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data8 = as.data.frame(data8) 
  test = data8[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
lodging = MEANS
data9 = merge(biomech, lodging, by = "Accession")
write.csv(data9, file = "Figure5A_03182021.csv", quote = F, row.names = F)
#Correlations
df_cor = data9[,c(2,3,5,6)]
cor(df_cor)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")


###Random Forest - Lodging susceptibility ####
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
plant_data = data4[,c(3,8:21)]
plant_data = data4[,c(3,8:15,17:21)] #Run this when excluding height as a predictor
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
tibble::glimpse(plant_test)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(susceptibility_avg_category2 ~ .) %>%
  recipes::step_string2factor(susceptibility_avg_category2)
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


###Random Forest - Lodging severity ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv(file = "LodgingData_2020.csv", header = TRUE, na.strings = "NA")
data1 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data)
colnames(data1)
data$Accession = data$Genotype
data2 = data %>%
  group_by(Accession) %>%
  summarise(severity_avg = mean(Lodging.Index))
data2$severity_avg2 = scale(data2$severity_avg, center=TRUE, scale=TRUE)
data2$severity_avg2 = as.numeric(data2$severity_avg2)
mean = mean(data2$severity_avg2)
sd = sd(data2$severity_avg2)
high = mean+sd
low = mean-sd
data2$severity_avg_cat = data2$severity_avg2
for (i in 1:length(data2$severity_avg_cat)){
  if (data2$severity_avg_cat[i] <= low) {
    data2$severity_avg_cat[i] = "low"
  } else if (data2$severity_avg_cat[i] < high) {
    data2$severity_avg_cat[i] = "average"
  } else if (data2$severity_avg_cat[i] >= high) {
    data2$severity_avg_cat[i] = "high"
  }
}
data2 = data2[,c(1,2,4)]
data3 = merge(data1, data2, by = "Accession")
colnames(data3)
data4 = data3[,c(1,21,6:19)]
plant_data = data4[,c(2:16)]
plant_data = data4[,c(2:10,12:16)] #Run this when excluding height as a predictor
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
rf_recipe = head(training(plant_split)) %>%
  recipes::recipe(severity_avg_cat ~ .) %>%
  recipes::step_string2factor(severity_avg_cat)
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


