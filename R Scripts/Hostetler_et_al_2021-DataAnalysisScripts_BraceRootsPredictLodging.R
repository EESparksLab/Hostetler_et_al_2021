#Brace root phenotypes predict root lodging
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler et al., 2021/R Input Files/")

#Library Load
library(arsenal)
library(lme4)
library(car)
library(agricolae)
library(dplyr)
library(lsmeans)
library(tidyverse)
library(vegan)
library(ggplot2)
library(goeveg)
library(PerformanceAnalytics)
library(rcompanion)
library(corrplot)
library(tidyverse)
library(tidymodels)
library(ranger)
library(kknn)
library(janitor)

###Process and Prep for Table 1####
#Process EI and None/All Ratio Data for other analysis; remove genotypes not included in both years; remove EI and Ratio < 0
cat("\014")
rm(list=ls()) 
ls() 
#Remove individuals that were not tested in both years
#EI
MzHyb <- read.csv(file = "Inbred_Subpop_2Years_EI_All.csv", header = TRUE, na.strings = "NA")
MzHyb_2018 = subset(MzHyb, Year=="2018")
MzHyb_2019 = subset(MzHyb, Year=="2019")
summary = summary(comparedf(MzHyb_2018, MzHyb_2019, by="Accession"))
notShared = summary$obs.table
MzHyb1 = MzHyb[!MzHyb$Accession %in% notShared$Accession,]
MzHyb1 = MzHyb1[MzHyb1$Raw.Slope.IMU..N.m...Calculated.from.Load.and.Unload. > 0,]
MzHyb1 %>% count(Accession)
MzHyb1 = subset(MzHyb1, select=-c(Geographic.Location, Species, Load.Cell.Height, Raw.Data.Label, X, Raw.Slope.IMU..N.m...Calculated.from.Load.and.Unload.,X.1,Additional.Notes,Publication.s.))
MzHyb1 = subset(MzHyb1, Brace.Root.Status=="A") #Subset to only the plants with brace roots intact (A)
MzHyb1 = MzHyb1[,c(1:9)]
str(MzHyb1)
MzHyb1$Plot.ID = as.factor(MzHyb1$Plot.ID)
MzHyb1$Year = as.factor(MzHyb1$Year)
MzHyb1$ID = paste(MzHyb1$Accession, MzHyb1$Plot.ID, sep="_")
write.csv(MzHyb1, file = "Inbred_Subpop_2Years_EI_Processed.csv", quote = F, row.names = F)

#Ratio None/All
MzHyb_BR <- read.csv(file = "Inbred_Subpop_2Years_Broot_Ratio.csv", header = TRUE, na.strings = "NA")
MzHyb_2018_BR = subset(MzHyb_BR, Year=="2018")
MzHyb_2019_BR = subset(MzHyb_BR, Year=="2019")
summary_BR = summary(comparedf(MzHyb_2018_BR, MzHyb_2019_BR, by="Accession"))
notShared_BR = summary_BR$obs.table
MzHyb1_BR = MzHyb_BR[!MzHyb_BR$Accession %in% notShared_BR$Accession,]
MzHyb1_BR = MzHyb1_BR[MzHyb1_BR$Ratio..None.All....IMU > 0,]
MzHyb1_BR %>% count(Accession)
MzHyb1_BR = subset(MzHyb1_BR, select=-c(Geographic.Location,Species,Load.Cell.Height,Raw.Data.Label,X,X.1,Additional.Notes,Publication.s.))
str(MzHyb1_BR)
MzHyb1_BR$Plot.ID = as.factor(MzHyb1_BR$Plot.ID)
MzHyb1_BR$Year = as.factor(MzHyb1_BR$Year)
MzHyb1_BR$ID = paste(MzHyb1_BR$Accession, MzHyb1_BR$Plot.ID, sep="_")
MzHyb1_BR$ID = as.factor(MzHyb1_BR$ID)
write.csv(MzHyb1_BR, file = "Inbred_Subpop_2Years_Broot_Ratio_Processed.csv", quote = F, row.names = F)

#2018 and 2019 biomechanics data, merge and calculate genotype average within a year
data = MzHyb1
data2 = MzHyb1_BR
data %>% count(Accession)
data2 %>% count(Accession)
colnames(data)
data$ID_Year_Accession_PlotID_PlantNumber = paste(data$Year, data$Accession, data$Plot.ID, data$Plant.Number, sep="_")
data = data[,c(11,8)]
colnames(data2)
data2$ID_Year_Accession_PlotID_PlantNumber = paste(data2$Year, data2$Accession, data2$Plot.ID, data2$Plant.Number, sep="_")
data2 = data2[,c(10,1,2,3,4,7,8)]
data3 = merge(data, data2, by = "ID_Year_Accession_PlotID_PlantNumber")
data3 %>% count(Accession)
colnames(data3)
data3 = data3[,c(4,3,2,7)]
#Calculate means for within a genotype and year 
colnames(data3)
MEANS = data3[,c(1,2)]
MEANS<-unique(MEANS[,c(1,2)])  
name<-colnames(data3) 
name 
name<-name[c(-1,-2)]
for(i in 1:length(name)){
  data4 = data3[c(1,2,i+2)]
  names(data4)[3] = paste("trait")
  data5 = data4 %>%
    group_by(Accession, Year) %>%
    summarise(trait_avg = mean(trait))
  data5 = as.data.frame(data5) 
  test = data5[,c(1,2,3)]
  names(test)[3] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession","Year"))
}
head(MEANS)
biomech = MEANS
colnames(biomech)
biomech2 = pivot_wider(data=biomech,
                       id=Accession,
                       names_from = Year,
                       values_from = c("Flexural.Rigidity..EI....IMU", "Ratio..None.All....IMU"))

#Merge with lodging data, calculate means of lodging data across all plots and merge with mechanics means
data6 <- read.csv(file = "LodgingData_2020_long_02162021.csv", header = TRUE, na.strings = "NA")
data6 = na.omit(data6)
head(data6)
#Calculate means for within a genotype 
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
data9 = merge(biomech2, lodging, by = "Accession")
data9 %>% count(Accession)

#Merge with 2019 phenotypes and calculate genotypes means for 2019 phenotypes
data10 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data10)
data10 = data10[,c(5,2,3,4,8:19)]
colnames(data10)
data10 = data10[,c(1,5:16)]
MEANS = data10[,c(1)]
MEANS = as.data.frame(MEANS)
MEANS<-unique(MEANS[,c(1)]) 
MEANS = as.data.frame(MEANS)
MEANS$Accession = MEANS$MEANS
MEANS$MEANS = NULL
MEANS = as.data.frame(MEANS)
name<-colnames(data10) 
name 
name<-name[c(-1)] 
for(i in 1:length(name)){
  data11 = data10[c(1,i+1)]
  names(data11)[2] = paste("trait")
  data12 = data11 %>%
    group_by(Accession) %>%
    summarise(trait_avg = mean(trait))
  data12 = as.data.frame(data12) 
  test = data12[,c(1,2)] 
  names(test)[2] = paste(name[i])
  MEANS = merge(MEANS, test, by.x =c("Accession"))
}
pheno = MEANS
data9 = merge(data9, pheno, by = "Accession")
data9 %>% count(Accession)
write.csv(data9, file = "Table1_03292021.csv", quote = F, row.names = F)



#Multiple Regression Modeling #####
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
data6 = data5[,c(2:3,10:21)]
#Regression for lodging susceptibility ####
plant_data = data6[,c(1,3:8,10:14)] #Not including height in regression model
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
  recipes::update_role(num_whorls, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
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

#Regression for lodging severity 
colnames(data6)
plant_data = data6[,c(2:8,10:14)] #Not including height in regression model
plant_data = janitor::clean_names(plant_data)
colnames(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_recipe = head(plant_train) %>%
  recipes::recipe() %>%
  recipes::update_role(lodging_index, new_role = "outcome") %>%
  recipes::update_role(brace_root_whorls_in_the_soil, new_role = "predictor") %>%
  recipes::update_role(num_whorls, new_role = "predictor") %>%
  recipes::update_role(w1, new_role = "predictor") %>%
  recipes::update_role(w2, new_role = "predictor") %>%
  recipes::update_role(w3, new_role = "predictor") %>%
  recipes::update_role(root_angle, new_role = "predictor") %>%
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
plant_metrics = yardstick::metrics(plant_test, truth = lodging_index, estimate = .pred)
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
data2 <- read.csv(file = "Averaged2019_Pheno_and_Combined_Mech_Data_04132021.csv", header = TRUE, na.strings = "NA")
colnames(data2)
data1$Accession = data1$Genotype
colnames(data1)
data3 = data1[,c(5,2,4)]
data4 = merge(data3, data2, by = "Accession")
colnames(data4)
plant_data = data4[,c(3,10:21)]
plant_data = data4[,c(3,10:15,17:21)] #Use this when rerunning model without height included
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
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
plant_data = data3[,c(21,8:19)]
plant_data = data3[,c(21,8:13,15:19)] #When rerunning and excluding height
plant_data = janitor::clean_names(plant_data)
plant_split = rsample::initial_split(plant_data, prop = 0.9, strata = "brace_root_whorls_in_the_soil")
plant_test = rsample::testing(plant_split)
plant_train = rsample::training(plant_split)
plant_folds = rsample::vfold_cv(training(plant_split), v=5, repeats = 5)
tibble::glimpse(plant_test)
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



