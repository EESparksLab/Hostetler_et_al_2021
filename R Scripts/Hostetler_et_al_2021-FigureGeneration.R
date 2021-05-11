###Input File Info####
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/Figure Files/")

####Libary Load####
library(ggplot2)
library(arsenal)
library(dplyr)
library(reshape2)
library(tidyr)
library(factoextra)
library(ggfortify)
library(corrplot)
library(lattice)
library(cowplot)

###Figure 2####
#Figure 2A-B
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv("Figure2A-B_04122021.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")

Fig2A=ggplot(lodging, aes(x=Lodging.Index)) +
  geom_histogram(bins=5, color="black")+
  xlab("Lodging severity") + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        plot.title=element_text(size=12, vjust=-3), 
        axis.text=element_text(size=12), 
        axis.title = element_text(size=12), 
        axis.title.y= element_text(vjust=2.5), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.5,"cm"),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 12, colour = "blue"),
        strip.text.y = element_text(size = 12, colour = "black"))+
  scale_y_continuous(name="Number of plants", breaks=seq(0,850,100))
Fig2A
pdf("Figure2A.pdf", width = 11, height = 8.5)
plot(Fig2A)
dev.off()

lodging$Lodged2 = as.numeric(lodging$Lodged2)
Fig2B=ggplot(lodging, aes(x=Lodged2)) +
  geom_histogram(bins=2,color="black")+
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        plot.title=element_text(size=12, vjust=-3), 
        axis.text=element_text(size=12), 
        axis.title = element_text(size=12), 
        axis.title.y= element_text(vjust=2.5), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.5,"cm"),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 12, colour = "blue"),
        strip.text.y = element_text(size = 12, colour = "black"))+
  scale_y_continuous(name="Number of plants", breaks=seq(0,900,100))+
  scale_x_continuous(name="Lodging susceptibility", breaks=seq(0,2,1))
Fig2B
pdf("Figure2B.pdf", width = 11, height = 8.5)
plot(Fig2B)
dev.off()

#Figure 2C-E
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv("Figure2C-E.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
accessions <- read.csv("Unique_55Genotypes_Accession.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
colnames(lodging)
colnames(accessions)
accessions$Genotype = accessions$x
accessions$x = NULL
lodging$Replicate = as.character(lodging$Replicate)
lodging$Plot = as.character(lodging$Plot)
summary = summary(comparedf(lodging, accessions, by="Genotype"))
notShared = summary$obs.table
lodging = lodging[!lodging$Genotype %in% notShared$Genotype,]
head(lodging)
lodging$Accession = lodging$Genotype

data_2C = lodging %>%
  dplyr::filter(Accession %in% c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71","A632","Mo17","PHK46","PHHB4","Ki3","CML322","W64A","Oh43"))
data_2C %>% count(Accession)

data_2D = lodging %>%
  dplyr::filter(Accession %in% c("CML52","SA24","CML10","Tzi9","Tzi8","Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103","Tx303","B97","R84"))
data_2D %>% count(Accession)

data_2E = lodging %>%
  dplyr::filter(Accession %in% c("R78","PHB47","Ky21","CML258","B73","CML341","Shoe Peg","LH123Ht","2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
data_2E %>% count(Accession)

#Figure 2C
colnames(data_2C)
lodging3 = data_2C[,c(3,13,8,9,12)]
lodging6 = lodging3[,c(1,3,4,5)]
lodging3_long = melt(lodging6, id = c("Plot","Lodging.Index.Cat"))
lodging3_long$Lodging.Index.Cat = as.character(lodging3_long$Lodging.Index.Cat)
lodging3_long$Lodging.Index.Cat[lodging3_long$variable == "Total.Plants"] = "Total"

lodging3$Accession = factor(lodging3$Accession,levels = c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71","A632","Mo17","PHK46","PHHB4","Ki3","CML322","W64A","Oh43"))
plotC=ggplot(data=lodging3_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  facet_grid(. ~ lodging3$Accession, scales="free", space="free")
plotC
pdf("Figure2C.pdf", width = 11, height = 8.5)
plot(plotC)
dev.off()
#Figure 2D
colnames(data_2D)
lodging4 = data_2D[,c(3,13,8,9,12)]
lodging7 = lodging4[,c(1,3,4,5)]
lodging4_long = melt(lodging7, id = c("Plot","Lodging.Index.Cat"))
lodging4_long$Lodging.Index.Cat = as.character(lodging4_long$Lodging.Index.Cat)
lodging4_long$Lodging.Index.Cat[lodging4_long$variable == "Total.Plants"] = "Total"
lodging4$Accession = factor(lodging4$Accession,levels = c("CML52","SA24","CML10","Tzi9","Tzi8","Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103","Tx303","B97","R84"))
plotD=ggplot(data=lodging4_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  facet_grid(. ~ lodging4$Accession, scales="free", space="free")
plotD
pdf("Figure2D.pdf", width = 11, height = 8.5)
plot(plotD)
dev.off()

#Figure 2E
colnames(data_2E)
lodging2 = data_2E[,c(3,13,8,9,12)]
colnames(lodging2)
lodging5 = lodging2[,c(1,3,4,5)]
lodging2_long = melt(lodging5, id = c("Plot","Lodging.Index.Cat"))
lodging2_long$Lodging.Index.Cat = as.character(lodging2_long$Lodging.Index.Cat)
lodging2_long$Lodging.Index.Cat[lodging2_long$variable == "Total.Plants"] = "Total"
lodging2$Accession = factor(lodging2$Accession,levels = c("R78","PHB47","Ky21","CML258","B73","CML341","Shoe Peg","LH123Ht","2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
plotE=ggplot(data=lodging2_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  facet_grid(. ~ lodging2$Accession, scales="free", space="free")
plotE
pdf("Figure2E.pdf", width = 11, height = 8.5)
plot(plotE)
dev.off()

####Figure 3####
cat("\014")
rm(list=ls()) 
ls() 
df <- read.csv("Figure3_03302021.csv", header = TRUE, na.strings = "NA")
pca <- prcomp(df, scale. = TRUE)
var <- get_pca_var(pca)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
pdf("Figure3.pdf", width = 8.5, height = 11)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
dev.off()

####Figure 4####
#Panel A
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure4A_03302021.csv", header = TRUE, na.strings = "NA")
data$Year = as.factor(data$Year)
data$PLOT_ID2 = as.factor(data$PLOT_ID2)
data$PLOT_ID3 = as.factor(data$PLOT_ID3)
data$Accession = factor(data$Accession,levels = c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71",
                                                   "A632","Mo17","PHK46","PHHB4","Ki3","CML322","CML52","SA24","CML10","Tzi9","Tzi8",
                                                   "Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103",
                                                   "Tx303","B97","R84","R78","PHB47","Ky21","CML258","B73","CML341","Shoe Peg","LH123Ht",
                                                   "2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
Fig4A=ggplot(data=data,aes(x=Accession, y=Flexural.Rigidity..EI....IMU, color=Year))+
  geom_boxplot(size = 0.25, color = "black", fill="gray80") + 
  geom_point(size = .50) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("red4","darkslategray4"))+
  xlab("Genotype") + 
  ylab("EI")
Fig4A
pdf("Figure4A.pdf", width = 8, height = 5)
plot(Fig4A)
dev.off()

#Panel B
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure4B_03302021.csv", header = TRUE, na.strings = "NA")
data$Year = as.factor(data$Year)
data$PLOT_ID2 = as.factor(data$PLOT_ID2)
data$PLOT_ID3 = as.factor(data$PLOT_ID3)
data$Accession = factor(data$Accession,levels = c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71",
                                                  "A632","Mo17","PHK46","PHHB4","Ki3","CML322","CML52","SA24","CML10","Tzi9","Tzi8",
                                                  "Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103",
                                                  "Tx303","B97","R84","R78","PHB47","Ky21","CML258","B73","CML341","Shoe Peg","LH123Ht",
                                                  "2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
Fig4B=ggplot(data=data,aes(x=Accession, y=Ratio..None.All....IMU, color=Year))+
  geom_boxplot(size = 0.25, color = "black",fill="gray80") + 
  geom_point(size = .75) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("red4","darkslategray4"))+
  xlab("Genotype") + 
  ylab("Ratio=None/All")
Fig4B
pdf("Figure4B.pdf", width = 8, height = 5)
plot(Fig4B)
dev.off()

####Figure 5######
#Figure 5A
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("Figure5A_03182021.csv", header = TRUE, na.strings = "NA")
cor(data)
corrplot(cor(data), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
pdf("BiomechsLodging_GenoMeans_Corr.pdf", width = 11, height = 8.5)
corrplot(cor(data), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
dev.off()

#Figure 5B
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("Figure5B_03302021.csv", header = TRUE, na.strings = "NA")
head(data)
df1 = subset(data, ModelNumber=="Model10")
df2 = subset(data, ModelNumber=="Model6")
df3 = subset(data, ModelNumber=="Model16")
df4 = subset(data, ModelNumber=="Model14")
data1 = rbind(df1, df2, df3, df4)
data1$Predictor = factor(data1$Predictor,levels = c("flexural_rigidity_ei_imu","ratio_none_all_imu","height_cm","stalk_width_cm","single_root_width_cm","w1","root_heightonstalk_cm","root_angle","stalk_to_rootgrounding_cm","spread_width_cm","w2","brace_root_whorls_in_the_soil","num_whorls","w3"))
data1$ModelNumber = factor(data1$ModelNumber, levels = c("Model10","Model6","Model16","Model14"))
Fig5B=ggplot(data1) +
  geom_segment( aes(x=Outcome, xend=Outcome, y=0, yend=LevelOfImportance), color="grey") +
  geom_point( aes(x=Outcome, y=LevelOfImportance, color=Outcome), size=3 ) +
  coord_flip()+
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  theme_classic()+
  xlab("") +
  ylab("Predictor rank") +
  scale_y_continuous(breaks=seq(0,11,1))+
  facet_wrap(~Predictor, ncol=1, scale="free_y")
Fig5B
pdf("Figure5B.pdf", width = 11, height = 8.5)
plot(Fig5B)
dev.off()

#Figure 5C - Generate each heatmap independently and then merge in Adobe Illustrator 
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure5C-S9_03302021.csv", header = TRUE, na.strings = "NA")
df14 = subset(data, ModelNumber == "Model14")
df16 = subset(data, ModelNumber == "Model16")
df6 = subset(data, ModelNumber == "Model6")
df8 = subset(data, ModelNumber == "Model8")
df10 = subset(data, ModelNumber == "Model10")
df12 = subset(data, ModelNumber == "Model12")
df14$Variable = factor(df14$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df16$Variable = factor(df16$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df6$Variable = factor(df6$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df8$Variable = factor(df8$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df10$Variable = factor(df10$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df12$Variable = factor(df12$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
m14=ggplot(df14, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m14
pdf("Figure5C_EI.pdf", width = 11, height = 8.5)
plot(m14)
dev.off()

m16=ggplot(df16, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m16
pdf("Figure5C_Ratio.pdf", width = 11, height = 8.5)
plot(m16)
dev.off()

m6=ggplot(df6, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m6
pdf("Figure5C_Susceptibility1.pdf", width = 11, height = 8.5)
plot(m6)
dev.off()

m8=ggplot(df8, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m8
pdf("Figure5C_Susceptibility2.pdf", width = 11, height = 8.5)
plot(m8)
dev.off()

m10=ggplot(df10, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m10
pdf("Figure5C_Severity1.pdf", width = 11, height = 8.5)
plot(m10)
dev.off()

m12=ggplot(df12, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m12
pdf("Figure5C_Severity2.pdf", width = 11, height = 8.5)
plot(m12)
dev.off()


####Figure S3 ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("FigureS3_03302021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(1:10,15:28)]
cor(df_cor)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
pdf("FigureS3.pdf", width = 11, height = 8.5)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
dev.off()

####Figure S5 ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("FigureS5_03302021.csv", header = TRUE, na.strings = "NA")
data$susceptibility_avg_category2 = factor(data$susceptibility_avg_category2, levels=c("low", "average", "high"))
ggplot (data, aes(cat, NumWhorls))+
  geom_jitter(size = 1)
pdf("FigureS5.pdf", width = 11, height = 8.5)
ggplot (data, aes(cat, NumWhorls))+
  geom_jitter(size = 1)
dev.off()

####Figure S6 ####
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv("FigureS6_03302021.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
lodging$Replicate = as.character(lodging$Replicate)
lodging$Plot = as.character(lodging$Plot)
A1=levelplot(Percent.Lodging ~ Column*Row, data=lodging, main="Percent lodging", col.regions = heat.colors(100))
A2=levelplot(Lodging.Index ~ Column*Row, data=lodging, main="Lodging Index", col.regions = heat.colors(100))
A3=levelplot(LR_1.LS_2.or.E_3 ~ Column*Row, data=lodging, main="Genetic vs. Variable", col.regions = heat.colors(100))
plot = plot_grid(A1, A2, A3, labels=c('A', 'B', 'C'), label_size = 12, ncol = 3)
pdf("FigureS6.pdf", width = 11, height = 8.5)
plot(plot)
dev.off()

####Figure S7 ####
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("FigureS7-S8_03302021.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
data$Accession = factor(data$Accession,levels = c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71",
                                                  "A632","Mo17","PHK46","PHHB4","Ki3","CML322","CML52","SA24","CML10","Tzi9","Tzi8",
                                                  "Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103",
                                                  "Tx303","B97","R84","R78","PHB47","Ky21","CML258","B73","CML341","Shoe Peg","LH123Ht",
                                                  "2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
data$PLOT_ID2 = as.factor(data$PLOT_ID2)
S7A=ggplot(data=data,aes(x=Accession, y=w1, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of roots on W1 (count)")
S7A
pdf("FigureS7A.pdf", width = 11, height = 8.5)
plot(S7A)
dev.off()

S7B=ggplot(data=data,aes(x=Accession, y=w2, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of roots on W2 (count)")
S7B
pdf("FigureS7B.pdf", width = 11, height = 8.5)
plot(S7B)
dev.off()

S7C=ggplot(data=data,aes(x=Accession, y=w3, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of roots on W3 (count)")
S7C
pdf("FigureS7C.pdf", width = 11, height = 8.5)
plot(S7C)
dev.off()

S7D=ggplot(data=data,aes(x=Accession, y=single_root_width.cm, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Single root width (cm)")
S7D
pdf("FigureS7D.pdf", width = 11, height = 8.5)
plot(S7D)
dev.off()

S7E=ggplot(data=data,aes(x=Accession, y=stalk_width.cm, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Stalk width (cm)")
S7E
pdf("FigureS7E.pdf", width = 11, height = 8.5)
plot(S7E)
dev.off()

S7F=ggplot(data=data,aes(x=Accession, y=root_heightonstalk.cm, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Root height on stalk (cm)")
S7F
pdf("FigureS7F.pdf", width = 11, height = 8.5)
plot(S7F)
dev.off()

S7G=ggplot(data=data,aes(x=Accession, y=stalk_to_rootgrounding.cm, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Stalk to root grounding (cm)")
S7G
pdf("FigureS7G.pdf", width = 11, height = 8.5)
plot(S7G)
dev.off()

S7H=ggplot(data=data,aes(x=Accession, y=root_angle, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Root angle (degrees)")
S7H
pdf("FigureS7H.pdf", width = 11, height = 8.5)
plot(S7H)
dev.off()

S7I=ggplot(data=data,aes(x=Accession, y=spread_width.cm, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Spread width(cm)")
S7I
pdf("FigureS7I.pdf", width = 11, height = 8.5)
plot(S7I)
dev.off()

S7J=ggplot(data=data,aes(x=Accession, y=num_whorls, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of whorls in the soil - GUI (count)")
S7J
pdf("FigureS7J.pdf", width = 11, height = 8.5)
plot(S7J)
dev.off()

S7K=ggplot(data=data,aes(x=Accession, y=Brace.Root.Whorls.in.the.Soil, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of whorls in the soil - Manual (count)")
S7K
pdf("FigureS7K.pdf", width = 11, height = 8.5)
plot(S7K)
dev.off()

S7L=ggplot(data=data,aes(x=Accession, y=height.cm., color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Height (cm)")
S7L
pdf("FigureS7L.pdf", width = 11, height = 8.5)
plot(S7L)
dev.off()


####Figure S8####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("FigureS7-S8_03302021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(5,7:17)]
cor(df_cor)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
pdf("Corr_2019_AllRepsPheno.pdf", width = 11, height = 8.5)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
dev.off()

####Figure S9####
#Generate each heatmap independently and then merge in Adobe Illustrator; Models without height were generated in Figure 5C script 
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure5C-S9_03302021.csv", header = TRUE, na.strings = "NA")
df5 = subset(data, ModelNumber == "Model5")
df7 = subset(data, ModelNumber == "Model7")
df9 = subset(data, ModelNumber == "Model9")
df11 = subset(data, ModelNumber == "Model11")
df13 = subset(data, ModelNumber == "Model13")
df15 = subset(data, ModelNumber == "Model15")
df5$Variable = factor(df5$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df7$Variable = factor(df7$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df9$Variable = factor(df9$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df11$Variable = factor(df11$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df13$Variable = factor(df13$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
df15$Variable = factor(df15$Variable,levels = c("w3","num_whorls","brace_root_whorls_in_the_soil","w2","spread_width_cm","stalk_to_rootgrounding_cm","root_angle","root_heightonstalk_cm","w1","single_root_width_cm","stalk_width_cm","height_cm","ratio_none_all_imu","flexural_rigidity_ei_imu"))
m13=ggplot(df13, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m13
pdf("FigureS9_EI.pdf", width = 11, height = 8.5)
plot(m13)
dev.off()

m15=ggplot(df15, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m15
pdf("FigureS9_Ratio.pdf", width = 11, height = 8.5)
plot(m15)
dev.off()

m5=ggplot(df5, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m5
pdf("FigureS9_Susceptibility1.pdf", width = 11, height = 8.5)
plot(m5)
dev.off()

m7=ggplot(df7, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m7
pdf("FigureS9_Susceptibility2.pdf", width = 11, height = 8.5)
plot(m7)
dev.off()

m9=ggplot(df9, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m9
pdf("FigureS9_Severity1.pdf", width = 11, height = 8.5)
plot(m9)
dev.off()

m11=ggplot(df11, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
m11
pdf("FigureS9_Severity2.pdf", width = 11, height = 8.5)
plot(m11)
dev.off()

###Table 1 ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("Table1_03292021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(2:19)]
df=cor(df_cor)
df = as.data.frame(df)
write.csv(df, "Table1_03292021.csv", quote = F, row.names = T)