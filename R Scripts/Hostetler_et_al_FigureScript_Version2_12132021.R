#Updated Figure script for Hostetler_et_al_2021 PCE-Version 2-2021/12
#Scripts are seperated into associate Figure for the manuscript
setwd(dir = "/Users/ashleyhostetler/Desktop/Hostetler_et_al_2021/R Input Files/Figure Files/")
####Libary Load
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
#Figure 1####
##Determining Genotype Order - High to low averages of BRC ####
#Ratio data - Figure 1
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure1_12142021.csv", header = TRUE, na.strings = "NA")
x = aggregate(data$Ratio..None.All....IMU, list(data$Accession), FUN=mean) 
x
#write.table(x, file="RatioData_AccessionOrder.csv", sep=",", row.names = FALSE)
data$Accession = factor(data$Accession,levels = c("CML52","SA24","Ky21","Oh7B","Ki11","CML10","CML69","P39","A632","NC350",
                                                  "CML228","NC262","Tzi8","PHT60","B104","Oh43","CML277","F118","CML341","B97",
                                                  "CML373","Ki3","CML333","Tzi9","2369","Mo18W","CML322","CML247","M37W","Ms71",
                                                  "W64A","CML103","R84","PHB47","M162W","Tx303","R78","6M502","PHK46","Hickory King",
                                                  "Shoe Peg","L163","LH123Ht","PHR58","B73","LH211","Mo17","CML258","LH252","PHHB4",
                                                  "Hp301","GT112"))
##Figure Generation####
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure1_12142021.csv", header = TRUE, na.strings = "NA")
data$Year = as.factor(data$Year)
data$PLOT_ID2 = as.factor(data$PLOT_ID2)
data$PLOT_ID3 = as.factor(data$PLOT_ID3)
data$Accession = factor(data$Accession,levels = c("CML52","SA24","Ky21","Oh7B","Ki11","CML10","CML69","P39","A632","NC350",
                                                  "CML228","NC262","Tzi8","PHT60","B104","Oh43","CML277","F118","CML341","B97",
                                                  "CML373","Ki3","CML333","Tzi9","2369","Mo18W","CML322","CML247","M37W","Ms71",
                                                  "W64A","CML103","R84","PHB47","M162W","Tx303","R78","6M502","PHK46","Hickory King",
                                                  "Shoe Peg","L163","LH123Ht","PHR58","B73","LH211","Mo17","CML258","LH252","PHHB4",
                                                  "Hp301","GT112"))
Fig1=ggplot(data=data,aes(x=Accession, y=Ratio..None.All....IMU, color=Year))+
  geom_boxplot(size = 0.25, color = "black",fill="gray80") + 
  geom_point(size = .75) +
  theme_minimal()+
  scale_colour_manual(values=c("red4","darkslategray4"))+
  xlab("Genotype") + 
  ylab("Ratio=None/All")+
  theme(axis.text.x = element_text(size=11,angle=90), 
        axis.text.y = element_text(size=11), 
        plot.title=element_text(size=11, vjust=3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_text(vjust=2.5), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"),
        legend.title=element_blank())
Fig1
pdf("Figure1.pdf", width = 8.5, height = 4.5)
plot(Fig1)
dev.off()
#Figure 2####
##Figure 2A####
#Generate each heatmap independently and then merge in Adobe Illustrator
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure2A-4_12152021.csv", header = TRUE, na.strings = "NA")

df15 = subset(data, ModelNumber == "Model15")
df3 = subset(data, ModelNumber == "Model3")

df15$Variable = factor(df15$Variable,levels = c("w3","w2","root_angle","spread_width_cm","stalk_to_rootgrounding_cm","w1","root_heightonstalk_cm","stalk_width_cm","height_cm","brace_root_whorls_in_the_soil","single_root_width_cm"))
df3$Variable = factor(df3$Variable,levels = c("w3","w2","root_angle","spread_width_cm","stalk_to_rootgrounding_cm","w1","root_heightonstalk_cm","stalk_width_cm","height_cm","brace_root_whorls_in_the_soil","single_root_width_cm"))

m3=ggplot(df3, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=11, vjust=3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_blank(), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"),
        legend.title=element_blank())
m3

m15=ggplot(df15, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=11, vjust=3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_blank(), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"),
        legend.title=element_blank())
m15
Figure2A=plot_grid(m3, m15, labels = c(), nrow = 1)
Figure2A
pdf("Figure2A.pdf", width = 4, height = 2.75)
plot(Figure2A)
dev.off()

##Figure 2B####
cat("\014")
rm(list=ls()) 
ls() 
df <- read.csv("Figure2B_12142021.csv", header = TRUE, na.strings = "NA")
df <- df[,c(1,3:11)]
pca <- prcomp(df, scale. = TRUE)
var <- get_pca_var(pca)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
pdf("Figure3B.pdf", width = 5, height = 4)
ggplot2::autoplot(pca, data = df, loadings = TRUE, loadings.label= TRUE)
dev.off()

#Figure 3####
##Figure 3A####
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv("Figure3A_12142021.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
lodging$Lodged2 = as.numeric(lodging$Lodged2)
Fig3A=ggplot(lodging, aes(x=Lodged2)) +
  geom_histogram(bins=2,color="black")+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11), 
        plot.title=element_text(size=11, vjust=-3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_text(vjust=2.5), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.5,"cm"),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"))+
  scale_y_continuous(name="Number of plants", breaks=seq(0,900,100))+
  scale_x_continuous(name="Lodging susceptibility", breaks=seq(0,2,1))
Fig3A
pdf("Figure3A.pdf", width = 2.25, height = 2.25)
plot(Fig3A)
dev.off()

##Figure 3B-C####
cat("\014")
rm(list=ls()) 
ls() 
lodging <- read.csv("Figure3B-C_12142021.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
accessions <- read.csv("../Unique_55Genotypes_Plots.csv", header = TRUE, strip.white=TRUE,na.strings = "NA")
colnames(lodging)
colnames(accessions)
accessions$Genotype = accessions$Accession
accessions$Accession = NULL
lodging$Replicate = as.character(lodging$Replicate)
lodging$Plot = as.character(lodging$Plot)
summary = summary(comparedf(lodging, accessions, by="Genotype"))
notShared = summary$obs.table
lodging = lodging[!lodging$Genotype %in% notShared$Genotype,]
head(lodging)
lodging$Accession = lodging$Genotype

data_3B = lodging %>%
  dplyr::filter(Accession %in% c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71","A632","Mo17","PHK46","PHHB4","Ki3","CML322"))
data_3B %>% count(Accession)

data_3C = lodging %>%
  dplyr::filter(Accession %in% c("CML52","SA24","CML10","Tzi9","Tzi8","Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103","Tx303","B97","R84","R78","PHB47","Ky21","CML258","B73"))
data_3C %>% count(Accession)

data_3C2 = lodging %>%
  dplyr::filter(Accession %in% c("CML341","Shoe Peg","LH123Ht","2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
data_3C2 %>% count(Accession)

#Figure 3B
colnames(data_3B)
lodging3 = data_3B[,c(3,13,8,9,12)]
lodging6 = lodging3[,c(1,3,4,5)]
lodging3_long = melt(lodging6, id = c("Plot","Lodging.Index.Cat"))
lodging3_long$Lodging.Index.Cat = as.character(lodging3_long$Lodging.Index.Cat)
lodging3_long$Lodging.Index.Cat[lodging3_long$variable == "Total.Plants"] = "Total"
lodging3$Accession = factor(lodging3$Accession,levels = c("B104","CML228","6M502","CML277","Hp301","PHR58","CML373","L163","Ki11","Ms71","A632","Mo17","PHK46","PHHB4","Ki3","CML322"))
plotB=ggplot(data=lodging3_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("red4","red1","gold2","skyblue4","gray60","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=9, vjust=3), 
        axis.text=element_text(size=9), 
        axis.title = element_text(size=9), 
        axis.title.y= element_blank(), 
        axis.title.x= element_blank(), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 9, colour = "black"),
        strip.text.y = element_text(size = 9, colour = "black"),
        legend.title=element_blank())+
  facet_grid(. ~ lodging3$Accession, scales="free", space="free")
plotB
pdf("Figure3B.pdf", width = 5.25, height = 2.25)
plot(plotB)
dev.off()

#Figure 3C
colnames(data_3C)
lodging4 = data_3C[,c(3,13,8,9,12)]
lodging7 = lodging4[,c(1,3,4,5)]
lodging4_long = melt(lodging7, id = c("Plot","Lodging.Index.Cat"))
lodging4_long$Lodging.Index.Cat = as.character(lodging4_long$Lodging.Index.Cat)
lodging4_long$Lodging.Index.Cat[lodging4_long$variable == "Total.Plants"] = "Total"
lodging4$Accession = factor(lodging4$Accession,levels = c("CML52","SA24","CML10","Tzi9","Tzi8","Mo18W","CML69","Oh7B","Hickory King","NC350","P39","CML247","CML333","LH252","CML103","Tx303","B97","R84","R78","PHB47","Ky21","CML258","B73"))
plotC=ggplot(data=lodging4_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("skyblue4","skyblue4","skyblue4","skyblue4","skyblue4","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("skyblue4","skyblue4","skyblue4","skyblue4","skyblue4","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=9, vjust=3), 
        axis.text=element_text(size=9), 
        axis.title = element_text(size=9), 
        axis.title.y= element_blank(), 
        axis.title.x= element_blank(), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 9, colour = "black"),
        strip.text.y = element_text(size = 9, colour = "black"),
        legend.title=element_blank())+
  facet_grid(. ~ lodging4$Accession, scales="free", space="free")
plotC
pdf("Figure3C.pdf", width = 8.5, height = 2.25)
plot(plotC)
dev.off()

#Figure 3C2
colnames(data_3C2)
lodging2 = data_3C2[,c(3,13,8,9,12)]
colnames(lodging2)
lodging5 = lodging2[,c(1,3,4,5)]
lodging2_long = melt(lodging5, id = c("Plot","Lodging.Index.Cat"))
lodging2_long$Lodging.Index.Cat = as.character(lodging2_long$Lodging.Index.Cat)
lodging2_long$Lodging.Index.Cat[lodging2_long$variable == "Total.Plants"] = "Total"
lodging2$Accession = factor(lodging2$Accession,levels = c("CML341","Shoe Peg","LH123Ht","2369","NC262","M37W","GT112","F118","PHT60","M162W","LH211","Oh43","W64A"))
plotC2=ggplot(data=lodging2_long,aes(x=Plot, y=value, fill=Lodging.Index.Cat, color=Lodging.Index.Cat, alpha=variable)) +
  geom_bar(stat="identity",position ="identity") +
  scale_colour_manual(values=c("skyblue4","skyblue4","skyblue4","skyblue4","skyblue4","gray60"),
                      breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_fill_manual(values=c("skyblue4","skyblue4","skyblue4","skyblue4","skyblue4","gray60"),
                    breaks=c("Severe","Moderate","SlightModerate","Slight","NoEffect","Total")) +
  scale_alpha_manual(values=c(.8,.5)) +
  geom_text(aes(label=value), color="black", vjust=1.6)+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=9, vjust=3), 
        axis.text=element_text(size=9), 
        axis.title = element_text(size=9), 
        axis.title.y= element_blank(), 
        axis.title.x= element_blank(), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 9, colour = "black"),
        strip.text.y = element_text(size = 9, colour = "black"),
        legend.title=element_blank())+
  facet_grid(. ~ lodging2$Accession, scales="free", space="free")
plotC2
pdf("Figure3C2.pdf", width = 6.5, height = 2.25)
plot(plotC2)
dev.off()
##Figure3D####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("Figure3D_12142021.csv", header = TRUE, na.strings = "NA")
head(data)
data1 = data %>%
  group_by(Accession) %>%
  summarise(NumWhorls = mean(NumWhorls))
data2 = data %>%
  group_by(Accession) %>%
  summarise(susceptibility_avg = mean(susceptibility_avg))
data3 = merge(data1, data2, by=c("Accession"))
plot3D=ggplot (data3, aes(x=susceptibility_avg, y=NumWhorls))+
  geom_point()+
  geom_smooth(method='lm')+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        plot.title=element_text(size=10, vjust=3), 
        axis.text=element_text(size=12), 
        axis.title = element_text(size=10), 
        axis.title.y= element_text(vjust=2.5), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        legend.title=element_blank())
plot3D
summary(lm(data3$NumWhorls ~ data3$susceptibility_avg))
lm(data3$NumWhorls ~ data3$susceptibility_avg)

pdf("Figure3D.pdf", width = 3, height = 2.25)
plot(plot3D)
dev.off()

#Figure 4####
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("Figure2A-4_12152021.csv", header = TRUE, na.strings = "NA")
df5 = subset(data, ModelNumber == "Model5")
df7 = subset(data, ModelNumber == "Model7")

df5$Variable = factor(df5$Variable,levels = c("w3","brace_root_whorls_in_the_soil","spread_width_cm","w1","w2","stalk_to_rootgrounding_cm","root_heightonstalk_cm","single_root_width_cm","root_angle","stalk_width_cm","height_cm","ratio_none_all_imu"))
df7$Variable = factor(df7$Variable,levels = c("w3","brace_root_whorls_in_the_soil","spread_width_cm","w1","w2","stalk_to_rootgrounding_cm","root_heightonstalk_cm","single_root_width_cm","root_angle","stalk_width_cm","height_cm","ratio_none_all_imu"))

m5=ggplot(df5, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=11, vjust=3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_blank(), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"),
        legend.title=element_blank())
m5

m7=ggplot(df7, aes(ModelNumber, Variable, fill=Importance)) + 
  geom_tile()+
  geom_text(aes(label = round(Importance, 3)))+
  scale_fill_continuous(low = "azure", high = "dodgerblue3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=11, vjust=3), 
        axis.text=element_text(size=11), 
        axis.title = element_text(size=11), 
        axis.title.y= element_blank(), 
        axis.title.x= element_text(vjust=-1.4), 
        axis.ticks.length = unit(.2,"cm"),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="grey"),
        strip.text.x = element_text(size = 11, colour = "blue"),
        strip.text.y = element_text(size = 11, colour = "black"),
        legend.title=element_blank())
m7

Figure4=plot_grid(m5, m7, labels = c(), nrow = 1)
Figure4
pdf("Figure4.pdf", width = 4, height = 2.75)
plot(Figure4)
dev.off()

#Figure S3 ####
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

#Figure S5####
cat("\014")
rm(list=ls()) 
ls() 
data = read.csv("FigureS5-S6_12152021.csv", header = TRUE, strip.white = TRUE, na.strings = "NA")
data$Accession = factor(data$Accession,levels = c("CML52","SA24","Ky21","Oh7B","Ki11","CML10","CML69","P39","A632","NC350",
                                                  "CML228","NC262","Tzi8","PHT60","B104","Oh43","CML277","F118","CML341","B97",
                                                  "CML373","Ki3","CML333","Tzi9","2369","Mo18W","CML322","CML247","M37W","Ms71",
                                                  "W64A","CML103","R84","PHB47","M162W","Tx303","R78","6M502","PHK46","Hickory King",
                                                  "Shoe Peg","L163","LH123Ht","PHR58","B73","LH211","Mo17","CML258","LH252","PHHB4",
                                                  "Hp301","GT112"))
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

S7J=ggplot(data=data,aes(x=Accession, y=Brace.Root.Whorls.in.the.Soil, color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Number of whorls in the soil (count)")
S7J
pdf("FigureS7J.pdf", width = 11, height = 8.5)
plot(S7J)
dev.off()

S7K=ggplot(data=data,aes(x=Accession, y=height.cm., color=PLOT_ID2))+
  geom_boxplot(size = 0.25, color = "black", fill="gray60") + 
  geom_point(size = 1.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  scale_colour_manual(values=c("skyblue4","gold2","red4","red1"))+
  xlab("Genotype") + 
  ylab("Height (cm)")
S7K
pdf("FigureS7K.pdf", width = 11, height = 8.5)
plot(S7K)
dev.off()

#Figure S6####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("FigureS5-S6_12152021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(5,8:17)]
cor(df_cor)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
pdf("FigureS6.pdf", width = 11, height = 8.5)
corrplot(cor(df_cor), method = "number", type = "upper", order = "FPC", tl.cex = 0.6,  number.cex = 0.6, diag = FALSE, tl.pos = "td")
dev.off()


#Table 1 ####
cat("\014")
rm(list=ls()) 
ls() 
data <- read.csv("Table1_10212021.csv", header = TRUE, na.strings = "NA")
colnames(data)
df_cor = data[,c(2:4,6:16)]
df=cor(df_cor)
df = as.data.frame(df)
write.csv(df, "Table1_Results.csv", quote = F, row.names = T)
