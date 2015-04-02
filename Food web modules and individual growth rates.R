#Full code to run all analyses in Granados et al. (2015) RNA:DNA ratios reveal higher consumer growth rates in food webs with omnivory 

####packages required 
library(plyr)
library(reshape)
library(RCurl)
library(ggplot2)


####load data from GitHub
RFU.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/RFU.csv")
RFU<-read.csv(text=RFU.URL)

SC.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/StandardCurves.csv")
SC<-read.csv(text=SC.URL)

Artemia.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/Artemia.csv")
Artemia<-read.csv(text=Artemia.URL)

Ratios.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/RNADNARatios.csv")
Ratios<-read.csv(text=Ratios.URL)


#R, N, Z, M and C are omnivory, food-chain, exploitative competition and control respectively 


####phytoplankton analyses 

#read in standard curves
SC1<-subset(SC, Trial==1)
SC2<-subset(SC, Trial==2)
SC3<-subset(SC, Trial==3)
SC4<-subset(SC, Trial==4)
SC5<-subset(SC, Trial==5)
SC6<-subset(SC, Trial==6)
SC7<-subset(SC, Trial==7)

#get equation for standard curve 
Fit1<-lm(SC1$RFU~SC1$Density)
Fit2<-lm(SC2$RFU~SC2$Density)
Fit3<-lm(SC3$RFU~SC3$Density)
Fit4<-lm(SC4$RFU~SC4$Density)
Fit5<-lm(SC5$RFU~SC5$Density)
Fit6<-lm(SC6$RFU~SC6$Density)
Fit7<-lm(SC7$RFU~SC7$Density)

Trial1<-ddply(.data=RFU[RFU$Trial == 1,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit1$coefficients[1])/Fit1$coefficients[2])
Trial2<-ddply(.data=RFU[RFU$Trial == 2,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit2$coefficients[1])/Fit2$coefficients[2])
Trial3<-ddply(.data=RFU[RFU$Trial == 3,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit3$coefficients[1])/Fit3$coefficients[2])
Trial4<-ddply(.data=RFU[RFU$Trial == 4,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit4$coefficients[1])/Fit4$coefficients[2])
Trial5<-ddply(.data=RFU[RFU$Trial == 5,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit5$coefficients[1])/Fit5$coefficients[2])
Trial6<-ddply(.data=RFU[RFU$Trial == 6,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit6$coefficients[1])/Fit6$coefficients[2])
Trial7<-ddply(.data=RFU[RFU$Trial == 7,], .variables=.(ID, Trial, RFU), transform, Density = (RFU-Fit7$coefficients[1])/Fit7$coefficients[2])

#remove mussels which stopped feeding across entire set of experiments 
Trial1<-Trial1[ ! Trial1$ID %in% c("R1"), ] #Remove Trial 1 R1
Trial2<-Trial2[ ! Trial2$ID %in% c("R2","Z2","M1","M2"), ] #Remove Trial 2 R2, Z2, M1, M2
Trial4<-Trial4[ ! Trial4$ID %in% c("Z2"), ]  #Remove Trial 4 Z2
Trial5<-Trial5[ ! Trial5$ID %in% c("R1"), ]  #Remove Trial 5 R1
Trial6<-Trial6[ ! Trial6$ID %in% c("R1", "R2"), ]  #Remove Trial 6 R1, R2 

#make new data frames for the 6 and 24-hour experiments 
PhytoSix<-rbind(Trial4, Trial5, Trial6)
PhytoTwentyFour<-rbind(Trial1, Trial2, Trial3, Trial7)

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
PhytoTwentyFour<-PhytoTwentyFour[ ! PhytoTwentyFour$ID %in% c("N"), ]

#6-hour experiments 
Mesocosm<-strtrim(PhytoSix$ID, 1)
PhytoSix$ID<-Mesocosm

#calculate means
PhytoSixMean1<-ddply(.data=PhytoSix, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density))

PhytoSixMean2<-ddply(.data=PhytoSixMean1, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
PhytoSixMean2$Configuration<-c("Control", "Consumer-resource", "Food chain", "Omnivory","Exploitative competition")

#order factors
PhytoSixMean2$Configuration<- factor(PhytoSixMean2$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Food chain", "Omnivory"))

####plot Figure S1a 
FigureS1a<-ggplot(PhytoSixMean2, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (cells/ml)")+
  xlab("Food web configuration")+scale_y_continuous(limits=c(0,160000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 6-hour phytoplankton data 
PhytoAnova6<-aov(mean1 ~ ID, Phyto6Mean1)
summary(PhytoAnova6)

#Tukey HSD to assess differences between configurations 
PhytoTukey6<-TukeyHSD(PhytoAnova6)
PhytoTukey6Results<-PhytoTukey6$ID

#24-hour experiments 
Mesocosm<-strtrim(PhytoTwentyFour$ID, 1)
PhytoTwentyFour$ID<-Mesocosm

#calculate means
PhytoTwentyFourMean1<-ddply(.data=PhytoTwentyFour, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density))

PhytoTwentyFourMean2<-ddply(.data=PhytoTwentyFourMean1, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
PhytoTwentyFourMean2$Configuration<-c("Control", "Consumer-resource", "Omnivory","Exploitative competition")

#order factors
PhytoTwentyFourMean2$Configuration<- factor(PhytoTwentyFourMean2$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition","Omnivory"))

####plot Figure S1b 
FigureS1b<-ggplot(AlgaeTwentyFourMean2, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (cells/ml)")+
  xlab("Configuration")+scale_y_continuous(limits=c(0,170000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 24-hour phytoplankton data 
PhytoAnova24<-aov(mean1 ~ ID, PhytoTwentyFourMean1)
TukeyHSD(PhytoAnova24)
summary(PhytoAnova24)
PhytoTukey24<-TukeyHSD(PhytoAnova24)
PhytoTukey24Results<-PhytoTukey24$ID

####Artemia analyses

#remove mussels which stopped feeding across entire set of experiments 

Artemia<-Artemia[!Artemia$Trial == 1 | !Artemia$ID == "R1",] #Remove Trial 1 R1
Artemia<-Artemia[!Artemia$Trial == 2 | !Artemia$ID %in% c("R2", "Z2", "M1", "M2"),] 
##Remove Trial 2 R2, Z2, M1, M2
Artemia<-Artemia[!Artemia$Trial == 4 | !Artemia$ID == "Z2",] #Remove Trial 4 Z2
Artemia<-Artemia[!Artemia$Trial == 5 | !Artemia$ID == "R1",] #Remove Trial 5 R1
Artemia<-Artemia[!Artemia$Trial == 6 | !Artemia$ID %in% c("R1", "R2"),]  #Remove Trial 6 R1, R2 

#6-hour experiments 
ArtemiaSix<-Artemia[Artemia$Trial%in%c(4,5,6),]

Mesocosm<-strtrim(ArtemiaSix$ID, 1)
ArtemiaSix$ID<-Mesocosm

#calculate mean 
ArtemiaSixSum<-ddply(.data=ArtemiaSix, .variables=.(ID), .fun= summarise, mean = mean(Density), se=sd(Density)/sqrt(length(Density)))

#convert codes to food web configurations 
ArtemiaSixSum$Configuration<-c("Control", "Consumer-resource", "Food chain", "Omnivory", "Exploitative competition")

#order factors
ArtemiaSixSum$Configuration<- factor(ArtemiaSixSum$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Food chain", "Omnivory"))

####plot Figure S2a 
FigureS2a<-ggplot(ArtemiaSixSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("Density (ind/L)")+
  xlab("Configuration")+scale_y_continuous(limits=c(0,225))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#one-way ANOVA for 6-hour Artemia data
ArtemiaANOVA6<-aov(Density ~ ID, ArtemiaSix)
ArtemiaTukey6<-TukeyHSD(ArtemiaANOVA6, which="ID")
ArtemiaTukey6Results<-ArtemiaTukey6$ID

#24-hour experiments
ArtemiaTwentyFour<-Artemia[Artemia$Trial%in%c(1,2,3,7),]

Mesocosm<-strtrim(ArtemiaTwentyFour$ID, 1)
ArtemiaTwentyFour$ID<-Mesocosm

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
ArtemiaTwentyFour<-ArtemiaTwentyFour[ ! ArtemiaTwentyFour$ID %in% c("N"), ]

#calculate mean
ArtemiaTwentyFourSum<-ddply(.data=ArtemiaTwentyFour, .variables=.(ID), .fun= summarise, mean = mean(Density), se=sd(Density)/sqrt(length(Density)))

#convert codes to food web configurations 
ArtemiaTwentyFourSum$Configuration<-c("Control", "Consumer-resource", "Omnivory", "Exploitative competition")

#order factors
ArtemiaTwentyFourSum$Configuration<- factor(ArtemiaTwentyFourSum$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Omnivory"))

####plot Figure S2b 
FigureS2b<-ggplot(ArtemiaTwentyFourSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("Density (ind/L)")+
  xlab("Configuration")+scale_y_continuous(limits=c(0,225))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 24-hour Artemia data
ArtemiaANOVA24<-aov(Density ~ ID, ArtemiaTwentyFour)
ArtemiaTukey24<-TukeyHSD(ArtemiaANOVA24, which="ID")
ArtemiaTukey24Results<-ArtemiaTukey24$ID

####RNA:DNA Analyses 

#remove mussels which stopped feeding across entire set of experiments 
Ratios<-Ratios[!Ratios$Trial == 1 | !Ratios$ID == "R1",] #Remove Trial 1 R1
Ratios<-Ratios[!Ratios$Trial == 2 | !Ratios$ID %in% c("R2", "Z2", "M1", "M2"),] 
##Remove Trial 2 R2, Z2, M1, M2
Ratios<-Ratios[!Ratios$Trial == 4 | !Ratios$ID == "Z2",] #Remove Trial 4 Z2
Ratios<-Ratios[!Ratios$Trial == 5 | !Ratios$ID == "R1",] #Remove Trial 5 R1
Ratios<-Ratios[!Ratios$Trial == 6 | !Ratios$ID %in% c("R1", "R2"),]  #Remove Trial 6 R1, R2

#6-hour experiments 
RatioSix<-Ratios[Ratios$Trial%in%c(4,5,6),]

Mesocosm<-strtrim(RatioSix$ID, 1)

RatioSix$ID<-Mesocosm

#calculate mean
RatioSixSum<-ddply(.data=RatioSix, .variables=.(ID), .fun= summarise, mean = mean(Ratio), se=sd(Ratio)/sqrt(length(Ratio)))

#convert codes to food web configurations 
RatioSixSum$Configuration<-c("Control", "Consumer-Resource", "Food chain", "Omnivory", "Exploitative Competition")

#order factors
RatioSixSum$Configuration<- factor(RatioSixSum$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Food chain", "Omnivory"))

####plot Figure 3a 
Figure3a<-ggplot(RatioSixSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("RNA:DNA ratio")+
  xlab("Configuration")+
  scale_y_continuous(limits=c(0,8))

#one-way ANOVA for 6-hour RNA:DNA ratios
RatioANOVA6<-aov(Ratio ~ ID, RatioSix)
TukeyHSD(RatioANOVA6)
RatioTukey6<-TukeyHSD(RatioANOVA6, which="ID")
RatioTukey6Results<-RatioTukey6$ID


#24-hour experiments 
RatioTwentyFour<-Ratios[Ratios$Trial%in%c(1,2,3,7),]

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
RatioTwentyFour<-RatioTwentyFour[ ! RatioTwentyFour$ID %in% c("N"), ]

Mesocosm<-strtrim(RatioTwentyFour$ID, 1)

RatioTwentyFour$ID<-Mesocosm

#calculate mean
RatioTwentyFourSum<-ddply(.data=RatioTwentyFour, .variables=.(ID), .fun= summarise, mean = mean(Ratio), se=sd(Ratio)/sqrt(length(Ratio)))

#convert codes to food web configurations 
RatioTwentyFourSum$Configuration<-c("Control", "Consumer-Resource", "Omnivory","Exploitative Competition")

#order factors
RatioTwentyFourSum$Configuration<- factor(RatioTwentyFourSum$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Omnivory"))

####plot Figure 3b 
Figure3b<-ggplot(RatioTwentyFourSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("RNA:DNA ratio")+
  xlab("Configuration")+
  scale_y_continuous(limits=c(0,8))

#one-way ANOVA for 24-hour RNA:DNA ratios
RatioANOVA24<-aov(Ratio ~ ID, RatioTwentyFour)
TukeyHSD(RatioANOVA24)
RatioTukey24<-TukeyHSD(RatioANOVA24, which="ID")
RatioTukey24Results<-RatioTukey24$ID







