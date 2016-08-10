#Full code to run all analyses in Granados et al. (2016) Size and variation in individual growth rates among food web modules

####packages required 
library(plyr)
library(reshape)
library(RCurl)
library(ggplot2)
library(gridExtra)
library(car)


####load data from GitHub
RFU.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/RFU.csv")
RFU<-read.csv(text=RFU.URL)

SC.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/StandardCurves.csv")
SC<-read.csv(text=SC.URL)

Artemia.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/Artemia.csv")
Artemia<-read.csv(text=Artemia.URL)

Ratios.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/RNADNARatios.csv")
Ratios<-read.csv(text=Ratios.URL)

#C-control
#R-omnivory
#Z-exploitative competition
#N-food-chain
#M-consumer-resource

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
PhytoSixTrial<-ddply(.data=PhytoSix, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density))
PhytoSixMean<-ddply(.data=PhytoSixTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
PhytoSixMean$Configuration<-c("Control", "Consumer-resource", "Food chain", "Omnivory","Exploitative competition")

#order factors
PhytoSixMean$Configuration<- factor(PhytoSixMean$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Food chain", "Omnivory"))

####plot Figure S1a 
FigureS1a<-ggplot(PhytoSixMean, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (cells/ml)")+
  xlab("Treatments")+scale_y_continuous(limits=c(0,160000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 6-hour phytoplankton data 
PhytoAnova6<-aov(mean1 ~ ID, PhytoSixTrial)
summary(PhytoAnova6)

#Tukey HSD to assess differences between configurations 
PhytoTukey6<-TukeyHSD(PhytoAnova6)
PhytoTukey6Results<-PhytoTukey6$ID

#24-hour experiments 
Mesocosm<-strtrim(PhytoTwentyFour$ID, 1)
PhytoTwentyFour$ID<-Mesocosm

#calculate means
PhytoTwentyFourTrial<-ddply(.data=PhytoTwentyFour, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density))
PhytoTwentyFourMean<-ddply(.data=PhytoTwentyFourTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
PhytoTwentyFourMean$Configuration<-c("Control", "Consumer-resource", "Omnivory","Exploitative competition")

#order factors
PhytoTwentyFourMean$Configuration<- factor(PhytoTwentyFourMean$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition","Omnivory"))

####plot Figure S1b 
FigureS1b<-ggplot(PhytoTwentyFourMean, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (cells/ml)")+
  xlab("Treatments")+scale_y_continuous(limits=c(0,170000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 24-hour phytoplankton data 
PhytoAnova24<-aov(mean1 ~ ID, PhytoTwentyFourTrial)
summary(PhytoAnova24)

#Tukey HSD to assess differences between configurations 
PhytoTukey24<-TukeyHSD(PhytoAnova24)
PhytoTukey24Results<-PhytoTukey24$ID

#Stich graphs together 
FigureS1<-grid.arrange(arrangeGrob(FigureS1a, FigureS1b, ncol=1, nrow=2))


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
ArtemiaSixTrial<-ddply(.data=ArtemiaSix, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density), se=sd(Density)/sqrt(length(Density)))
ArtemiaSixSum<-ddply(.data=ArtemiaSixTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
ArtemiaSixSum$Configuration<-c("Control", "Consumer-resource", "Food chain", "Omnivory", "Exploitative competition")

#order factors
ArtemiaSixSum$Configuration<- factor(ArtemiaSixSum$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Food chain", "Omnivory"))

####plot Figure S2a 
FigureS2a<-ggplot(ArtemiaSixSum, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (ind/L)")+
  xlab("Treatments")+scale_y_continuous(limits=c(0,300))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#one-way ANOVA for 6-hour Artemia data
ArtemiaANOVA6<-aov(mean1 ~ ID, ArtemiaSixTrial)
summary(ArtemiaANOVA6)

#Tukey HSD to assess differences between configurations 
ArtemiaTukey6<-TukeyHSD(ArtemiaANOVA6, which="ID")
ArtemiaTukey6Results<-ArtemiaTukey6$ID

#24-hour experiments
ArtemiaTwentyFour<-Artemia[Artemia$Trial%in%c(1,2,3,7),]

Mesocosm<-strtrim(ArtemiaTwentyFour$ID, 1)
ArtemiaTwentyFour$ID<-Mesocosm

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
ArtemiaTwentyFour<-ArtemiaTwentyFour[ ! ArtemiaTwentyFour$ID %in% c("N"), ]

#calculate mean
ArtemiaTwentyFourTrial<-ddply(.data=ArtemiaTwentyFour, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Density))
ArtemiaTwentyFourSum<-ddply(.data=ArtemiaTwentyFourTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
ArtemiaTwentyFourSum$Configuration<-c("Control", "Consumer-resource", "Omnivory", "Exploitative competition")

#order factors
ArtemiaTwentyFourSum$Configuration<- factor(ArtemiaTwentyFourSum$Configuration, levels=c("Control", "Consumer-resource", "Exploitative competition", "Omnivory"))

####plot Figure S2b 
FigureS2b<-ggplot(ArtemiaTwentyFourSum, aes(x = Configuration, y = mean2))+
  geom_point()+geom_errorbar(aes(ymin=mean2-se, ymax=mean2+se, width=0.2))+
  theme_bw()+
  ylab("Density (ind/L)")+
  xlab("Treatments")+scale_y_continuous(limits=c(0,300))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#one-way ANOVA for 24-hour Artemia data
ArtemiaANOVA24<-aov(mean1 ~ ID, ArtemiaTwentyFourTrial)
summary(ArtemiaANOVA24)

#Tukey HSD to assess differences between configurations
ArtemiaTukey24<-TukeyHSD(ArtemiaANOVA24, which="ID")
ArtemiaTukey24Results<-ArtemiaTukey24$ID

#Stich graphs together 
FigureS2<-grid.arrange(arrangeGrob(FigureS2a, FigureS2b, ncol=1, nrow=2))

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

#calculate means 
RatioSixTrial<-ddply(.data=RatioSix, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Ratio))
RatioSixSum<-ddply(.data=RatioSixTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
RatioSixSum$Configuration<-c("Control", "Consumer-Resource", "Food chain", "Omnivory", "Exploitative Competition")

#order factors
RatioSixSum$Configuration<- factor(RatioSixSum$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Food chain", "Omnivory"))

####plot Figure 3a 
Figure3a<-ggplot(RatioSixSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("RNA:DNA ratio")+
  xlab("Treatments")+
  scale_y_continuous(limits=c(0,10))

#one-way ANOVA for 6-hour RNA:DNA ratios
RatioANOVA6<-aov(mean1 ~ ID, RatioSixTrial)
summary(RatioANOVA6)

#Tukey HSD to assess differences between configurations
RatioTukey6<-TukeyHSD(RatioANOVA6, which="ID")
RatioTukey6Results<-RatioTukey6$ID

#24-hour experiments 
RatioTwentyFour<-Ratios[Ratios$Trial%in%c(1,2,3,7),]

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
RatioTwentyFour<-RatioTwentyFour[ ! RatioTwentyFour$ID %in% c("N"), ]

Mesocosm<-strtrim(RatioTwentyFour$ID, 1)

RatioTwentyFour$ID<-Mesocosm

#calculate mean
RatioTwentyFourTrial<-ddply(.data=RatioTwentyFour, .variables=.(ID, Trial), .fun= summarise, mean1 = mean(Ratio))
RatioTwentyFourSum<-ddply(.data=RatioTwentyFourTrial, .variables=.(ID), .fun= summarise, mean2 = mean(mean1), se=sd(mean1)/sqrt(length(mean1)))

#convert codes to food web configurations 
RatioTwentyFourSum$Configuration<-c("Control", "Consumer-Resource", "Omnivory","Exploitative Competition")

#order factors
RatioTwentyFourSum$Configuration<- factor(RatioTwentyFourSum$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Omnivory"))

####plot Figure 3b 
Figure3b<-ggplot(RatioTwentyFourSum, aes(x = Configuration, y = mean))+
  geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))+
  theme_bw()+
  ylab("RNA:DNA ratio")+
  xlab("Treatments")+
  scale_y_continuous(limits=c(0,10))

#one-way ANOVA for 24-hour RNA:DNA ratios
RatioANOVA24<-aov(mean1 ~ ID, RatioTwentyFourTrial)
summary(RatioANOVA24)

#Tukey HSD to assess differences between configurations
RatioTukey24<-TukeyHSD(RatioANOVA24, which="ID")
RatioTukey24Results<-RatioTukey24$ID

#Stich graphs together 
Figure3<-grid.arrange(arrangeGrob(Figure3a, Figure3b, ncol=1, nrow=2))

#calculate variance 
Ratio6Var<-ddply(.data=RatioSix, .variables=.(ID, Trial), .fun= summarise, var1 = var(Ratio))

#convert codes to food web configurations 
Ratio6Var$Configuration<-rep(c("Control", "Consumer-Resource", "Food chain", "Omnivory", "Exploitative Competition"), each=3)

#remove NA
Ratio6Var<-na.omit(Ratio6Var)

#order factors
Ratio6Var$Configuration<- factor(Ratio6Var$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Food chain", "Omnivory"))

#calculate variance
Ratio24Var<-ddply(.data=RatioTwentyFour, .variables=.(ID, Trial), .fun= summarise, var = var(Ratio))

#convert codes to food web configurations 
Ratio24Var$Configuration<-NA
Ratio24Var[1:8,4]<-rep(c("Control", "Consumer-Resource"), each=4)
Ratio24Var[9:10,4]<-rep(c("Omnivory"), each=2)
Ratio24Var[11:14,4]<-rep(c("Exploitative Competition"), each=4)

#order factors
Ratio24Var$Configuration<- factor(Ratio24Var$Configuration, levels=c("Control", "Consumer-Resource", "Exploitative Competition", "Omnivory"))

#remove NA
Ratio24Var<-na.omit(Ratio24Var)

#Bartlett’s test of homogeneity of variance
Ratio6Var<-Ratio6Var[ ! Ratio6Var$ID %in% c("R"), ] #Remove single omnivory entry 
bartlett.test(var1 ~ Configuration, data=Ratio6Var)

bartlett.test(var ~ Configuration, data=Ratio24Var)

####plot Figure 4a
Figure4a<-ggplot(Ratio6Var, aes(x = Configuration, y = var1))+
  geom_boxplot()+ylab("Variance")+xlab("Treatments")+theme_bw()

Figure4b<-ggplot(Ratio24Var, aes(x = Configuration, y = var))+
  geom_boxplot()+ylab("Variance")+xlab("Treatments")+theme_bw()

#Stich graphs together 
Figure4<-grid.arrange(arrangeGrob(Figure4a, Figure4b, ncol=1, nrow=2))






