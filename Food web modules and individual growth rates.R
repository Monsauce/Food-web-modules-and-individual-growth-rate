#Full code to run all analyses in Granados et al. (2015) RNA:DNA ratios reveal higher consumer growth rates in food webs with omnivory 

####packages required 
library(plyr)
library(reshape)
library(RCurl)


####load data from GitHub
RFU.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/RFU.csv")
RFU<-read.csv(text=RFU.URL)

SC.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Food-web-modules-and-individual-growth-rate/master/StandardCurves.csv")
SC<-read.csv(text=SC.URL)


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
Trial1<-Trial1[ ! Trial1$ID %in% c("R1"), ]
Trial2<-Trial2[ ! Trial2$ID %in% c("R2","Z2","M1","M2"), ]
Trial4<-Trial4[ ! Trial4$ID %in% c("Z2"), ]
Trial5<-Trial5[ ! Trial5$ID %in% c("R1"), ]
Trial6<-Trial6[ ! Trial6$ID %in% c("R1", "R2"), ]

#make new data frames for the 6 and 24-hour experiments 
AlgaeSix<-rbind(Trial4, Trial5, Trial6)
AlgaeTwentyFour<-rbind(Trial1, Trial2, Trial3, Trial7)

#remove food-chain from the 24-hour experiment since we have no data for predation and RNA:DNA ratios
AlgaeTwentyFour<-AlgaeTwentyFour[ ! AlgaeTwentyFour$ID %in% c("N"), ]


