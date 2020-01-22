setwd("C:/Users/Seyni DIOP/Documents/ETUDES/ITS3/Projet/")
library(readr)
library(urca)
library(stats)
library(ggplot2)
library(ggthemes)
library(tseries)
library(knitr)
library(vars)
library(stats)
library(grDevices)
library(forecast)
library(timeSeries)
library(foreign)
library(gridExtra)
library(tsDyn)

IPI <- read_delim("IPI.csv", ";", escape_double = FALSE, 
                  col_types = cols(Date = col_datetime(format = "%d/%m/%Y")), 
                  trim_ws = TRUE)


ICA <- read_delim("ICA.csv", ";", escape_double = FALSE, 
                  col_types = cols(Date = col_datetime(format = "%d/%m/%Y")), 
                  trim_ws = TRUE)


IPC <- read_delim("IPC.csv", ";", escape_double = FALSE, 
                  col_types = cols(Date = col_datetime(format = "%d/%m/%Y")), 
                  trim_ws = TRUE)


CE <- read_delim("CE.csv", ";", escape_double = FALSE, 
                  col_types = cols(Date = col_datetime(format = "%d/%m/%Y")), 
                  trim_ws = TRUE)



IPI=ts(IPI$IPI,start=c(2008,1),frequency = 12)
ICA=ts(ICA$ICA,start=c(2008,1),frequency = 12)
IPC=ts(IPC$IPC,start=c(2008,1),frequency = 12)
CE=ts(CE$CE,start=c(2008,1),frequency = 12)


IPI=window(IPI,start=c(2008,1),end=c(2016,9))
ICA=window(ICA,start=c(2008,1),end=c(2016,9))
IPC=window(IPC,start=c(2008,1),end=c(2016,9))
CE =window(CE,start=c(2008,1),end=c(2016,9))

data <- read_delim("data.csv", ";", escape_double = FALSE, 
                   col_types = cols(Date = col_datetime(format = "%d/%m/%Y ")), 
                   trim_ws = TRUE)
data.ts=ts(data[,2:5],start = c(2008,1),frequency = 12,names = c("IPI","IPC","ICA","CE") )



# 3.1.1 Descriptive
summary(IPI)
summary(ICA)
summary(IPC)
summary(CE)

#********************* Représentation graphique *********************#

p1=autoplot(IPI, xlab="Année", ylab="IPI")+theme_igray()
p2=autoplot(ICA, xlab="Année", ylab="ICA")+theme_igray()
p3=autoplot(IPC,xlab="Année", ylab="IPC")+theme_igray()
p4=autoplot(CE ,xlab="Année", ylab="CE")+theme_igray()
grid.arrange(p1,p2,p3,p4)


#--------------------------------------------------------------------------
#************************TEST DE STATIONNARITE********************#
IPI_log_cent<-log(IPI)-mean(log(IPI))
ICA_log_cent<-log(ICA)-mean(log(ICA))
IPC_log_cent<-log(IPC)-mean(log(IPC))
CE_log_cent<-log(CE)-mean(log(CE))
#********** Variable : Indice de la Production Industrielle

autoplot(IPI_log_cent,xlab="Année", ylab="Indice de la Production Industrielle")
kpss.test(IPI_log_cent)

#KPSS Test for Level Stationarity
#
#data:  IPI_log_cent
#KPSS Level = 1.0223, Truncation lag parameter = 2, p-value = 0.01
#le test de kpss donne une p-value de 0.01437 donc on rejette la stationnarité

#*********** Variable : Indice des chiffres d'affaires

autoplot(ICA_log_cent, xlab="Année", ylab="Indice des chiffres d'affaires")+theme_cowplot()
kpss.test(ICA_log_cent)

#data:  ICA_log_cent
#KPSS Level = 0.69098, Truncation lag parameter = 2, p-value = 0.01437
#le test de kpss donne une p-value de 0.01437 donc on rejette la stationnarité


#*********** Variable :Indice des Prix à la consommation

autoplot(IPC_log_cent,xlab="Année", ylab="Indice des Prix à la consommation")
kpss.test(IPC_log_cent)

#data:  IPC_log_cent
#KPSS Level = 2.5675, Truncation lag parameter = 2, p-value = 0.01
#le test de kpss donne une p-value de 0.01437 donc on rejette la stationnarité

#*********** Variable :Crédit à l'économie****************

autoplot(CE_log_cent,xlab="Année", ylab="Crédits à l'économie")
kpss.test(CE_log_cent)

#data:  CE_log_cent
#KPSS Level = 3.588, Truncation lag parameter = 2, p-value = 0.01


#****************** Stationnarisation des données *****************#

#********* Pour l'IPI **********#
IPI_log_cent_tdf<-diff(IPI_log_cent,lag=1)
ts.plot(IPI_log_cent_tdf,xlab="Année", ylab="Indice des chiffres d'affaires")
kpss.test(IPI_log_cent_tdf)

#KPSS Test for Level Stationarity

#data:  IPI_log_cent_tdf
#KPSS Level = 0.025554, Truncation lag parameter = 2, p-value = 0.1


#********* Pour l'ICA **********#

ICA_log_cent_tdf<-diff(ICA_log_cent,lag=1)
ts.plot(ICA_log_cent_tdf,xlab="Année", ylab="Indice des chiffres d'affaires")
kpss.test(ICA_log_cent_tdf)

#data:  ICA_log_cent_tdf
#KPSS Level = 0.041926, Truncation lag parameter = 2, p-value = 0.1
#On constate qu'aprés différenciation la série devient stationnaire



#********* Pour l'IPC **********#

IPC_log_cent_tdf<-diff(IPC_log_cent,lag=1)
ts.plot(IPC_log_cent_tdf,xlab="Année", ylab="Indice des Prix à la consommation")
kpss.test(IPC_log_cent_tdf)

#data:  IPC_log_cent_tdf
#KPSS Level = 0.041457, Truncation lag parameter = 2, p-value = 0.1
#On constate qu'aprés différenciation la série devient stationnaire

#********* Pour le CE **********#

CE_log_cent_tdf<-diff(CE_log_cent,lag=1)
ts.plot(CE_log_cent_tdf,xlab="Année", ylab="Indice du Chiffre d'Affaire")
kpss.test(CE_log_cent_tdf)

#data:  CE_log_cent_tdf
#KPSS Level = 0.080667, Truncation lag parameter = 2, p-value = 0.1




#******

new_data<-matrix(ncol=4,nrow=105)
new_data[1:105,1]<-IPI
new_data[1:105,2]<-IPC
new_data[1:105,3]<-ICA
new_data[1:105,4]<-CE/10000
new_data<-ts(new_data,start=c(2008,1),frequency = 12,names=c("IPI","IPC","ICA","CE"))

#***********Test de cointégration de Johansen*************#
#Réestimation de l'ordre p
VARselect(new_data)

#***********Nouvelles Estimation des paramètres
#On choisit p=2
new_coint<-ca.jo(new_data,type="trace",K=2)
summary(new_coint)
NEW_ESTI=VECM(new_data,lag=1,r=1,estim = "ML",include = "const")
summary(NEW_ESTI)
toLatex(summary(NEW_ESTI),parenthese = "Pvalue")


#************VALIDATION****************#

#ACF des résidus
k1<-ggAcf(NEW_ESTI$residuals[,1],lag.max=10)+labs(title="ACF_RESIDU1")+theme_igray()
k2<-ggAcf(NEW_ESTI$residuals[,2],lag.max=10)+labs(title="ACF_RESIDU2")+theme_igray()
k3<-ggAcf(NEW_ESTI$residuals[,3],lag.max=10)+labs(title="ACF_RESIDU3")+theme_igray()
k4<-ggAcf(NEW_ESTI$residuals[,4],lag.max=10)+labs(title="ACF_RESIDU4")+theme_igray()

grid.arrange(k1,k2,k3,k4,nrow=2,ncol=2)

#Test d'autocorrélations des résidus 
library(portes)

portest(NEW_ESTI$residuals,ncores = 8)

normality.test(vec2var(new_coint,r=1))

###########  CAUSALITE  ############

#/********POUR IPC
IPC_cause_IPI<-matrix(ncol=2,nrow=105)
IPC_cause_IPI[1:105,1]<-IPC
IPC_cause_IPI[1:105,2]<-IPI
IPC_cause_IPI<-ts(IPC_cause_IPI,start=c(2008,1),frequency = 12,names=c("IPC","IPI"))
causality(VAR(y=IPC_cause_IPI,p=2),cause="IPC")

#/********POUR ICA
ICA_cause_IPI<-matrix(ncol=2,nrow=105)
ICA_cause_IPI[1:105,1]<-ICA
ICA_cause_IPI[1:105,2]<-IPI
ICA_cause_IPI<-ts(ICA_cause_IPI,start=c(2008,1),frequency = 12,names=c("ICA","IPI"))
causality(VAR(y=ICA_cause_IPI,p=2),cause="ICA")

#/********POUR CE
CE_cause_IPI<-matrix(ncol=2,nrow=105)
CE_cause_IPI[1:105,1]<-CE
CE_cause_IPI[1:105,2]<-IPI
CE_cause_IPI<-ts(CE_cause_IPI,start=c(2008,1),frequency = 12,names=c("CE","IPI"))
causality(VAR(y=CE_cause_IPI,p=2),cause="CE")

###PREVISION DE L'IPI AVEC LE MODELE MCE

IPI_pred<-predict(NEW_ESTI,n.ahead = 12)
IPI_pred


IPI_reconst=matrix(ncol=2,nrow=117)
IPI_reconst[1:105,2]=IPI
IPI_reconst[1:105,1]=IPI
IPI_reconst[106:117,1]=IPI_pred[,1]
IPI_reconst<-ts(IPI_reconst,start = c(2008,1),frequency = 12,names=c("Prévision","IPI"))
autoplot(IPI_reconst)+labs(title="IPI",ylab="IPI")+theme_igray()

#Fonction d'impulsion
par(mfrow = c(2,2))
plot(irf(NEW_ESTI))

#Décomposition de la variance


k=12
VAR_DECOMP<-fevd(NEW_ESTI,n.ahead = k)$IPI

DECOMPOSITION<-matrix(nrow = 4*k,ncol=3)

colnames(DECOMPOSITION)<-c("Index","Valeurs","Variables")

DECOMPOSITION[1:12,1]<-as.numeric(seq(1,12,1))
DECOMPOSITION[1:12,2]<-VAR_DECOMP[,1]
DECOMPOSITION[1:12,3]<-"IPI"

DECOMPOSITION[13:24,1]<-as.numeric(seq(1,12,1))
DECOMPOSITION[13:24,2]<-VAR_DECOMP[,2]
DECOMPOSITION[13:24,3]<-"IPC"

DECOMPOSITION[25:36,1]<-as.numeric(seq(1,12,1))
DECOMPOSITION[25:36,2]<-VAR_DECOMP[,3]
DECOMPOSITION[25:36,3]<-"ICA"

DECOMPOSITION[37:48,1]<-as.numeric(seq(1,12,1))
DECOMPOSITION[37:48,2]<-VAR_DECOMP[,4]
DECOMPOSITION[37:48,3]<-"CE"

DECOMPOSITION<-as.data.frame(DECOMPOSITION)
DECOMPOSITION$Index<-as.numeric(as.character(DECOMPOSITION$Index))
DECOMPOSITION$Valeurs<-as.numeric(as.character(DECOMPOSITION$Valeurs))

p1 <- ggplot(DECOMPOSITION, aes(x=Index, y=Valeurs,color=Variables))+
  geom_line(size=1.5) +ggtitle("DECOMPOSITION DE LA VARIANCE DE PREVISION DE L'IPI")+xlab("TIME")+ylab("Part de variance")
p1

