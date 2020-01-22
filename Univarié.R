setwd("C:/Users/Seyni DIOP/Documents/ETUDES/ITS3/Projet/")
## ----CHARGEMENT DES PACKAGES----------------------------------------------
library(readr)
library(ggplot2)
library(ggfortify)
library(forecast) 
library(TSA)
library(Kendall)
library(cowplot)
library(gridExtra)
library(ggthemes)
library(trend)
library(astsa)
library(scales)
library(lmtest)
library(normtest)
library(urca)

## ----CHARGEMENT DE LA BASE-------------------------------------------------
BASE <- read_delim("IPI.csv", 
                       ";", escape_double = FALSE, 
                   col_types = cols(Date = col_datetime(format = "%d/%m/%Y")),
                       trim_ws = TRUE)
BASE.ini=ts(BASE$IPI,start=c(2008,1),frequency=12)
## ----CHRONOGRAMME---------
IPI1<-ggplot(data = BASE, aes(x = Date, y = IPI))+
  geom_line(color = "blue", size = 0.5)+ stat_smooth(
    color = "#FC4E07", fill = "#FC4E07",
    method = "loess",se=F
  )+theme_gray()


IPI1

## ----DECOMPOSITION -------------------------------------------
decomposition<-autoplot(stl(BASE.ini,s.window="periodic"))

decomposition+theme_gray()

## ----Construction ACF1 et PACF1

ACF1<-ggAcf(BASE.ini, lag.max = 50,
            plot = TRUE, na.action = na.contiguous, demean = TRUE)+ 
  labs(title="Fonction d'autocorrélation",
       x ="Retard", y = "Autocorrélation")+theme_igray()
  
PACF1<-ggPacf(BASE.ini, lag.max = 50,
              plot = TRUE, na.action = na.contiguous, demean = TRUE)+ 
  labs(title="Fonction d'Autocorrélation Partielle",
       x ="Retard", y = "Autocorrélation Partielle")+
  theme_igray()

grid.arrange(ACF1,PACF1,ncol=1,nrow=2)


csmk.test(BASE.ini,alternative=c("greater"))
#Le test de Mann-Kendall dit qu'il n'y a pas de tendance de à 5%
## ---- PERIODOGRAMME -------------------------
p.star = periodogram (BASE.ini, log = "yes", plot = FALSE)
periodo<-autoplot(p.star,color="red")+geom_vline(xintercept=1/12,show.legend = T,color="red")+
  labs(title="Périodogramme", x ="Fréquence", y = "Spectre")+theme_igray()

r=which.max(p.star$spec)
T=120
p=1/p.star$freq[r] # ou encore p=T/r

periodo

plot(log(p.star$spec), type="l",lwd=2, xlab="Indice (r) de fréquence de Fourier",
     ylab="Spectre discret ", col=rgb(0.5, 0.7, 0.2),
     main = "Périodogramme des observations de l'IPI par indice de fréquence")
abline(v=10)

## ----CONSTRUCTION D'UN ECHANTILLON TEST -------------------------
data0<-window(BASE.ini,end=c(2015,12))
sdata<-diff(data0,lag = 12)
dsdata<-diff(sdata,lag = 1)
data<-dsdata

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
#### Test de racines unité et de stationnarité
adf.test(dsdata)
PP.test(dsdata)
kpss.test(dsdata)

## ----include=FALSE-------------------------------------------------------
##IPI DESAISONNALISEE
IPI2<-autoplot(data, ts.colour = "blue")+ labs(title="Série stationnarisée",x ="Date", y = "IPI*")+
  theme(plot.title = element_text(size=20, face="bold.italic")) +theme_gray()   #représentons la série sans saisonnalité

IPI2

## ----Construction de ACF2 et PACF2 ----
ACF2<-ggAcf(dsdata,lag.max = 50,
       plot = TRUE, na.action = na.contiguous, demean = TRUE)+ 
  labs(title="Fonction d'Autocorrélation",
       x ="Retard", y = "Autocorrélation")+theme_igray()
  

PACF2<-ggPacf(dsdata, lag.max = 50,
       plot = TRUE, na.action = na.contiguous, demean = TRUE)+ 
  labs(title="Fonction d'Autocorrélation Partielle",
       x ="Retard", y = "Autocorrélation partielle")+theme_igray()

grid.arrange(ACF2,PACF2,nrow=2,ncol=1)

## ----Déclaration modèle candidats, message=FALSE, warning=FALSE, include=FALSE----
#Déclaration des modèles candidats
MODELE1<-astsa::sarima(data0,p=0,d=1,q=1,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE2<-astsa::sarima(data0,p=0,d=1,q=1,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE3<-astsa::sarima(data0,p=0,d=1,q=1,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE4<-astsa::sarima(data0,p=0,d=1,q=1,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE5<-astsa::sarima(data0,p=0,d=1,q=1,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE6<-astsa::sarima(data0,p=0,d=1,q=1,P=2,D=1,Q=1,S=12 , Model = TRUE)
#***

MODELE7<-astsa::sarima(data0,p=1,d=1,q=0,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE8<-astsa::sarima(data0,p=1,d=1,q=0,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE9<-astsa::sarima(data0,p=1,d=1,q=0,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE10<-astsa::sarima(data0,p=1,d=1,q=0,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE11<-astsa::sarima(data0,p=1,d=1,q=0,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE12<-astsa::sarima(data0,p=1,d=1,q=0,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***
MODELE13<-astsa::sarima(data0,p=1,d=1,q=1,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE14<-astsa::sarima(data0,p=1,d=1,q=1,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE15<-astsa::sarima(data0,p=1,d=1,q=1,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE16<-astsa::sarima(data0,p=1,d=1,q=1,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE17<-astsa::sarima(data0,p=1,d=1,q=1,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE18<-astsa::sarima(data0,p=1,d=1,q=1,P=2,D=1,Q=1,S=12 , Model = TRUE)
  
#***

MODELE19<-astsa::sarima(data0,p=2,d=1,q=0,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE20<-astsa::sarima(data0,p=2,d=1,q=0,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE21<-astsa::sarima(data0,p=2,d=1,q=0,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE22<-astsa::sarima(data0,p=2,d=1,q=0,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE23<-astsa::sarima(data0,p=2,d=1,q=0,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE24<-astsa::sarima(data0,p=2,d=1,q=0,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***

MODELE25<-astsa::sarima(data0,p=2,d=1,q=1,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE26<-astsa::sarima(data0,p=2,d=1,q=1,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE27<-astsa::sarima(data0,p=2,d=1,q=1,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE28<-astsa::sarima(data0,p=2,d=1,q=1,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE29<-astsa::sarima(data0,p=2,d=1,q=1,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE30<-astsa::sarima(data0,p=2,d=1,q=1,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***

MODELE31<-astsa::sarima(data0,p=3,d=1,q=0,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE32<-astsa::sarima(data0,p=3,d=1,q=0,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE33<-astsa::sarima(data0,p=3,d=1,q=0,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE34<-astsa::sarima(data0,p=3,d=1,q=0,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE35<-astsa::sarima(data0,p=3,d=1,q=0,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE36<-astsa::sarima(data0,p=3,d=1,q=0,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***

MODELE37<-astsa::sarima(data0,p=3,d=1,q=1,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE38<-astsa::sarima(data0,p=3,d=1,q=1,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE39<-astsa::sarima(data0,p=3,d=1,q=1,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE40<-astsa::sarima(data0,p=3,d=1,q=1,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE41<-astsa::sarima(data0,p=3,d=1,q=1,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE42<-astsa::sarima(data0,p=3,d=1,q=1,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***

MODELE43<-astsa::sarima(data0,p=0,d=1,q=0,P=0,D=1,Q=0,S=12 , Model = TRUE)
MODELE44<-astsa::sarima(data0,p=0,d=1,q=0,P=0,D=1,Q=1,S=12 , Model = TRUE)

MODELE45<-astsa::sarima(data0,p=0,d=1,q=0,P=1,D=1,Q=0,S=12 , Model = TRUE)
MODELE46<-astsa::sarima(data0,p=0,d=1,q=0,P=1,D=1,Q=1,S=12 , Model = TRUE)

MODELE47<-astsa::sarima(data0,p=0,d=1,q=0,P=2,D=1,Q=0,S=12 , Model = TRUE)
MODELE48<-astsa::sarima(data0,p=0,d=1,q=0,P=2,D=1,Q=1,S=12 , Model = TRUE)

#***
## ----Vérification de convergence, message=FALSE, warning=FALSE, include=FALSE----
#Vérification de la convergence des modèles
MODELE1$fit$code
MODELE2$fit$code
MODELE3$fit$code
MODELE4$fit$code
MODELE5$fit$code
MODELE6$fit$code
MODELE7$fit$code
MODELE8$fit$code
MODELE9$fit$code
MODELE10$fit$code
MODELE11$fit$code
MODELE12$fit$code
MODELE13$fit$code
MODELE14$fit$code
MODELE15$fit$code
MODELE16$fit$code
MODELE17$fit$code
MODELE18$fit$code
MODELE19$fit$code
MODELE20$fit$code
MODELE21$fit$code
MODELE22$fit$code
MODELE23$fit$code
MODELE24$fit$code
MODELE25$fit$code
MODELE26$fit$code
MODELE27$fit$code
MODELE28$fit$code
MODELE29$fit$code
MODELE30$fit$code
MODELE31$fit$code
MODELE32$fit$code
MODELE33$fit$code
MODELE34$fit$code
MODELE35$fit$code
MODELE36$fit$code
MODELE37$fit$code
MODELE38$fit$code
MODELE39$fit$code
MODELE40$fit$code
MODELE41$fit$code
MODELE42$fit$code
MODELE43$fit$code
MODELE44$fit$code
MODELE45$fit$code
MODELE46$fit$code
MODELE47$fit$code
MODELE48$fit$code
## ----Significativité des paramètres, message=FALSE, warning=FALSE, include=FALSE----
#Significativité des paramètres
abs(MODELE1$ttable[,1])>1.96*MODELE1$ttable[,2]
abs(MODELE2$ttable[,1])>1.96*MODELE2$ttable[,2]
abs(MODELE3$ttable[,1])>1.96*MODELE3$ttable[,2]
abs(MODELE4$ttable[,1])>1.96*MODELE4$ttable[,2]
abs(MODELE5$ttable[,1])>1.96*MODELE5$ttable[,2]
abs(MODELE6$ttable[,1])>1.96*MODELE6$ttable[,2]
abs(MODELE7$ttable[,1])>1.96*MODELE7$ttable[,2]
abs(MODELE8$ttable[,1])>1.96*MODELE8$ttable[,2]
abs(MODELE9$ttable[,1])>1.96*MODELE9$ttable[,2]
abs(MODELE10$ttable[,1])>1.96*MODELE10$ttable[,2]
abs(MODELE11$ttable[,1])>1.96*MODELE11$ttable[,2]
abs(MODELE12$ttable[,1])>1.96*MODELE12$ttable[,2]
abs(MODELE13$ttable[,1])>1.96*MODELE13$ttable[,2]
abs(MODELE14$ttable[,1])>1.96*MODELE14$ttable[,2]
abs(MODELE15$ttable[,1])>1.96*MODELE15$ttable[,2]
abs(MODELE16$ttable[,1])>1.96*MODELE16$ttable[,2]
abs(MODELE17$ttable[,1])>1.96*MODELE17$ttable[,2]
abs(MODELE18$ttable[,1])>1.96*MODELE18$ttable[,2]
abs(MODELE19$ttable[,1])>1.96*MODELE19$ttable[,2]
abs(MODELE20$ttable[,1])>1.96*MODELE20$ttable[,2]
abs(MODELE21$ttable[,1])>1.96*MODELE21$ttable[,2]
abs(MODELE22$ttable[,1])>1.96*MODELE22$ttable[,2]
abs(MODELE23$ttable[,1])>1.96*MODELE23$ttable[,2]
abs(MODELE24$ttable[,1])>1.96*MODELE24$ttable[,2]
abs(MODELE25$ttable[,1])>1.96*MODELE25$ttable[,2]
abs(MODELE26$ttable[,1])>1.96*MODELE26$ttable[,2]
abs(MODELE27$ttable[,1])>1.96*MODELE27$ttable[,2]
abs(MODELE28$ttable[,1])>1.96*MODELE28$ttable[,2]
abs(MODELE29$ttable[,1])>1.96*MODELE29$ttable[,2]
abs(MODELE30$ttable[,1])>1.96*MODELE30$ttable[,2]
abs(MODELE31$ttable[,1])>1.96*MODELE31$ttable[,2]
abs(MODELE32$ttable[,1])>1.96*MODELE32$ttable[,2]
abs(MODELE33$ttable[,1])>1.96*MODELE33$ttable[,2]
abs(MODELE34$ttable[,1])>1.96*MODELE34$ttable[,2]
abs(MODELE35$ttable[,1])>1.96*MODELE35$ttable[,2]
abs(MODELE36$ttable[,1])>1.96*MODELE36$ttable[,2]
abs(MODELE37$ttable[,1])>1.96*MODELE37$ttable[,2]
abs(MODELE38$ttable[,1])>1.96*MODELE38$ttable[,2]
abs(MODELE39$ttable[,1])>1.96*MODELE39$ttable[,2]
abs(MODELE40$ttable[,1])>1.96*MODELE40$ttable[,2]

abs(MODELE41$ttable[,1])>1.96*MODELE41$ttable[,2]
abs(MODELE42$ttable[,1])>1.96*MODELE42$ttable[,2]
abs(MODELE43$ttable[,1])>1.96*MODELE43$ttable[,2]
abs(MODELE44$ttable[,1])>1.96*MODELE44$ttable[,2]
abs(MODELE45$ttable[,1])>1.96*MODELE45$ttable[,2]
abs(MODELE46$ttable[,1])>1.96*MODELE46$ttable[,2]
abs(MODELE47$ttable[,1])>1.96*MODELE47$ttable[,2]
abs(MODELE48$ttable[,1])>1.96*MODELE48$ttable[,2]
# Les modèles valides sont 1 , 2, 3,5,7,8,9,11,19,20,21,23,31,32,33,35,44,45,46


## ----Affectation résidus, message=FALSE, warning=FALSE, include=FALSE----
#Affectation résidus
RESIDU.MODELE1<-MODELE1$fit$residuals
RESIDU.MODELE2<-MODELE2$fit$residuals
RESIDU.MODELE3<-MODELE3$fit$residuals
RESIDU.MODELE5<-MODELE5$fit$residuals

RESIDU.MODELE7<-MODELE7$fit$residuals
RESIDU.MODELE8<-MODELE8$fit$residuals

RESIDU.MODELE9<-MODELE9$fit$residuals
RESIDU.MODELE11<-MODELE11$fit$residuals
RESIDU.MODELE19<-MODELE19$fit$residuals


RESIDU.MODELE20<-MODELE20$fit$residuals
RESIDU.MODELE21<-MODELE21$fit$residuals
RESIDU.MODELE23<-MODELE23$fit$residuals
RESIDU.MODELE31<-MODELE31$fit$residuals
RESIDU.MODELE32<-MODELE32$fit$residuals

RESIDU.MODELE33<-MODELE33$fit$residuals
RESIDU.MODELE35<-MODELE35$fit$residuals

RESIDU.MODELE44<-MODELE44$fit$residuals
RESIDU.MODELE45<-MODELE45$fit$residuals
RESIDU.MODELE46<-MODELE46$fit$residuals


## ----Test de nullité de la moyenne des résidus, message=FALSE, warning=FALSE, include=FALSE----
#Test de nullité de la moyenne des résidus
abs(mean(RESIDU.MODELE1))<1.96*sqrt(var(RESIDU.MODELE1))
abs(mean(RESIDU.MODELE2))<1.96*sqrt(var(RESIDU.MODELE2))
abs(mean(RESIDU.MODELE3))<1.96*sqrt(var(RESIDU.MODELE3))
abs(mean(RESIDU.MODELE5))<1.96*sqrt(var(RESIDU.MODELE5))

abs(mean(RESIDU.MODELE7))<1.96*sqrt(var(RESIDU.MODELE7))
abs(mean(RESIDU.MODELE8))<1.96*sqrt(var(RESIDU.MODELE8))

abs(mean(RESIDU.MODELE9))<1.96*sqrt(var(RESIDU.MODELE9))
abs(mean(RESIDU.MODELE11))<1.96*sqrt(var(RESIDU.MODELE11))

abs(mean(RESIDU.MODELE19))<1.96*sqrt(var(RESIDU.MODELE19))
abs(mean(RESIDU.MODELE20))<1.96*sqrt(var(RESIDU.MODELE20))
abs(mean(RESIDU.MODELE21))<1.96*sqrt(var(RESIDU.MODELE21))
abs(mean(RESIDU.MODELE23))<1.96*sqrt(var(RESIDU.MODELE23))
abs(mean(RESIDU.MODELE31))<1.96*sqrt(var(RESIDU.MODELE31))
abs(mean(RESIDU.MODELE32))<1.96*sqrt(var(RESIDU.MODELE32))
abs(mean(RESIDU.MODELE33))<1.96*sqrt(var(RESIDU.MODELE33))
abs(mean(RESIDU.MODELE35))<1.96*sqrt(var(RESIDU.MODELE35))

abs(mean(RESIDU.MODELE44))<1.96*sqrt(var(RESIDU.MODELE44))
abs(mean(RESIDU.MODELE45))<1.96*sqrt(var(RESIDU.MODELE45))
abs(mean(RESIDU.MODELE46))<1.96*sqrt(var(RESIDU.MODELE46))

## ----Homoscédasticité des résidus, message=FALSE, warning=FALSE, include=FALSE----
#Test d'homoscédasticité des résidus

white.test(RESIDU.MODELE1)
white.test(RESIDU.MODELE2)
white.test(RESIDU.MODELE3)
white.test(RESIDU.MODELE5)
white.test(RESIDU.MODELE7)
white.test(RESIDU.MODELE8)
white.test(RESIDU.MODELE9)
white.test(RESIDU.MODELE11)
white.test(RESIDU.MODELE19)
white.test(RESIDU.MODELE20)
white.test(RESIDU.MODELE21)
white.test(RESIDU.MODELE23)
white.test(RESIDU.MODELE31)
white.test(RESIDU.MODELE32)
white.test(RESIDU.MODELE33)
white.test(RESIDU.MODELE35)

white.test(RESIDU.MODELE44)
white.test(RESIDU.MODELE45)
white.test(RESIDU.MODELE46)

## ----Test autocorrélation----
Box.test(RESIDU.MODELE1,lag = 1/3*length(RESIDU.MODELE1),type = c("Box-Pierce"),fitdf = 1 )

Box.test(RESIDU.MODELE2,lag = 1/3*length(RESIDU.MODELE2),type = c("Box-Pierce"),fitdf =2)

Box.test(RESIDU.MODELE3,lag = 1/3*length(RESIDU.MODELE3),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE5,lag = 1/3*length(RESIDU.MODELE5),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE7,lag = 1/3*length(RESIDU.MODELE7),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE8,lag = 1/3*length(RESIDU.MODELE8),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE9,lag = 1/3*length(RESIDU.MODELE9),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE11,lag = 1/3*length(RESIDU.MODELE11),type = c("Box-Pierce"),fitdf =1)

Box.test(RESIDU.MODELE19,lag = 1/3*length(RESIDU.MODELE19),type = c("Box-Pierce"),fitdf =2)

Box.test(RESIDU.MODELE20,lag = 1/3*length(RESIDU.MODELE20),type = c("Box-Pierce"),fitdf =2)

Box.test(RESIDU.MODELE21,lag = 1/3*length(RESIDU.MODELE21),type = c("Box-Pierce"),fitdf =2)

Box.test(RESIDU.MODELE23,lag = 1/3*length(RESIDU.MODELE23),type = c("Box-Pierce"),fitdf =2)

Box.test(RESIDU.MODELE31,lag = 1/3*length(RESIDU.MODELE31),type = c("Box-Pierce"),fitdf =3)

Box.test(RESIDU.MODELE32,lag = 1/3*length(RESIDU.MODELE32),type = c("Box-Pierce"),fitdf =3)

Box.test(RESIDU.MODELE33,lag = 1/3*length(RESIDU.MODELE33),type = c("Box-Pierce"),fitdf =3)

Box.test(RESIDU.MODELE35,lag = 1/3*length(RESIDU.MODELE35),type = c("Box-Pierce"),fitdf =3)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
AIC_MODELE1<-MODELE1$AIC
AICc_MODELE1<-MODELE1$AICc
BIC_MODELE1<-MODELE1$BIC
RMSE.MODELE1<-sqrt(mean(RESIDU.MODELE1^2))
MAPE.MODELE1<-1/length(RESIDU.MODELE1)*sum(abs(RESIDU.MODELE1/data0))
HQ.MODELE1<-log(var(RESIDU.MODELE1))+3*log(log(length(data0))/length(data0))


AIC_MODELE3<-MODELE3$AIC
AICc_MODELE3<-MODELE3$AICc
BIC_MODELE3<-MODELE3$BIC
RMSE.MODELE3<-sqrt(mean(RESIDU.MODELE3^2))
MAPE.MODELE3<-1/length(RESIDU.MODELE3)*sum(abs(RESIDU.MODELE3/data0))
HQ.MODELE3<-log(var(RESIDU.MODELE3))+3*log(log(length(data0))/length(data0))

AIC_MODELE5<-MODELE5$AIC
AICc_MODELE5<-MODELE5$AICc
BIC_MODELE5<-MODELE5$BIC
RMSE.MODELE5<-sqrt(mean(RESIDU.MODELE5^2))
MAPE.MODELE5<-1/length(RESIDU.MODELE5)*sum(abs(RESIDU.MODELE5/data0))
HQ.MODELE5<-log(var(RESIDU.MODELE5))+3*log(log(length(data0))/length(data0))


AIC_MODELE31<-MODELE31$AIC
AICc_MODELE31<-MODELE31$AICc
BIC_MODELE31<-MODELE31$BIC
RMSE.MODELE31<-sqrt(mean(RESIDU.MODELE31^2))
MAPE.MODELE31<-1/length(RESIDU.MODELE31)*sum(abs(RESIDU.MODELE31/data0))
HQ.MODELE31<-log(var(RESIDU.MODELE31))+3*(3+0)*log(log(length(data0))/length(data0))

AIC_MODELE35<-MODELE35$AIC
AICc_MODELE35<-MODELE35$AICc
BIC_MODELE35<-MODELE35$BIC
RMSE.MODELE35<-sqrt(mean(RESIDU.MODELE35^2))
MAPE.MODELE35<-1/length(RESIDU.MODELE35)*sum(abs(RESIDU.MODELE35/data0))
HQ.MODELE35<-log(var(RESIDU.MODELE35))+3*(3+0)*log(log(length(data0))/length(data0))

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
AIC_MODELE1
AIC_MODELE3
AIC_MODELE5
AIC_MODELE31
AIC_MODELE35

AICc_MODELE1
AICc_MODELE3
AICc_MODELE5
AICc_MODELE31
AICc_MODELE35

BIC_MODELE1
BIC_MODELE3
BIC_MODELE5
BIC_MODELE31
BIC_MODELE35

RMSE.MODELE1
RMSE.MODELE3
RMSE.MODELE5
RMSE.MODELE31
RMSE.MODELE35

MAPE.MODELE1
MAPE.MODELE3
MAPE.MODELE5
MAPE.MODELE31
MAPE.MODELE35

HQ.MODELE1
HQ.MODELE3
HQ.MODELE5
HQ.MODELE31
HQ.MODELE35

## ------------------------------------------------------------------------
jb.norm.test(RESIDU.MODELE5)
shapiro.test(RESIDU.MODELE5)

## ------------------------------------------------------------------------
###Graphique :histogramme
layout(matrix(2:1, ncol=2))

qqnorm(RESIDU.MODELE5,datax=TRUE,main="QQ-plot")
qqline(RESIDU.MODELE5,datax=TRUE) 

curve(dnorm(x,mean(RESIDU.MODELE5),sd(RESIDU.MODELE5)),col="green",add=T)
hist(RESIDU.MODELE5,prob=T,breaks=30,col="pink",main="Test graphique de la normaité des résidus",add=T)
lines(density(RESIDU.MODELE5),col="blue")

## Prévision
MOD5<-forecast::Arima(data0,order = c(0,1,1),seasonal = list(order=c(2,1,0),period=12))
MODELE5_forecast<-forecast(MOD5,h=12)

autoplot(MODELE5_forecast)+theme_gray()
 