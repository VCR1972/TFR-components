### This program decomposes the TFR into its fertility vs childlessness parts
### we use  life table 

library(tidyr)
library(patchwork)
library(dplyr)
options(scipen = 999)

LifeTableMx<-function(mx){ 
  
  N<-length(mx)
  
  ax<-rep(0.5,N) 
  
  qx<-mx/(1+(1-ax)*mx) 
  
  ### For mortality with the last value of qx = 1... here NOT
  ## qx[N] <- 1
  
  px<-1-qx
  
  lx<-100000
  
  for(y in 1:(N-1)){          
    lx[y+1]<-lx[y]*px[y]
  }
  
  dx<-lx*qx
  
  Lx<-lx[-1]+ax[-N]*dx[-N] 
  
  Lx[N]<-lx[N] ### this reduces to lx in the last age-group since nobody is having children at ages 55+                
  
  Tx<-c() 
  for(y in 1:N){
    Tx[y]<-sum(Lx[y:N]) 
  }
  
  Tx<-Tx-Tx[N] 
  ### since our interest is of life expectancy between ages 12 and 55, 
  ### we subtract the last value
  
  ex<-Tx/lx 
  
  Age<-12:55    
  
  ALL<-data.frame(Age,mx,lx,dx,Lx,Tx,ex)
  return(ALL)
}


Fertility1<-function(LTs,As){
  B<-As
  PC<-(LTs$lx/100000)
  AFR<-B$ASFR[1]
  ASFR1<-B$ASFR[1]/(1-PC[1])
  ASFR1[(is.nan(ASFR1))|(is.infinite(ASFR1))]<-0
  AFR1<-ASFR1
  for (x in 2:length(A1$Age)){
    AFR[x]<-sum(B$ASFR[1:x])
    ASFR1[x]<-0
    if (PC[x]<1){
      ASFR1[x]<-(AFR[x]/(1-PC[x]))-AFR1[x-1]}
    AFR1[x]<-sum(ASFR1[1:(x)])
  }
  ## to test
  ## AFR/(1-PC)
  return(ASFR1)
}


AgeDecomp<-function(lt1,b1,a1,lt2,b2,a2){ 
  ### checking for time 1
  
  PC<-lt1$lx/100000
  afr<-b1$AFR
  afr1<-b1$AFR1
  
  n<-length(afr)
  
  f1<-(afr[-1]/(1-PC[-1]))-afr1[-n]
  f1[(1-PC[-1])==0]<-0
  
  t1<-f1*(1-PC[-1])
  t2<-afr1[-n]*(PC[-n]-PC[-1])
  
  #  cbind(A1$ASFR[-1],t1+t2) 
  
  ### checking for time 2
  
  PC2<-lt2$lx/100000
  afr2<-b2$AFR
  afr12<-b2$AFR1
  
  n<-length(afr2)
  
  f12<-(afr2[-1]/(1-PC2[-1]))-afr12[-n]
  f12[(1-PC2[-1])==0]<-0
  
  t12<-f12*(1-PC2[-1])
  t22<-afr12[-n]*(PC2[-n]-PC2[-1])
  
  #   cbind(A2$ASFR[-1],t12+t22)   
  ### very good matching except for the first ages
  
  
  ## now the decomposition kitagawa type
  
  termPC<-(afr12[-n]+afr1[-n])/2*(PC2[-n]-PC2[-1]-(PC[-n]-PC[-1]))+
    (f12+f1)*((1-PC2[-1])-(1-PC[-1]))/2
  
  ## here we use instead differences to avoid issues with the negative f1
  termf1<-(f12-f1)*((1-PC[-1])+(1-PC2[-1]))/2 +
    (afr12[-n]-afr1[-n])*((PC2[-n]-PC2[-1])+(PC[-n]-PC[-1]))/2
  
  
  ASFR<-a1$ASFR[-1]
  ASFR2<-a2$ASFR[-1]
  
  CASFR<-ASFR2-ASFR
  D0<-cbind(ASFR,ASFR2,sqrt(ASFR2*ASFR),CASFR)
  
  Decomp<-termf1+termPC
  age<-Age[-1]
  TFR2<-sum(ASFR2)
  TFR1<-sum(ASFR)
  CTFR<-TFR2-TFR1
  c(colSums(D0),CTFR) 
  D<-cbind(age,termf1,termPC)
  
  Tot<-  colSums(D)
  
  D_d <- as.data.frame(D)
  
  return(D_d)
}


Name<-c("AUT","BLR","BEL","BGR","CAN","CHL",
        "HRV","CZE","DNK","EST","FIN",
        "DEUTNP","HUN","ISL","IRL",
        "ITA","JPN","LTU","NLD","NOR","POL","PRT",
        "KOR","RUS","SVK","SVN","ESP","SWE","CHE",
        "TWN","GBR_NP","GBRTENW","GBR_NIR","GBR_SCO",
        "UKR","USA")

Name2<-c("Austria","Belarus","Belgium","Bulgaria",
         "Canada","Chile","Croatia","Czechia",
         "Denmark","Estonia","Finland",
         "Germany",
         "Hungary","Iceland","Ireland","Italy",
         "Japan","Lithuania","The Netherlands","Norway",
         "Poland","Portugal","Republic of Korea",
         "Russia","Slovakia","Slovenia","Spain",
         "Sweden","Switzerland","Taiwan","United Kingdom",
         "England and Wales","Northern Ireland","Scotland",
         "Ukraine","USA") 

#CAN=5  DNK=9 IRL=17 ITA=18 
#JPN=19 NLD=21 KOR=25   RUS=26
#GBR_NP=33 GBRENW=34 GBR_SCO=36
#SWE=30  TWN=32  USA=38


N<-c(1:21,23:36)






####### Fig 1  Age-Components of TFR

B_all <- c()

for (w in 1:length(N)){

  i<-N[w]

setwd(paste("C:/Users/u1019088/DATA/HFD/",Name[i],sep=""))

A0<-read.table(paste(Name[i],"asfrRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]

A0<-A0[A0$Year>1999,]

Rg<-range(A0$Year)

A1<-A0[(A0$Year==Rg[1])|(A0$Year==Rg[2]),]
A1$Pop<-rep(Name[i],dim(A1)[1])

ASFR1<-A1

A1<-ASFR1[ASFR1$Year==Rg[1],]
A2<-ASFR1[ASFR1$Year==Rg[2],]

A1$Age <- as.numeric(gsub("[^0-9]", "", A1$Age))
A2$Age <- as.numeric(gsub("[^0-9]", "", A2$Age))

LT1<-LifeTableMx(A1$ASFR1)
LT2<-LifeTableMx(A2$ASFR1)

Age<-A1$Age

AFR11<-cumsum(Fertility1(LT1,A1))
AFR12<-cumsum(Fertility1(LT2,A2))

AFR1<-cumsum(A1$ASFR)
AFR2<-cumsum(A2$ASFR)

B1 <- data.frame(Age = Age,AFR1 = AFR11,AFR = AFR1)
B2 <- data.frame(Age = Age,AFR1 = AFR12,AFR = AFR2)

B2$PC<-LT2$lx/100000
B2$Year <- Rg[2] 
B2$Country <- Name[i]

B1$PC<-LT1$lx/100000
B1$Year <- Rg[1] 
B1$Country <- Name[i]

B_all <- rbind(B_all,B1,B2)
}

Components<-B_all

setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC/Data")

write.csv(B_all, "AgeComponents.csv", row.names = FALSE)





####### Fig 1b  Components of TFR

F <- c()

for (w in 1:length(N)){
  
  i<-N[w]
  
  setwd(paste("C:/Users/u1019088/DATA/HFD/",Name[i],sep=""))
  
  A0<-read.table(paste(Name[i],"asfrRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]
  
  setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC")
  
  A0<-A0[A0$Year>1999,]
   
  Rg<-range(A0$Year)
  F1<-c()
  
  for (t in Rg[1]:Rg[2]){
    
    if((Name[i]=="POL")&((t==2018))){t<-2019}
  A1<-A0[(A0$Year==t),]
  
  TFR<-sum(A1$ASFR)
  LT<-LifeTableMx(A1$ASFR1)
  PC<-LT$lx[44]/100000
  TFR1<-sum(Fertility1(LT1,A1))
  
  Country<-Name[i]
  Year<-t
  
  F1<- rbind(F1,cbind(Country,Year,TFR,PC,TFR1))
}
  
F<-rbind(F,F1)

}

F <- as.data.frame(F)

F$Year  <- as.numeric(F$Year)
F$TFR   <- as.numeric(F$TFR)
F$TFR1  <- as.numeric(F$TFR1)
F$PC    <- as.numeric(F$PC)

setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC/Data")

write.csv(F, "TFRComponents.csv", row.names = FALSE)



########## Figure Time trends



D <- c()

for (w in 1:length(N)){
  
  i<-N[w]
  
  setwd(paste("C:/Users/u1019088/DATA/HFD/",Name[i],sep=""))
  
  A0<-read.table(paste(Name[i],"asfrRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]
  
  setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC")
  
  A0<-A0[A0$Year>1999,]
  
  Rg<-range(A0$Year)
  
Cont<-c()
ASFR1<-A0 

for (t in 1:(Rg[2]-Rg[1])){
  
  if((Name[i]=="POL")&((t==18)|(t==19))){t<-20}
  A1<-ASFR1[ASFR1$Year==(Rg[1]+(t-1)),]
  A2<-ASFR1[ASFR1$Year==(Rg[1]+t),]
  
  LT1<-LifeTableMx(A1$ASFR1)
  LT2<-LifeTableMx(A2$ASFR1)
  
  TFR1<-sum(A1$ASFR)
  PC1<-LT1$lx[44]/100000
  
  TFR2<-sum(A2$ASFR)
  PC2<-LT2$lx[44]/100000
  
  
  TFR1n<-(TFR1)/(1-PC1)
  TFR2n<-(TFR2)/(1-PC2)
  
  
  ## changes continuous
  
  CTFR<-log(TFR2/TFR1)*sqrt(TFR1*TFR2)
  
  Term1<-log(TFR2n/TFR1n)*sqrt(TFR1n*TFR2n)*sqrt((1-PC1)*(1-PC2))
  Term2<-log(PC2/PC1)*sqrt(PC1*PC2)*sqrt(TFR1n*TFR2n)
  
  Cont<-rbind(Cont,c(CTFR,Term1,-Term2,Term1-Term2))
  
  
}

Year<-Rg[1]:(Rg[2]-1)
Eq2<-c("CTFR","Term1","Term2","Term1+2")


df2 <- reshape2::melt(Cont[,c(2,3)], c("Year","Eq3"), value.name = "values")
df2$Year<-df2$Year+Rg[1]-1
df2$Eq3<-as.character(df2$Eq3)


df_agg <- df2 %>%
  filter(Eq3 %in% c(1, 2)) %>%
  mutate(
    Period = case_when(
      Year >= 2000 & Year <= 2004 ~ "2000-2005",
      Year >= 2005 & Year <= 2009 ~ "2005-2010",
      Year >= 2010 & Year <= 2014 ~ "2010-2015",
      Year >= 2015 & Year <= 2019 ~ "2015-2020",
      Year >= 2020                ~ "2020+"
    )
  ) %>%
  group_by(Period, Eq3) %>%
  summarise(
    values = sum(values, na.rm = TRUE),
    .groups = "drop"
  )

df_agg$Country<-Name[i]
  
D<-rbind(D,df_agg)
}

setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC/Data")

write.csv(D,"TimeTrends.csv", row.names = FALSE)



############### Figure Age-decomposition



E <- c()


for (w in 1:length(N)){

  i<-N[w]

setwd(paste("C:/Users/u1019088/DATA/HFD/",Name[i],sep=""))

A0<-read.table(paste(Name[i],"asfrRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]

setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC")

A0<-A0[A0$Year>1999,]


Rg<-range(A0$Year)

A1<-A0[(A0$Year==Rg[1])|(A0$Year==Rg[2]),]
A1$Pop<-rep(Name[i],dim(A1)[1])

ASFR1<-A1

A1<-ASFR1[ASFR1$Year==Rg[1],]
A2<-ASFR1[ASFR1$Year==Rg[2],]

A1$Age <- as.numeric(gsub("[^0-9]", "", A1$Age))
A2$Age <- as.numeric(gsub("[^0-9]", "", A2$Age))

LT1<-LifeTableMx(A1$ASFR1)
LT2<-LifeTableMx(A2$ASFR1)

Age<-A1$Age

AFR11<-cumsum(Fertility1(LT1,A1))
AFR12<-cumsum(Fertility1(LT2,A2))

AFR1<-cumsum(A1$ASFR)
AFR2<-cumsum(A2$ASFR)

B1 <- data.frame(Age = Age,AFR1 = AFR11,AFR = AFR1)
B2 <- data.frame(Age = Age,AFR1 = AFR12,AFR = AFR2)

D_df<-AgeDecomp(LT1,B1,A1,LT2,B2,A2)

colSums(D_df)

D_long <- D_df |>
  pivot_longer(
    cols = c(termf1, termPC),
    names_to = "component",
    values_to = "value"
  )

D_long$component <- factor(
  D_long$component,
  levels = c("termf1", "termPC"),
  labels = c("f+1", "PC")
)
 
D_long$Country<-Name[i]

E<-rbind(E,D_long)
}


setwd("C:/Users/u1019088/OneDrive - Australian National University/Articles/TFR & PC/Data")

write.csv(E,"AgeDecomposition.csv", row.names = FALSE)
