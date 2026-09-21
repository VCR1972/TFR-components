### This program decomposes the TFR into its fertility vs childlessness parts
### we use  life table 
library(tidyr)
library(patchwork)
library(dplyr)
options(scipen = 999)


########################################################################################
########################     FUNCTIONS     #############################################
########################################################################################


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
  
  ex<-Tx/lx 
  
  Age<-15:55    
  
  ALL<-data.frame(Age,mx,lx,dx,Lx,Tx,ex)
  return(ALL)
}


Fertility<-function(As,Y,i){
  
  A1<-As[As$Year==Rg[Y],]
  
  Age<-A1$Age
  Year<-Rg[Y]
  Country <- Name[i]
  AFR<-cumsum(A1$ASFR)
  AFR1p<-cumsum(A1$ASFR1p)
  LT<-LifeTableMx(A1$ASFR1)$lx/100000
  
  B <- data.frame(Country = Country,Year = Year,Age = Age,
                  AFR1 = AFR1p,AFR = AFR,PC = LT)
  return(B)
}


AgeDecomp<-function(b1,b2,Y){ 
  ### checking for time 1
  
  PC<-b1$PC
  AFR<-b1$AFR
  afr1<-b1$AFR1
  afr<-afr1*(1-PC)
  
  
  n<-length(afr)
  
  f1<-(afr[-1]/(1-PC[-1]))-afr1[-n]
  
  
  t1<-f1*(1-PC[-1])
  t2<-afr1[-n]*(PC[-n]-PC[-1])
    
  PC2<-b2$PC
  AFR2<-b2$AFR
  afr12<-b2$AFR1
  afr2<-afr12*(1-PC2)
  
  n<-length(afr2)
  
  f12<-(afr2[-1]/(1-PC2[-1]))-afr12[-n]
  f12[(1-PC2[-1])==0]<-0
  
  t12<-f12*(1-PC2[-1])
  t22<-afr12[-n]*(PC2[-n]-PC2[-1])
    
  ## now the decomposition kitagawa type
  
  termPC<-(afr12[-n]+afr1[-n])/2*(PC2[-n]-PC2[-1]-(PC[-n]-PC[-1]))+
    (f12+f1)*((1-PC2[-1])-(1-PC[-1]))/2
  
  ## here we use instead differences to avoid issues with the negative f1
  termf1<-(f12-f1)*((1-PC[-1])+(1-PC2[-1]))/2 +
    (afr12[-n]-afr1[-n])*((PC2[-n]-PC2[-1])+(PC[-n]-PC[-1]))/2
  
  
  Decomp<-termf1+termPC
  age<-b1$Age[-1]
  Country <- rep(Name[i],length(age))
  Year <- rep(Y,length(age))
  D<-cbind(Country,Year,age,Decomp,termf1,termPC)
  
  D_d <- as.data.frame(D)
  
  return(D_d)
}

 

Name<-c("AUT","BLR","BEL","BGR","CAN","CHL",
        "HRV","CZE","DNK","EST","FIN",
        "DEUTNP","DEUTE","DEUTW","HUN","ISL","IRL",
        "ITA","JPN","LTU","NLD","NOR","POL","PRT",
        "KOR","RUS","SVK","SVN","ESP","SWE","CHE",
        "TWN","GBR_NP","GBRTENW","GBR_NIR","GBR_SCO",
        "UKR","USA")

Name2<-c("Austria","Belarus","Belgium","Bulgaria",
         "Canada","Chile","Croatia","Czechia",
         "Denmark","Estonia","Finland",
         "Germany","Germany East","Germany West",
         "Hungary","Iceland","Ireland","Italy",
         "Japan","Lithuania","The Netherlands","Norway",
         "Poland","Portugal","Republic of Korea",
         "Russia","Slovakia","Slovenia","Spain",
         "Sweden","Switzerland","Taiwan","United Kingdom",
         "England and Wales","Northern Ireland","Scotland",
         "Ukraine","United States") 


ALL<-c()
## selected countries without data issues 
N<-c(1,3:5,7:12,15:25,27:36,38)


########################################################################################
########################        DATA       #############################################
########################################################################################



####   first preparing the data, and creating the new first birth age-specific fertility rates


for (w in 1:length(N)){
   i<-N[w]
   
  setwd(paste(".../DATA/HFD/",Name[i],sep=""))
  
  A0<-read.table(paste(Name[i],"asfrRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]
  
  B1<-read.table(paste(Name[i],"birthsRRbo.txt",sep=""),header=TRUE, skip=2)[,c(1:4)]
  
  E1<-read.table(paste(Name[i],"exposRR.txt",sep=""),header=TRUE, skip=2)
  
  setwd(".../Articles/TFR & PC/Data/NewASFR1")
  
## we fill the data with the information from the first year  
    
  B1$Age <- as.character( B1$Age)
  B1$Age[ B1$Age=="12-"] <- "12"
  B1$Age[ B1$Age=="55+"] <- "55"
  B1$Age <- as.numeric( B1$Age)
  
  
  E1$Age <- as.character( E1$Age)
  E1$Age[ E1$Age=="12-"] <- "12"
  E1$Age[ E1$Age=="55+"] <- "55"
  E1$Age <- as.numeric( E1$Age)
  
  BE1 <- inner_join(
    B1,
    E1,
    by = c("Year", "Age")
  ) 
   
  BE1<-BE1[BE1$Year>1900,]  # instead this gives you all the twentieth century 
  BE1<-BE1[BE1$Age>14,]     # which can be useful to compare period and cohort
  
  BE1 <- BE1 %>%
    group_by(Year) %>%
    mutate( asf1 = B1/Exposure,
      B1_cum_lag = lag(cumsum(asf1), default = 0)) %>%
    ungroup()
  
  
  #and create a new database that calculates the
  # ASFR for first births with the exposure 
  # the total number of women minus the number of 
  # first births accumulated and half of those in each age
  
  
  YA <- BE1 %>% 
    group_by(Year) %>%
    mutate(
      Ratio = B1_cum_lag + (B1/(Exposure*2)),
      ASFR1 = B1 / (Exposure * (1 - Ratio)),
      ASFR = Total/Exposure,
      lx = LifeTableMx(ASFR1)$lx,
      ASFR1pp = if_else(Exposure * Ratio>0,((Total - B1) / (Exposure * Ratio)),0),
      ASFR1p = ASFR1pp*(sum(ASFR)/(1-(lx[41]/100000)))/sum(ASFR1pp)
    ) %>%
    ungroup()
  
   
  
  YA[is.nan(YA$ASFR1p),]$ASFR1p<-0
  YA[YA$ASFR1p == Inf, "ASFR1p"] <- 0
  YA[YA$Age==12,]$ASFR1p<-0
  
      

TFR<-colSums(matrix(YA$ASFR,(55-15+1)))
TFR1<-colSums(matrix(YA$ASFR1p,(55-15+1)))
PC<-1-(TFR/TFR1)
Year<-unique(YA$Year)
Population<-rep(Name[i],length(unique(YA$Year)))
PCb<-YA[YA$Age==55,]$lx/100000
TFR1b<-TFR/(1-PCb)

ALL<-rbind(ALL,cbind(Population,Year,TFR,TFR1,TFR1b,PC,PCb))


  
YAs<-data.frame(cbind(YA$Year,YA$Age,YA$ASFR,YA$ASFR1,YA$ASFR1p))
colnames(YAs)<-c("Year","Age","ASFR","ASFR1","ASFR1p")
B1<-c()
E1<-c()
A0<-c()
write.csv(YAs,paste(Name[i],"NewASFR1b.txt",sep=""), row.names = FALSE)
}








########################################################################################
########################  NOW THE FIGURES  #############################################
########################################################################################




####### Figure   Components of TFR

F <- c()


for (w in 1:length(N)){
  
  i<-N[w]
  
  setwd("..... /Articles/TFR & PC/Data/NewASFR1")
  
  A0<-read.table(paste(Name[i],"NewASFR1.txt",sep=""),header=TRUE,sep=",")
  
  setwd("...../Articles/TFR & PC/Results")
  
  A0<-A0[A0$Year>1999,]
   
  Rg<-range(A0$Year)
  
 
  F1<-c()
  for (t in Rg[1]:Rg[2]){
     if((Name[i]=="POL")&((t==2018))){t<-2019}
 
    A1<-A0[(A0$Year==t),]
    
    TFR<-sum(A1$ASFR)
    LT<-LifeTableMx(A1$ASFR1)
    PC<-LT$lx[41]/100000
    TFR1<-TFR/(1-PC)
    
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

setwd("..... /Articles/TFR & PC/shiny")

write.csv(F, "TFRComponents.csv", row.names = FALSE)









########## Figure Time trends


D <- c()

for (w in 1:length(N)){
  
  i<-N[w]
  
  setwd(".... /Articles/TFR & PC/Data/NewASFR1")
  
  A0<-read.table(paste(Name[i],"NewASFR1.txt",sep=""),header=TRUE,sep=",")
  
  setwd("..../Articles/TFR & PC/Results")
  
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
    PC1<-LT1$lx[41]/100000
    
    TFR2<-sum(A2$ASFR)
    PC2<-LT2$lx[41]/100000
    
    
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
  
  
  df2 <- reshape2::melt(Cont[,-1], c("Year","Eq2"), value.name = "values")
  df2$Year<-df2$Year+Rg[1]-1
  df2$Eq2<-as.character(df2$Eq2)
  
  
  df_agg <- df2 %>%
    filter(Eq2 %in% c(1, 2)) %>%
    mutate(
      Period = case_when(
        Year >= 2000 & Year <= 2004 ~ "2000-2005",
        Year >= 2005 & Year <= 2009 ~ "2005-2010",
        Year >= 2010 & Year <= 2014 ~ "2010-2015",
        Year >= 2015 & Year <= 2019 ~ "2015-2020",
        Year >= 2020                ~ "2020+"
      )
    ) %>%
    group_by(Period, Eq2) %>%
    summarise(
      values = sum(values, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_agg$Country<-Name[i]
  
  D<-rbind(D,df_agg)
}


D <- as.data.frame(D)

D$Eq2  <- as.numeric(D$Eq2)

setwd("..../Articles/TFR & PC/shiny")

write.csv(D,"TimeTrends.csv", row.names = FALSE)











####### Figure  Age-Components of TFR

B_all <- c()

for (w in 1:length(N)){

  i<-N[w]

  setwd("..../Articles/TFR & PC/Data/NewASFR1")
  
  A0<-read.table(paste(Name[i],"NewASFR1.txt",sep=""),header=TRUE,sep=",")
  
  setwd("..../Articles/TFR & PC/Results")
  
A0<-A0[A0$Year>1999,]
Rg<-range(A0$Year)

B1<-Fertility(A0,1,i)
B2<-Fertility(A0,2,i)

B_all <- rbind(B_all,B1,B2)
}

Components<-as.data.frame(B_all)


setwd("..../Articles/TFR & PC/shiny")

write.csv(B_all, "AgeComponents.csv", row.names = FALSE)











############### Figure Age-decomposition



E <- c()


for (w in 1:length(N)){
  
  E2 <- c()
  
  i<-N[w]

  setwd("..../Articles/TFR & PC/Data/NewASFR1")
  
  A0<-read.table(paste(Name[i],"NewASFR1.txt",sep=""),header=TRUE,sep=",")
  
  setwd("..../Articles/TFR & PC/Results")
  
A0<-A0[A0$Year>1999,]

Rg<-range(A0$Year)
Rg<-c(Rg[1],2010)

B1<-Fertility(A0,1,i)
B2<-Fertility(A0,2,i)

E2<-rbind(E2,AgeDecomp(B1,B2,Rg[1]))


Rg<-range(A0$Year)
Rg<-c(2010,Rg[2])

B1<-Fertility(A0,1,i)
B2<-Fertility(A0,2,i)

E2<-rbind(E2,AgeDecomp(B1,B2,Rg[2]))

D_long <- E2 |>
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


E <- as.data.frame(E)

E$Year  <- as.numeric(E$Year)
E$age  <- as.numeric(E$age)
E$Decomp  <- as.numeric(E$Decomp)
E$value  <- as.numeric(E$value)


E <-E %>%
  mutate(
    Period = case_when(
      Year >= 2000 & Year < 2010 ~ paste(Year,"-2010",sep=""),
      Year >= 2010 & Year <= 2026 ~ paste("2010-",Year,sep="")
    ) )


setwd("..../Articles/TFR & PC/shiny")

write.csv(E,"AgeDecomposition.csv", row.names = FALSE)
