library(haven)
library(dplyr)
library(foreign)
library(tidyr)
###for muted function in ggplot2
library(scales)
##Mapping
library(rgeos)
library(maptools)
library("ggplot2")
library(broom) ## for converting a map to a data frame
library(glm2)
library(ResourceSelection) ## for hosmer and lemeshow testing
library(lme4)
library(boot)

admData <- haven::read_dta("mortality.dta" )
admData$timeR <- admData$time_to_readmission /30
admData$rainX <- admData$rain_mm /50
admData$age_yr <- admData$nagem/12

##Model

tmod <-   lme4::glmer(noutcome ~ whz_all + as.factor(nsex) + as.factor(severe_disease) +
                        EVI_VALUE + rainX +
                        total_admission +  admdays  + timeR  + age_yr +(1|sublocation),
                      family="binomial",
                      control=glmerControl(optimizer="bobyqa"),     data= admData)


summary(tmod)
tmod@beta[2]


dataSam <- admData %>% select(whz_all , noutcome,malnutKid,yr , sublocation)

dataSam$rr <- exp(tmod@beta[2] * dataSam$whz_all )

dataMalnut <- dataSam %>% 
  filter(dataSam$noutcome==1) %>% 
  select(malnutKid, rr, noutcome, sublocation  , whz_all)
#subset of those that died
N <- nrow(dataMalnut)         #total no. of death cases

lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion



## function utlised by the bootstap to calculate the confidence interval
AF_fraction2 <- function(data=data , indices){
  datax <- data[indices,] 
  #bTau <- bestTau.fit(datax)
  #datax$whz06tau <- (abs((datax$whz06new))^bTau )*sign(datax$whz06new)
  tmod <-   lme4::glmer(noutcome ~ whz_all + as.factor(nsex) + as.factor(severe_disease) +
                          EVI_VALUE + rainX +
                          total_admission +  admdays  + timeR  + age_yr +(1|sublocation),
                        family="binomial",
                        control=glmerControl(optimizer="bobyqa"),     data= datax)
    #RR
  datax$rr <- exp(tmod@beta[2] * datax$whz_all ) 
  dataMalnut <- datax %>% 
    filter(datax$noutcome==1) %>% 
    select(malnutKid, rr, noutcome,  whz_all , nlocation2)
  #subset of those that died
  N <- nrow(dataMalnut)         #total no. of death cases
  lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion
  print(lambda)
}

x <- boot(data=admData, statistic=AF_fraction2, R=300 ,  sim = "ordinary")
y <- boot.ci(x, type="basic" , index = 1) # intercept 





atrributable <- data.frame()
for (i in unique(dataMalnut$sublocation)) {
  dataMalnut2 <- dataMalnut %>% filter(sublocation==i)
  N <- nrow(dataMalnut2)  
 x <-  1/N*sum((dataMalnut2$rr-1)/dataMalnut2$rr, na.rm = T) 
 dat <- data.frame(af=x, subloc=i)
atrributable <- rbind(dat,  atrributable)
 print(paste0(i , "-" , x))
}
write.csv(atrributable , "atributable.csv")
# 



## location attributable fractions
atrribut_loc <- data.frame()
for (x in unique(admData$nlocation2)) {
  admData2 <- admData %>% filter(nlocation2==x)
  admData2$nlocation23 <- as.numeric(as.factor(admData2$nlocation2))
  if (length(table(admData2$sublocation))[1]>1) {
tmod2 <-   lme4::glmer(noutcome ~ whz_all + as.factor(nsex) + as.factor(severe_disease) +
                         EVI_VALUE + rainX +
                         total_admission +  admdays  + timeR  +(1|sublocation),
                       family="binomial",control=glmerControl(optimizer="bobyqa"),     data= admData2)




dataSam <- admData2 %>% select(whz_all , noutcome,malnutKid,yr , sublocation)

dataSam$rr <- exp(tmod2@beta[2] * dataSam$whz_all )

dataMalnut <- dataSam %>% 
  filter(dataSam$noutcome==1) %>% 
  select(malnutKid, rr, noutcome, sublocation  , whz_all)
N <- nrow(dataMalnut) 
lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion
  dat <- data.frame(af=lambda, loc=x)
atrribut_loc <- rbind(dat,  atrribut_loc)
  }
  else {
    print(x)
  }
}

write.csv(atrribut_loc , "atributable_loc.csv")


### fitting a winbugs model
  