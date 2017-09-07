set.seed(1221223)
rm(list = ls())

library(R2WinBUGS)
##Mapping
library(rgeos)
library(maptools)
library("ggplot2")
library(broom) ## for converting a map to a data frame
#library(glm2)
#library(ResourceSelection) ## for hosmer and lemeshow testing
library(dplyr)
library(INLA)
library(spdep)
## coloring the spplot
library(colorspace)

###reading and exporting the shape file
kilifi_sub <- maptools::readShapePoly ( "data/kilif_sub_loc_Shape/DSS_subloc_Arc.shp",
                                        IDvar="Adj_ID", proj4string=CRS("+proj=longlat +ellps=clrk66"))

temp <- spdep::poly2nb(kilifi_sub)
nb2INLA("data/kilif_sub_loc_Shape/DSS_subloc_Arc.graph", temp)
klf.adj <- paste(getwd(),"/data/kilif_sub_loc_Shape/DSS_subloc_Arc.graph",sep="")



### model for mortality
admData <- haven::read_dta("data/mortality.dta" )
admData$timeR <- admData$time_to_readmission /30
admData$rainX <- admData$rain_mm /50
admData$age_yr <- admData$nagem/12
admData$whz_all <- admData$whz06
admData$whz_all[is.na(admData$whz_all)] <- admData$`_zbfa`[is.na(admData$whz_all)]
#xx <- admData %>% select(whz06,whz_all , `_zbfa`, serialno)

#'         \\\_Model 0\\\         #  
###############MOdel 0#######################################################  
### spatial model structured and unstrustured  with the temporal component included
### fitting model 0 
### DIC 7690.56
##pd =26,37
formulaMorta0 <- noutcome ~ whz_all + as.factor(gender) + as.factor(severe_disease) +
  EVI_VALUE + rainX + total_admission +  admdays  + timeR + age_yr +
  f(Adj_ID, model = "iid",prior="normal",param=c(0, 0.001) , initial = 1)

# f(count_adm, model = "ar1", replicate = Adj_ID3)
resultMorta0 <- inla(formulaMorta0,family="binomial",
                    data=admData, control.compute=list(dic=TRUE,cpo=TRUE), 
                    control.predictor(compute=TRUE))
summary(resultMorta0)
write.csv(data.frame(resultMorta0$summary.fixed), "resultsMortality_7690_56.csv")


#'         \\\_Model 1\\\         #  
###############MOdel 1#######################################################  
### spatial model structured and unstrustured  with the temporal component included
### fitting model 1 
### DIC 7690.93
##pd 22.41
formulaMorta <- noutcome ~ whz_all + as.factor(gender) + as.factor(severe_disease) +
  EVI_VALUE + rainX + total_admission +  admdays  + timeR + age_yr +
f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
  hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),
             prec.spatial=list(prior="loggamma",param=c(1,0.001))))

# f(count_adm, model = "ar1", replicate = Adj_ID3)
resultMorta <- inla(formulaMorta,family="binomial",
                    data=admData, control.compute=list(dic=TRUE,cpo=TRUE), 
                    control.predictor(compute=TRUE))
summary(resultMorta)
plot( resultMorta, plot.fixed.effects = TRUE , constant=FALSE,plot.cpo = F,single =F)

exp(resultMorta$summary.fixed)
write.csv(data.frame(resultMorta$summary.fixed), "resultsMortality_7690.93.csv")

beta <- resultMorta$summary.fixed[2,1]
dataSam <- admData %>% select(whz_all , noutcome,malnutKid,yr , sublocation)
dataSam$rr <- exp(beta* dataSam$whz_all )


dataMalnut <- dataSam %>% 
  filter(dataSam$noutcome==1) %>% 
  select(malnutKid, rr, noutcome, sublocation  , whz_all)
#subset of those that died
N <- nrow(dataMalnut)         #total no. of death cases

lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion
atrributable <- data.frame()
for (i in unique(dataMalnut$sublocation)) {
  dataMalnut2 <- dataMalnut %>% filter(sublocation==i)
  N <- nrow(dataMalnut2)  
  x <-  1/N*sum((dataMalnut2$rr-1)/dataMalnut2$rr, na.rm = T) 
  dat <- data.frame(af=x, subloc=i)
  atrributable <- rbind(dat,  atrributable)
  print(paste0(i , "-" , x))
}
write.csv(atrributable , "data/atributable_inla.csv")
# 



 
## function utlised by the bootstap to calculate the confidence interval
AF_fraction2 <- function(data=data , indices){
  datax <- data[indices,] 
  #bTau <- bestTau.fit(datax)
  #datax$whz06tau <- (abs((datax$whz06new))^bTau )*sign(datax$whz06new)
  formulaMorta <- noutcome ~ whz_all + as.factor(gender) + as.factor(severe_disease) +
    EVI_VALUE + rainX + total_admission +  admdays  + timeR + age_yr +
    f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
      hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),
                 prec.spatial=list(prior="loggamma",param=c(1,0.001))))
  
  # f(count_adm, model = "ar1", replicate = Adj_ID3)
  resultMorta <- inla(formulaMorta,family="binomial",
                      data=datax, control.compute=list(dic=TRUE,cpo=TRUE), 
                      control.predictor(compute=TRUE))
  #RR
  beta <- resultMorta$summary.fixed[2,1]
  datax$rr <- exp(beta * datax$whz_all ) 
  dataMalnut <- datax %>% 
    filter(datax$noutcome==1) %>% 
    select(malnutKid, rr, noutcome,  whz_all , nlocation2)
  #subset of those that died
  N <- nrow(dataMalnut)         #total no. of death cases
  lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion
  print(lambda)
}
x <- boot(data=admData, statistic=AF_fraction2, R=20 ,  sim = "ordinary")
y <- boot.ci(x, type="basic" , index = 1) # intercept




## function utlised by the bootstap to calculate the confidence interval
AF_fraction2 <- function(data=data , indices){
  datax <- data[indices,] 
  #bTau <- bestTau.fit(datax)
  #datax$whz06tau <- (abs((datax$whz06new))^bTau )*sign(datax$whz06new)
  formulaMorta <- noutcome ~ whz_all + as.factor(gender) + as.factor(severe_disease) +
    EVI_VALUE + rainX + total_admission +  admdays  + timeR + age_yr +
    f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
      hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),
                 prec.spatial=list(prior="loggamma",param=c(1,0.001))))
  
  # f(count_adm, model = "ar1", replicate = Adj_ID3)
  resultMorta <- inla(formulaMorta,family="binomial",
                      data=datax, control.compute=list(dic=TRUE,cpo=TRUE), 
                      control.predictor(compute=TRUE))
  #RR
  beta <- resultMorta$summary.fixed[2,1]
  datax$rr <- exp(beta * datax$whz_all ) 
  dataMalnut <- datax %>% 
    filter(datax$noutcome==1) %>% 
    select(malnutKid, rr, noutcome,  whz_all , nlocation2)
  #subset of those that died
  N <- nrow(dataMalnut)         #total no. of death cases
  lambda <- 1/N*sum((dataMalnut$rr-1)/dataMalnut$rr, na.rm = T) #proportion of death cases attributable to malnutrion
  print(lambda)
}

x2 <- boot(data=admData, statistic=AF_fraction2, R=20 ,  sim = "ordinary")
y2 <- boot.ci(x, type="basic" , index = 1) # intercept 
