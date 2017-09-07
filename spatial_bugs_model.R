rm(list = ls())
library(R2WinBUGS)
##Mapping
library(rgeos)
library(maptools)
library("ggplot2")
library(broom) ## for converting a map to a data frame
library(glm2)
library(ResourceSelection) ## for hosmer and lemeshow testing
library(dplyr)
library(mcmc)
library(shapefiles)
library("R2BayesX")
library(spdep)


###reading and exporting the shape file
kilifi_sub <- maptools::readShapePoly ( "data/kilif_sub_loc_Shape/DSS_subloc_Arc.shp",
                                  IDvar="Adj_ID", proj4string=CRS("+proj=longlat +ellps=clrk66"))
plot(kilifi_sub, border="red", axes=TRUE, las=1 )
###export to Splus
###then open wibug as a text file import the map
maptools::sp2WB(map =as(kilifi_sub, "SpatialPolygons"),filename="maps/kilifi_sub" )

### adjacency matrix from Arcmap using Adjacency for Winbugs Addon
num <- c(2,4,2,3,3,2,4,3,5,5,5,7,5,7,2,4,4,6,5,2,6,1,4,5,7,2,4,6,6,6,4,5,6,3,5,5,4,3,2,2)
adj <- c(2,40,
         1,3,4,40,
         2,4,
         2,3,5,
         4,7,8,
         9,12,
         5,8,10,11,
         5,7,11,
         6,10,12,13,14,
         7,9,11,13,14,
         7,8,10,13,17,
         6,9,14,15,16,18,19,
         9,10,11,14,17,
         9,10,12,13,16,17,20,
         12,19,
         12,14,18,21,
         11,13,14,20,
         12,16,19,21,24,25,
         12,15,18,24,25,
         14,17,
         16,18,23,25,29,30,
         23,
         21,22,26,29,
         18,19,25,27,28,
         18,19,21,24,27,28,30,
         23,29,
         24,25,28,31,
         24,25,27,30,31,32,
         21,23,26,30,33,34,
         21,25,28,29,32,33,
         27,28,32,35,
         28,30,31,33,35,
         29,30,32,34,35,36,
         29,33,36,
         31,32,33,36,37,
         33,34,35,37,38,
         35,36,38,39,
         36,37,39,
         37,38,
         1,2)
sumNumNeigh <- 166


### load the admissions data
admData <- read.csv("data/morbidity.csv")
admData$rain_mm <- admData$rain_mm/50
admData2 <- admData %>%  select(Adj_ID , sublocation, yr, nagem, gender , severe_disease ,
                                 cumulitive_count,cumulitive_time , EVI_VALUE ,count_adm ,
                                rain_mm,total_admission, admdays, nweight)
admData2 <- admData2[complete.cases(admData2),]
write.csv(admData2 , "data/bayesian_data_morbidity.csv")
###modelling 
dataModel <- list("N"=length(admData2$cumulitive_count) , "adj"=adj,"num"=num ,
             "sumNumNeigh" =sumNumNeigh,"cumulitive_count" =admData2$cumulitive_count,
             "EVI_VALUE"=as.numeric(admData2$EVI_VALUE) ,rainX=as.numeric(admData2$rain_mm), 
             "sublocation"=admData2$Adj_ID,"nsex"=admData2$gender ,"severe_disease"=admData2$severe_disease , 
             "nagem"=admData2$nagem ,"total_adm"=admData2$total_admission ,"nweight"=admData2$nweight,
             "admDays"=admData2$admdays, "count_adm"=admData2$count_adm)
names(dataModel)

### 
####defining the initials
# Poisson Inits
psi_s <- rep(0, times=10, each = 1)
psi_s <- c(NA,psi_s)
inits_Vals <-  function(){
  list(alpha = c(0,-0.58,0.07,-0.36,0.02,-0.20,0.34,0.11,0.017,-0.41), 
       Psi=structure(.Data =c(NA , 0  ,0 , 0 , 0 , 0,  0 , 0,  0 , 0 ,0), .Dim=c( 1, 11 )))
}

###
####defining the parameters to monitor
###param convolution
params <- c("alpha[]" , "Psi[1,]", "Phi[]" ,"r")


####use of bugs function for simulation
###codaPkg=F to brinng results to R
sink(paste0(file.path(getwd(),"/analysis/code/str_unstru/poisson_s.txt")))
modelSp <-  R2OpenBUGS::bugs(data = dataModel, model.file="model_ust_st.txt",inits=inits_Vals,
                            parameters.to.save = params,
                            n.chains=1,n.iter=20000,n.burnin=1000,n.thin=10,codaPkg=F,digits = 4,  debug = F,
                            working.directory = file.path(getwd(),"/analysis/code/str_unstru"))
print(modelSp)
sink()
plot(modelSp)
mcmcplots::mcmcplot(modelSp, dir=file.path(getwd(),"/analysis/code/str_unstru/figures") , regex = "alpha")

print(modelSp)
saveList <- list(modelSp, dataModel, adj,num, params, sumNumNeigh, admData2)


### the reported model

####use of bugs function for simulation
###codaPkg=F to brinng results to R
params <- c("alpha[]" , "Phi[]" ,"r")
inits_Vals <-  function(){
  list(alpha = c(0,-0.58,0.07,-0.36,0.02,-0.20,0.34,0.11,0.017,-0.41))
}

sink(paste0(file.path(getwd(),"/analysis/code/str_unstru/poisson_s2_se.txt")))
modelSp2 <-  R2OpenBUGS::bugs(data = dataModel, model.file="model_ust_st2.txt",inits=inits_Vals,
                             parameters.to.save = params,
                             n.chains=1,n.iter=20000,n.burnin=1000,n.thin=10,codaPkg=F,digits = 4,  
                             debug = F,
                             working.directory = file.path(getwd(),"/analysis/code/str_unstru"))
print(modelSp2)
print(modelSp2 , digits.summary=5)
sink()
plot(modelSp2)
mcmcplots::mcmcplot(modelSp2, dir=file.path(getwd(),"/analysis/code/str_unstru/figures") , regex = "alpha")

print(modelSp2 , digits.summary=5)
saveList <- list(modelSp2, dataModel, adj,num, params, sumNumNeigh, admData2)





