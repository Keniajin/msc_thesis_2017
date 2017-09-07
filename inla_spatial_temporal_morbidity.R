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
  
  ### load the admissions data
  admData <- read.csv("data/morbidity.csv")
  admData$rain_mm <- admData$rain_mm/50
  admData$severe_disease <- factor(admData$severe_disease , levels=c(0,1,2,3))
  admData <- admData %>% mutate(gender2= ifelse(gender==1 ,0,1 ))

  admData$gender2 <- factor(admData$gender2 , levels=c(0,1) )
  admData$gender <- admData$gender2
  admData2 <- admData %>%  dplyr::select(Adj_ID , sublocation,mnth, nagem, gender , 
                                         severe_disease ,
                                  cumulitive_count,cumulitive_time , EVI_VALUE ,count_adm ,rain_mm,
                                  total_admission ,admdays ,nweight ,yr)
  
  ###generate othe variables to be used within INLa
  admData2$Adj_ID2 <- admData2$Adj_ID
  admData2$Adj_ID3 <- admData2$Adj_ID
  admData2$count_adm2 <- admData2$count_adm
  admData2$count_adm3 <- admData2$count_adm
  admData2$count_adm4 <- admData2$count_adm
  admData2$count_adm5 <- admData2$count_adm
  admData2$count_adm6 <- admData2$count_adm
  admData2$count_adm7 <- admData2$count_adm
  admData2$EVI_VALUE2 <- admData2$EVI_VALUE
  admData2$rain_mm2 <- admData2$rain_mm
  admData2$nagem2 <- admData2$nagem
  admData2$severe_disease2 <- as.factor(admData2$severe_disease)
  admData2$mnth2 <- admData2$mnth
  admData2$nweight2 <-   admData2$nweight
  admData2$admdays2 <-   admData2$admdays
  
  ###defining the priors
  prior.iid = c(1,0.01)
  prior.besag = c(1,0.001)
  initial.iid = 4
  initial.besag = 3

  #'         \\\_Model_1_\\\         #  
###############MOdel 1 #######################################################    
  ### replicate the final model in STATA
  ## spatial unstructured
  ## DIC 11229.11
  formulaUH0 <- cumulitive_count ~ EVI_VALUE + rain_mm  + gender + severe_disease +
    total_admission  + admdays + nweight +
    f(Adj_ID, model = "iid",prior="normal",param=c(0, 0.001) , initial = 1)
  
  resultUH0 <- inla(formulaUH0,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE,cpo=TRUE),E=log(nagem) ,
                   control.predictor = list(compute = TRUE))
  ##summary in 3 decimal places
  summary(resultUH0)
  exp(resultUH0$summary.fixed)
  write.csv(data.frame(resultUH0$summary.fixed), "results1_14504_36.csv")
  
  pdresultUH0 <- resultUH0$dic$p.eff
  
  ###resultUH0 22.03 (14504.36) pd and DIC
  
 #'         \\\_Model_1_\\\         #  
 #'         #
  
  #'         \\\_Model 1B\\\         #  
  ###############MOdel 1B#######################################################  
  ### spatial model structured and unstrustured  without 
  ### to comapare with Winbugs
  ### DIC 14498.08
  ## PD 25.03
  formulaUHB <- cumulitive_count ~ EVI_VALUE + rain_mm   + gender + severe_disease +
    total_admission + admdays + nweight +
    f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
      hyper=list(prec.unstruct=list(prior="loggamma",param=c(0.0111,0.001)),
                 prec.spatial=list(prior="loggamma",param=c(0.0011,0.001)))) 
  
  # f(count_adm, model = "ar1", replicate = Adj_ID3)
  resultUHB <- inla(formulaUHB,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE,cpo=TRUE),E=log(nagem)
                   ,control.predictor(compute=TRUE))
  summary(resultUHB)
  pdresultUHB <- resultUHB$dic$p.eff #25.03
  exp(resultUHB$summary.fixed)
  write.csv(data.frame(resultUHB$summary.fixed), "results2_14498.08.csv")
  
  #write.csv(data.frame(resultUHB$summary.fixed), "results_20.05_under5_10700.63.csv")
  
  ####The   computation of the posterior mean for the random effects 𝝃 is performed in two
  # steps as we have more than one parameter:
  # we extract the marginal posterior distribution for each element of the random effect
  csi <- resultUHB$marginals.random$Adj_ID[1:40]
  
  ## then apply the exponential transformation and calculate the posterior mean for each of   them using the lapply function.
  zeta <- lapply(csi,function(x) inla.emarginal(exp,x))
  ##define the cut offs for your risk ratio
  zeta.cutoff <- c(0.9, 0.95, 0.999 ,1.0,1.01,1.05, 1.1)
  
  #Transform zeta in categorical variable
  cat.zeta <- cut(unlist(zeta),breaks=zeta.cutoff,
                  include.lowest=TRUE )
  
  #Create a dataframe with all the information needed for the map
  maps.cat.zeta <- data.frame(unique(admData2$Adj_ID), cat.zeta=cat.zeta)
  
  #Add the categorized zeta to the kilifi spatial polygon
  ## 
  data.kilifi <- attr(kilifi_sub, "data")
  attr(kilifi_sub, "data") <- merge(data.kilifi, maps.cat.zeta,
                                    by.x="Adj_ID" , by.y="unique.admData2.Adj_ID.")
  
  ## mapping the risk ratio 
  #spplot(obj=kilifi_sub, zcol= "cat.zeta", col.regions=gray(seq(0.9,0.1,length=4)), asp=1)
  spplot(obj=kilifi_sub, zcol= "cat.zeta",col.regions=diverge_hsv(8), scales=list(draw = TRUE), asp=1)
  
  
#'         \\\_Model 2\\\         #  
###############MOdel 2#######################################################  
### spatial model structured and unstrustured  with the temporal component included
### fitting model 1 
### DIC 13640.21
## pD 35.50
admData2$nagem_int <- as.integer(admData2$nagem)
  
formulaUH <- cumulitive_count ~ EVI_VALUE + rain_mm   + gender + severe_disease + total_admission + admdays + nweight +
    f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,hyper=list(prec.unstruct=list(prior="loggamma",param=c(0.0111,0.001)),
                 prec.spatial=list(prior="loggamma",param=c(0.0011,0.001)))) + f(count_adm, model = "ar1")
  
# f(count_adm, model = "ar1", replicate = Adj_ID3)
resultUH <- inla(formulaUH,family="nbinomial",
                 data=admData2, control.compute=list(dic=TRUE,cpo=TRUE),E=log(nagem_int)
                 ,control.predictor(compute=TRUE))
summary(resultUH)
pdresultUH <- resultUH$dic$p.eff #35.50
exp(resultUH$summary.fixed)
write.csv(data.frame(resultUH$summary.fixed), "nm_results2_13640.2.csv")

#write.csv(data.frame(resultUH$summary.fixed), "nm_results28.75_under5_9896.31.csv")

####The   computation of the posterior mean for the random effects 𝝃 is performed in two
# steps as we have more than one parameter:
# we extract the marginal posterior distribution for each element of the random effect
csi <- resultUH$marginals.random$Adj_ID[1:40]

## then apply the exponential transformation and calculate the posterior mean for each of   them using the lapply function.
zeta <- lapply(csi,function(x) inla.emarginal(exp,x))
##define the cut offs for your risk ratio
zeta.cutoff <- c(0.83,0.9, 0.95, 0.999 ,1.0,1.01,1.05, 1.1 ,1.2)

#Transform zeta in categorical variable
cat.zeta <- cut(unlist(zeta),breaks=zeta.cutoff,
                include.lowest=TRUE )

#Create a dataframe with all the information needed for the map
maps.cat.zeta <- data.frame(unique(admData2$Adj_ID), cat.zeta=cat.zeta)

#Add the categorized zeta to the kilifi spatial polygon
## 
data.kilifi <- attr(kilifi_sub, "data")
attr(kilifi_sub, "data") <- merge(data.kilifi, maps.cat.zeta,
                                  by.x="Adj_ID" , by.y="unique.admData2.Adj_ID.")

## mapping the risk ratio 
#spplot(obj=kilifi_sub, zcol= "cat.zeta", col.regions=gray(seq(0.9,0.1,length=4)), asp=1)
png(filename=paste0("figure4A","img.png") , width = 19.45 , height =  22.40 , units = "cm" , res=300)

spplot(obj=kilifi_sub, zcol= "cat.zeta",col.regions=diverge_hsv(8), scales=list(draw = TRUE), asp=1)
dev.off()

### temporal graph
# csiz <- resultUH$marginals.random$count_adm[1:11]
# zetaz <- lapply(csiz,function(x) inla.emarginal(exp,x))
# x <- seq(1,11) # Years
# plot(x,zetaz, type="l", main="",xlab="Count of admissions",ylab=expression(beta*t))

plot( resultUH, plot.fixed.effects = TRUE, constant=FALSE,
      plot.lincomb = TRUE, 
      plot.random.effects = TRUE, 
      plot.hyperparameters = TRUE,
      plot.predictor = TRUE, 
      plot.q = TRUE, 
      plot.cpo = TRUE,
      single = TRUE)

  plot( resultUH, plot.fixed.effects = TRUE , constant=FALSE,plot.cpo = F,single =F)
  
save.image("stModel.RDA")  

#'         \\\_Model 3\\\         #
#'         ###############MOdel With variables changing over time#######################################################  
#### Fitting a SPATIAL Temporal Model
formulaUH2b <- cumulitive_count ~  EVI_VALUE  + gender + 
  severe_disease + total_admission + rain_mm + admdays + nweight +
  f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
    hyper=list(prec.unstruct=list(prior="loggamma",param=c(0.001,0.001)),
               prec.spatial=list(prior="loggamma",param=c(0.1,0.01))))+ 
  f(EVI_VALUE2 , count_adm2, model = "iid") +
  f(rain_mm2 , count_adm3, model = "iid") +
  f(nweight2 , count_adm5, model = "iid") +
  f(admdays2 , count_adm6, model = "iid") +
  f( count_adm7, model = "ar1")

### added due to heissan values errors 
##https://groups.google.com/forum/#!topic/r-inla-discussion-group/rTdjAnILdnM
resultUH2b <- inla(formulaUH2b,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE),control.predictor(compute=TRUE) ,
                   control.inla = list(tolerance = 1e-20, h = 1e-08),E=log(nagem))
pdresultH2 <- resultUH2b$dic$p.eff #447.2864


summary(resultUH2b)
pdresultUH2b <- resultUH2b$dic$p.eff
write.csv(data.frame(resultUH2b$summary.fixed), "resultsST_10296.73.csv")

csi2 <- resultUH2b$marginals.random$Adj_ID[1:40]

## then apply the exponential transformation and calculate the posterior mean for each of   them using the lapply function.
zeta2 <- lapply(csi2,function(x) inla.emarginal(exp,x))
##define the cut offs for your risk ratio
zeta.cutoff2 <- c(0.8,0.99, 1.0,1.001,1.1, 1.2)


#Transform zeta in categorical variable
cat.zeta2 <- cut(unlist(zeta2),breaks=zeta.cutoff2,
                include.lowest=TRUE)

#Create a dataframe with all the information needed for the map
maps.cat.zeta2 <- data.frame(unique(admData2$Adj_ID), cat.zeta2=cat.zeta2)

#Add the categorized zeta to the kilifi spatial polygon
## 
data.kilifi2 <- attr(kilifi_sub, "data")
attr(kilifi_sub, "data") <- merge(data.kilifi2, maps.cat.zeta2,
                                  by.x="Adj_ID" , by.y="unique.admData2.Adj_ID.")

## mapping the risk ratio 
#spplot(obj=kilifi_sub, zcol= "cat.zeta2",
       #col.regions=gray(seq(0.9,0.1,length=4)), asp=1)
spplot(obj=kilifi_sub, zcol= "cat.zeta2",col.regions=diverge_hsv(8), scales=list(draw = TRUE), asp=1)
### temporal graph
# csiz <- resultUH$marginals.random$count_adm[1:11]
# zetaz <- lapply(csiz,function(x) inla.emarginal(exp,x))
# x <- seq(1,11) # Years
# plot(x,zetaz, type="l", main="",xlab="Count of admissions",ylab=expression(beta*t))
# 


plot( resultUH2b, plot.fixed.effects = TRUE , constant=FALSE,plot.cpo = F,single =F)
plot( resultUH2b, plot.fixed.effects = TRUE, constant=FALSE,
      plot.lincomb = TRUE, 
      plot.random.effects = TRUE, 
      plot.hyperparameters = TRUE,
      plot.predictor = TRUE, 
      plot.q = TRUE, 
      plot.cpo = TRUE,
      single = TRUE)


### others
###
mat.marg <- matrix(NA, nrow=40, ncol=100000)
m <- resultUH2$marginals.random$Adj_ID

for (i in 1:40){
  #Remember that the first Nareas values of the random effects
  #are u+v, while u values are stored in the (Nareas+1) to
  #(2*Nareas) elements.
  u <- m[[40+i]]
  mat.marg[i,] <- inla.rmarginal(100000, u)
}

var.u <- apply(mat.marg, 2, var)

var.v <- inla.rmarginal(100000,inla.tmarginal(function(x) 1/x,
                                              resultUH2$marginals.hyper$"Precision for ID (iid component)"))


##default priors intercept ~ N(0,0) other fixed effects ~ N(0, 0.001) 
  formulaUH2 <- cumulitive_count ~ EVI_VALUE  + gender + severe_disease + total_admission + rain_mm + admdays + nweight +
    f(Adj_ID, model = "iid",prior="normal",param=c(0, 0.001) , initial = 0) +  
    f(Adj_ID2,model="besag",graph=klf.adj , initial = 1) + 
    f(count_adm, model = "ar1") +
    f(yr, EVI_VALUE2, model="rw2" , initial = -0.38 , prior = "normal", param = c(0.38,0.01) ) +
    f(mnth2, rain_mm2, model="rw2" ,  initial = 0.076 , prior = "normal", param = c(0.63,0.88)) +
    f(count_adm3, severe_disease, model="rw2" , initial = 0)
  
 
  
  
#### model 2 - fitting model 2
  ##default priors intercept ~ N(0,0) other fixed effects ~ N(0, 0.001) 
  formulaUH <- cumulitive_count ~ EVI_VALUE  + gender + severe_disease + total_admission + rain_mm + admdays + nweight +
    f(Adj_ID, model = "iid",prior="normal",param=c(0, 0.001) , initial = 0) +  
    f(Adj_ID2,model="besag",graph=klf.adj , initial = 1) + 
    f(count_adm, model = "ar1")  +
    f(yr, EVI_VALUE2, model="rw2" , initial = -0.38 , prior = "normal", param = c(0.38,0.01) ) +
    f(mnth2, rain_mm2, model="rw2" ,  initial = 0.076 , prior = "normal", param = c(0.63,0.88)) +
    f(count_adm3, severe_disease, model="rw2" , initial = 0)
  
  # f(count_adm, model = "ar1", replicate = Adj_ID3)
  resultUH <- inla(formulaUH,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE,cpo=TRUE),E=log(nagem))
  summary(resultUH)
  
  
  # 
  # ## Provide improved estimates for the hyperparameters
  # # hyper = inla.hyperpar(resultUH)
  # # summary(hyper)
  # # plot(hyper)
  # 
  # exp(resultUH$summary.fixed)
  # exp(resultUH)
  # names(resultUH)
  #   RE1<-resultUH$summary.random$Adj_ID[,2]
  # plot(density(RE1))
  # 
  # 
  # 
  # exp.fixed.maginal <- inla.tmarginal(function(x) exp(x), result1.bym
  #                                  $marginals.fixed[[1]]) 
  # plot(exp.fixed.maginal) 
  # 
  # 

  

  
  
  
  