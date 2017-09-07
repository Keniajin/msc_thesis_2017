for ( i in unique(admData2$count_adm)) {
  rm(list = ls())
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
  
  ###
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
  
  
  admData2x <- admData2 %>% filter(count_adm==i)
  formulaUH <- cumulitive_count ~ EVI_VALUE + rain_mm   + 
  gender + severe_disease + total_admission + admdays + nweight +
  f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
    hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),
               prec.spatial=list(prior="loggamma",param=c(1,0.001)))) 

# f(count_adm, model = "ar1", replicate = Adj_ID3)
  resultUH <- inla(formulaUH,family="nbinomial",
                 data=admData2x, control.compute=list(dic=TRUE,cpo=TRUE),E=log(nagem)
                 ,control.predictor(compute=TRUE))

#summary(resultUH)

#exp(resultUH$summary.fixed)
  write.csv(data.frame(exp(resultUH$summary.fixed)), paste0(i,"_results2_10504.53.csv"))
####The   computation of the posterior mean for the random effects 𝝃 is performed in two
# steps as we have more than one parameter:
# we extract the marginal posterior distribution for each element of the random effect
csi <- resultUH$marginals.random$Adj_ID[1:40]

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
png(filename=paste0("temp_", i,"_count.png") , width = 15.47 , height = 17.57 , units = "cm" , res=72)
spplot(obj=kilifi_sub, zcol= "cat.zeta",col.regions=diverge_hsv(8), scales=list(draw = TRUE), asp=1)
dev.off()
}



## interaction 2
ID.area.int <- admData2$Adj_ID
ID.year.int <- admData2$count_adm
temporalModel2 <- cumulitive_count ~  EVI_VALUE  + gender + 
  severe_disease + total_admission + rain_mm + admdays + nweight +
  f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
    hyper=list(prec.unstruct=list(prior="loggamma",param=c(0.001,0.001)),
               prec.spatial=list(prior="loggamma",param=c(0.1,0.01)))) +
  f( count_adm, model = "ar1") + f(ID.area.int,model="iid", group=ID.year.int,control.group=list(model="ar1")) 
  
 
result_tM2 <- inla(temporalModel2,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE),control.predictor(compute=TRUE) ,
                   control.inla = list(tolerance = 1e-20, h = 1e-08),E=log(nagem))

### 
delta.intII <- data.frame(delta=exp(result_tM2$summary.random$ID.area.int[,2]),
                              tempC=rep(1:11, each = 40) ,ID.area=result_tM2$summary.random$ID.area.int[,1])
delta.intII.matrix <- matrix(delta.intII[,1], 40,11,byrow=FALSE)
rownames(delta.intII.matrix)<- delta.intII[1:40,3]


# Check the absence of spatial trend for (intII)
 # cutoff.interaction <- c(-1,-0.01,0.01,1)
  cutoff.interaction <- c(0.20, 0.50,0.70, 0.999 ,1.0,1.01, 1.1,1.4,1.7)
  data.klf <- attr(kilifi_sub, "data")
    delta.intII.factor <- data.frame(NAME=data.klf$Adj_ID)
    for(i in 1:11){
      delta.factor.temp <- cut(delta.intII.matrix[,i],breaks=cutoff.interaction,include.lowest=TRUE) 
      delta.intII.factor <- cbind(delta.intII.factor,delta.factor.temp)
    }
    colnames(delta.intII.factor)<- c("NAME",seq(1,11))
    
# *** Code for Figure 7.5
attr(kilifi_sub, "data") <- data.frame(data.klf,  intII=delta.intII.factor)
trellis.par.set(axis.line=list(col=NA))
   
spplot(obj=kilifi_sub, zcol=c("intII.1","intII.2","intII.3",
                               "intII.4", "intII.5","intII.6",
                               "intII.7", "intII.8","intII.9",
                               "intII.10","intII.11"), 
       col.regions=diverge_hsv(8),
           names.attr=seq(1,11),main="")     
    

### interaction 3 - used this to report
ID.area.int <- admData2$Adj_ID
ID.year.int <- admData2$count_adm
temporalModel3 <- cumulitive_count ~  EVI_VALUE  + gender + 
  severe_disease +  rain_mm + admdays + nweight +
  f(Adj_ID, model = "bym"  ,graph=klf.adj , scale.model=TRUE,
    hyper=list(prec.unstruct=list(prior="loggamma",param=c(0.001,0.001)),
               prec.spatial=list(prior="loggamma",param=c(0.1,0.01)))) +
  f( count_adm, model = "ar1") + f(ID.year.int,model="iid", group=ID.area.int,
                                   control.group=list(model="besag",
                                                      graph=klf.adj))

result_tM3 <- inla(temporalModel3,family="nbinomial",
                   data=admData2, control.compute=list(dic=TRUE),control.predictor(compute=TRUE) ,
                   control.inla = list(tolerance = 1e-20, h = 1e-08),E=log(nagem_int))


delta.intIII <- data.frame(delta=exp(result_tM3$summary.random$ID.year.int[,2]),tempC=rep(1:11, each = 40),
                           ID.area=result_tM3$summary.random$ID.year.int[,1])
delta.intIII.matrix <- matrix(delta.intIII[,1], 40,11,byrow=FALSE)
rownames(delta.intIII.matrix)<- delta.intIII[1:40,3]

save.image("st_model3.RDA")

##load("st_model3.RDA")
cutoff.interaction <- c(0.50,0.8801.0, 1.3,1.9,2.5,3.4,7.4)
data.klf <- attr(kilifi_sub, "data")
delta.intIII.factor <- data.frame(NAME=data.klf$Adj_ID)
for(i in 1:11){
  delta.factor.temp <- cut(delta.intIII.matrix[,i],breaks=cutoff.interaction,include.lowest=TRUE ) 
  delta.intIII.factor <- cbind(delta.intIII.factor,delta.factor.temp)
}
colnames(delta.intIII.factor)<- c("NAME",seq(1,11))

# *** Code for Figure 7.6
attr(kilifi_sub, "data") <- data.frame(data.klf,  intIII=delta.intIII.factor)
trellis.par.set(axis.line=list(col=NA))

png(filename=paste0("temp_","img.png") , width = 25.47 , height = 27.57 , units = "cm" , res=300)

spplot(obj=kilifi_sub, zcol=c("intIII.1","intIII.2","intIII.3",
                              "intIII.4", "intIII.5","intIII.6",
                              "intIII.7", "intIII.8","intIII.9",
                              "intIII.10","intIII.11"), 
       col.regions=diverge_hsv(8),
       names.attr=seq(1,11),main="")    
dev.off()
