library("tidyverse")
library("randomForest")
library("caret")
library("skimr")
library("rsq")
library("ggplot2")
library("sp")
library("earlywarnings")
library("raster")
library("rgdal")
library("DFA")
library("foreach")
library("doParallel")
library("doSNOW")
library("trend")
library("corrplot")
library("relaimpo")
library("pdp")
library("RColorBrewer")
library("ggsci")
library("egg")
library("see")
library("epm")
library("patchwork")
library("ranger")
library("ggh4x")
library("scales")


######################################################################

#function to transform 3d data to 2d
trans3dto2d <- function(data3d){
  xl = nrow(data3d)
  yl = ncol(data3d)
  zl = length(data3d[1,1,])
  data2d = matrix(data = NA, nrow = xl*yl, ncol = zl, dimnames = NULL)
  for(x in 1:xl){
    for(y in 1:yl){
      loc = yl*(x-1)+y
      data2d[loc,] = as.vector(data3d[x,y,])
    }
  }
  return(data2d)
}

#function to transform 2d data to 3d
trans2dto3d <- function(data2d,xl,yl){
  if(xl*yl != length(data2d[,1])){
    return(NULL)
  }
  ln = colnames(data2d)
  zl = length(data2d[1,])
  data3d = array(data = NA, dim = c(xl,yl,zl), dimnames = NULL)
  for(x in 1:xl){
    for(y in 1:yl){
      loc = (x-1)*yl+y
      t = as.numeric(data2d[loc,])
      data3d[x,y,] = t
    }
  }
  return(list(data3d,ln))
}


list2num <- function(mlist){
  nro = nrow(mlist)
  nco = ncol(mlist)
  mnum = matrix(NA,nro,nco)
  for(i in 1:nro){
    for(j in 1:nco){
      if(is.numeric(mlist[[i,j]])){
        mnum[i,j] = mlist[[i,j]]
      }
    }
    
  }
  return(mnum)
}


wherestr <- function(allstr, keystrs, keyexclus = NULL){
  floc = 1:length(allstr)
  for (ks in keystrs) {
    loc = grep(ks, allstr)
    floc = intersect(loc,floc)
  }
  if(!is.null(keyexclus)){
    for (kes in keyexclus) {
      loc = grep(kes, allstr)
      floc = setdiff(floc,loc)
    }
  }
  
  return(floc)
}

wherestrs <- function(allstr,keystrs){
  rr = c()
  for (i in keystrs) {
    rr = c(rr,wherestr(allstr,i))
  }
  return(rr)
}




sen_mk <- function(tt){
  if(is.na(mean(tt))){
    result2 = c(NA,NA,NA)
  }else{
    result1 = sens.slope(tt)
    result2 = c(result1$estimates,result1$statistic,result1$p.value)
  }
  names(result2)<-c("SenSlope","z","p_value")
  return(result2)
}



remakeraster <- function(vecdata,locdata,refraster){
  if(length(vecdata)!=length(locdata)){
    return(print("Length of data is wrong!"))
  }
  
  mm = matrix(NA,nrow = nrow(refraster),ncol = ncol(refraster))
  mm[locdata] = vecdata
  rr = raster(mm)
  extent(rr) = extent(refraster)
  crs(rr) = crs(refraster)
  return(rr)
}


biostateloc <- function(biomemap, biometype, vegcovermap, thresholds){
  cond0 = (nrow(biomemap)==nrow(vegcovermap))&(ncol(biomemap)==ncol(vegcovermap))
  if(!cond0){
    return(print("Inconsistant input maps!"))
  }
  
  cond1 = which(as.matrix(biomemap==biometype)==T)
  cond2 = which(as.matrix((vegcovermap>=min(thresholds))&(vegcovermap<=max(thresholds)))==T)
  resultloc = intersect(cond1,cond2)
  return(resultloc)
}

standarizedf <- function(df00){
  df0 = df00
  df0 = as.matrix(df0)
  colnames(df0) = NULL
  maxs = apply(df0,2,max)
  mins = apply(df0,2,min)
  for (c in 1:ncol(df0)) {
    df0[,c] = (df0[,c]-mins[c])/(maxs[c]-mins[c])
  }
  df0 = as.data.frame(df0)
  colnames(df0) = colnames(df00)
  return(df0)
}



combinestrvec <- function(str1,sep = ""){
  r = ""
  for (i in 1:length(str1)) {
    if(i == 1){
      r = paste0(r,str1[i])
    }else{
      r = paste0(r,sep,str1[i])
    }
  }
  return(r)
}

delstr <- function(ostr, strs){
  if(length(strs)==0){
    return(ostr)
  }
  newstr = ostr
  for (st in strs) {
    newstr = gsub(st,"",newstr)
  }
  
  newstr = newstr[-which(newstr=="")]
  
  return(newstr)
}


rsq.ap <- function(actu,pred){
  meanactu = mean(actu,na.rm = T)
  rsq = 1-sum((actu-pred)^2)/sum((actu-meanactu)^2)
  return(rsq)
}



renamedr <- function(alldrivernames){
  newnames = c()
  for (nn in alldrivernames) {
    tailstr = ""
    if(substr(nn,1,1)=="t"){
      tailstr = "_t"
    }else if(substr(nn,1,1)=="m"){
      tailstr = "_m"
    }
    nn = sub("climsoil","climsoilwater",nn)
    nn = sub("clim","",nn)
    nn = sub("soil","",nn)
    nn = sub("fire","firefre",nn)
    nn = sub("hb","lsd",nn)
    nn = sub("tmean","temp",nn)
    nn = substring(nn,2)
    nn = toupper(nn)
    nn = paste0(nn,tailstr)
    newnames = c(newnames,nn)
  }
  return(newnames)
}

addunit <- function(strs){
  nounit = c("PDSI_t","PPT_t","SI_t","SRAD_t","TEMP_t","CO2_t","FIREFRE_t",
             "PDSI_m","PPT_m","SI_m","SRAD_m","TEMP_m","CO2_m","FIREFRE_m",
             "LSD_m","HUMANFP_m","CEC_m","SOC_m")
  drunits = c("(/yr)","(mm/yr)","(/yr)","(W/m^2/yr)","(℃/yr)","(ppm/yr)","(/yr)",
              "","(mm)","","(W/m^2)","(℃)","(ppm)","",
              "(TLU/km^2)","","(cmol/kg)","(g/kg)")
  
  
  newstr = c()
  for (i in 1:length(strs)) {
    str0 = strs[i]
    unit0 = drunits[which(nounit == str0)]
    newstr0 = paste0(str0," ",unit0)
    newstr = c(newstr,newstr0)
  }
  
  return(newstr)
}

mse.RF <- function(rfmodel,data0){
  var.name = rfmodel$dependent.variable.name
  predictions0 <- predict(rfmodel, data = data0)$predictions
  mse0 <- mean((data.test[[var.name]] - predictions0)^2)
  return(mse0)
}


########################################revised anadata



setwd(paste0("E:/AfricaResearch/StabilityTrend/data/AlignProcess/Pool_Aligned_data"))

map_afbs = raster("AFABSs.tif")

mapar1trend = raster("AF_istsg_AR1t_senslope_6yr.tif")
mapar1mean = raster("AF_istsg_AR1m.tif")


maptclimpdsi = raster("ClimateTrend_PDSI_SenSlope_2001_2020.tif")
maptclimppt = raster("ClimateTrend_ppt_SenSlope_2001_2020.tif")
maptclimsrad = raster("ClimateTrend_srad_SenSlope_2001_2020.tif")
maptclimtmean = raster("ClimateTrend_tmean_SenSlope_2001_2020.tif")
maptclimsi = raster("TerraClimate_SI_20ytrend_senslope.tif")

# maptclimtmax = raster("ClimateTrend_tmax_SenSlope_2001_2020.tif")
# maptclimtmin = raster("ClimateTrend_tmin_SenSlope_2001_2020.tif")
# maptclimai = raster("ClimateTrend_ai_SenSlope_2001_2020.tif")
# maptclimdef = raster("ClimateTrend_def_SenSlope_2001_2020.tif")
# maptclimpet = raster("ClimateTrend_pet_SenSlope_2001_2020.tif")
# maptclimsoil = raster("ClimateTrend_soil_SenSlope_2001_2020.tif")
# maptclimaet = raster("ClimateTrend_aet_SenSlope_2001_2020.tif")

# maptco2 = raster("CO2annualmean_20ytrend_senslope.tif")
# maptpop = raster("WB_Pop20ytrend_senslope.tif")
# maptgdp = raster("WB_GDP20ytrend_senslope.tif")
# maptfire = raster("FireFre_20ytrend_senslope.tif")
# maptndvi = raster("NDVITrend_SenSlope_2001_2020.tif")

mapmclimpdsi = raster("TerraClimate_20ym_PDSI.tif")
mapmclimppt = raster("TerraClimate_20ym_ppt.tif")
mapmclimsrad = raster("TerraClimate_20ym_srad.tif")
mapmclimtmean = raster("TerraClimate_20ym_tmean.tif")
mapmclimsi = raster("TerraClimate_SI_20ymean.tif")

# mapmclimtmax = raster("TerraClimate_20ym_tmax.tif")
# mapmclimtmin = raster("TerraClimate_20ym_tmin.tif")
# mapmclimai = raster("TerraClimate_20ym_ai.tif")
# mapmclimdef = raster("TerraClimate_20ym_def.tif")
# mapmclimpet = raster("TerraClimate_20ym_pet.tif")
# mapmclimsoil = raster("TerraClimate_20ym_soil.tif")
# mapmclimaet = raster("TerraClimate_20ym_aet.tif")

# mapmco2 = raster("CO2annualmean_20ymean.tif")
# mapmpop = raster("WB_Pop20ymean.tif")
# mapmgdp = raster("WB_GDP20ymean.tif")
# mapmfire = raster("FireFre_20ymean.tif")
# mapmndvi = raster("NDVI_20ymean.tif")

mapmhb = raster("AF_herbivoredens.tif")
mapmsoilcec = raster("iSDAsoil_cec.tif")
mapmsoilsoc = raster("iSDAsoil_soc.tif")
mapmhumanfp = raster("HFP2009.tif")


trendvars = ls()
trendvars = trendvars[wherestr(trendvars,"mapt")]

meanvars = ls()
meanvars = meanvars[wherestr(meanvars,"mapm")]

allvars = c(trendvars,meanvars)


biomefull = c('Closed Forest','Woody Savanna','Savanna',
              'Dense Grassland','Sparse Grassland',
              'Closed Shrubland','Open Shrubland')
biomeshort = c('CF','WS','SA','DG','SG','CS','OS')
lcnums = c(1,1,1,2,2,3,3)
ecmins = c(70,40,20,70,0,60,0)
ecmaxs = c(100,50,35,100,60,100,50)

typeinfo = data.frame(biomefull,biomeshort,lcnums,ecmins,ecmaxs)


datapath = "D:/Earthfuture_reanalysis/RF_ANA_data/"


# for(bsnum in c(1:7)){
#   fullname = biomefull[bsnum]
#   shortname = biomeshort[bsnum]
# 
#   bsloc = which(as.vector(map_afbs)==bsnum)
# 
#   AR1_trend = as.vector(mapar1trend)[bsloc]
#   AR1_mean = as.vector(mapar1mean)[bsloc]
# 
#   anadata = cbind(AR1_trend,AR1_mean)
# 
#   for(v in allvars){
#     ordertext1 = paste0("vdata = as.vector(", v,")[bsloc]")
#     class (ordertext1)
#     eval(parse(text = ordertext1))
#     anadata = cbind(anadata,vdata)
#   }
# 
#   colnames(anadata) = c("AR1_trend","AR1_mean",substring(allvars,4))
#   anadata = as.data.frame(anadata)
#   # ll = which(!is.na(apply(anadata, 1, sum)))
#   # anadata = anadata[ll,]
# 
#   saveRDS(anadata,file = paste0(datapath,"anadata_",
#                              shortname,".rds"),
#           compress=F)
# 
# }




####################################################################

setwd(paste0("D:/Earthfuture_reanalysis/RF_ANA_data"))

biomefull = c('Closed Forest','Woody Savanna','Savanna',
              'Dense Grassland','Sparse Grassland',
              'Closed Shrubland','Open Shrubland')
biomeshort = c('CF','WS','SA','DG','SG','CS','OS')

biometags = c('Closed\nForest','Woody\nSavanna','Savanna',
              'Dense\nGrassland','Sparse\nGrassland',
              'Closed\nShrubland','Open\nShrubland')

outpath.rf = "D:/Earthfuture_reanalysis/RF_ANA_data/result/"


datapath = "D:/Earthfuture_reanalysis/RF_ANA_data/"


################ 
#first loop for RF model training

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
              shortname,".rds"))
  
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  numoftree = 300
  sampsizerate = 0.6
  
  formu1 = as.formula(paste("AR1_trend ~", combinestrvec(alldrivernames,sep = " + ")))
  formu2 = as.formula(paste("AR1_mean ~", combinestrvec(meandrivernames,sep = " + ")))
  
  # set.seed(200+bsnum)
  # trains = createDataPartition(y = anadata1$AR1_trend,
  #                              p = 0.7,
  #                              list = F)
  # data.train = anadata1[trains,]
  # data.test = anadata1[-trains,]
  # saveRDS(data.train,file = paste0(datapath,"datatrain_",shortname,".rds"),
  #                    compress=F)
  # saveRDS(data.test,file = paste0(datapath,"datatest_",shortname,".rds"),
  #                    compress=F)
  
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  
  numoftree = 1000
  mn.size = 30
  sampsizerate = 0.6
  
  if(bsnum == 4){
    core.num = 8
  }else{
    core.num = 10
  }

  
  rf1 = ranger(
    formula = formu1,
    data = data.train,
    num.trees = numoftree,
    min.node.size = mn.size,
    mtry = floor(length(alldrivernames)/3),
    sample.fraction = sampsizerate,
    replace = F,
    num.threads = core.num,
    importance = "permutation",  
    probability = FALSE,  
    verbose = TRUE,
    seed = 200+bsnum              
  )
  
  print(Sys.time())
  
  rf2 = ranger(
    formula = formu2,
    data = data.train,
    num.trees = numoftree,
    min.node.size = mn.size,
    mtry = floor(length(meandrivernames)/3),
    sample.fraction = sampsizerate,
    replace = F,
    num.threads = core.num,
    importance = "permutation",  
    probability = FALSE,  
    verbose = TRUE,
    seed = 100+bsnum               
  )
  
  saveRDS(rf1,file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"),
                     compress=F)
  saveRDS(rf2,file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"),
                     compress=F)
  rm(data.test)
  rm(data.train)
  rm(rf1)
  rm(rf2)
  
  gc()
  print(paste0(fullname," is finished!"))
  print(Sys.time())
  
}



################ 
#second loop for partial relation 


for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
  rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
  
  
  rm(anadata)
  rm(anadata1)
  gc()
  
  
  resultframe1 = c()
  resultframe2 = c()
  
  
  nc = floor(as.numeric(10*(1024^3)/(object.size(rf1)+object.size(data.train)*1.5)))
  
  if(nc>10){
    nc = 10
  }
  
  print(paste0("Trend Model applying ",nc, " cores!"))
  
  cl <- makeCluster(nc,type = "SOCK")
  registerDoParallel(cl)
  
  for (xxx in 1:length(alldrivernames)) {
    
    partialrelation1 =  pdp::partial(rf1, pred.var = alldrivernames[xxx], grid.resolution = 20,
                                     type = "regression",parallel = T,
                                     paropts = list(.packages = "ranger"))
    
    resultframe01 = data.frame(xvar = partialrelation1[,1], 
                               yvar = partialrelation1[,2],
                               xname = rep(alldrivernames[xxx],length(partialrelation1[,1])),
                               yname = rep("AR1_trend",length(partialrelation1[,1])),
                               bsname = rep(shortname,length(partialrelation1[,1])))
  
    resultframe1 = rbind(resultframe1,resultframe01)
    
    print(paste0(xxx,"/",length(alldrivernames)))
  }
  
  stopCluster(cl)
  
  saveRDS(resultframe1,
          file = paste0(datapath,"Partialrela_",shortname,"_","AR1_trend.rds"),
          compress = F)
  
  rm(rf1)
  gc()
  
  #######
  
  nc = floor(as.numeric(10*(1024^3)/(object.size(rf2)+object.size(data.train)*1.5)))
  
  if(nc>10){
    nc = 10
  }
  
  print(paste0("Mean Model applying ",nc, " cores!"))
  
  cl <- makeCluster(nc,type = "SOCK")
  registerDoParallel(cl)
  
  for (xxx in 1:length(meandrivernames)) {
    
    partialrelation2 =  pdp::partial(rf2, pred.var = meandrivernames[xxx], grid.resolution = 20,
                                     type = "regression",parallel = T,
                                     paropts = list(.packages = "ranger"))
    

    resultframe02 = data.frame(xvar = partialrelation2[,1], 
                               yvar = partialrelation2[,2],
                               xname = rep(meandrivernames[xxx],length(partialrelation2[,1])),
                               yname = rep("AR1_mean",length(partialrelation2[,1])),
                               bsname = rep(shortname,length(partialrelation2[,1])))
    
    resultframe2 = rbind(resultframe2,resultframe02)
    
    print(paste0(xxx,"/",length(meandrivernames)))
  }
  
  stopCluster(cl)

  saveRDS(resultframe2,
          file = paste0(datapath,"Partialrela_",shortname,"_","AR1_mean.rds"),
          compress = F)
  


  rm(rf2)
  
  rm(resultframe1)
  rm(resultframe2)
  
  rm(data.test)
  rm(data.train)
  gc()
  
  print(paste0(fullname," is finished!"))
  print(Sys.time())
  
}



################ 
#third loop for plotting ntree

plotdata.ntree <- c()

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
  rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
  
  formu1 = as.formula(paste("AR1_trend ~", combinestrvec(alldrivernames,sep = " + ")))
  formu2 = as.formula(paste("AR1_mean ~", combinestrvec(meandrivernames,sep = " + ")))
  
  
  # ntreeseq = c(50,seq(100,600,100),800)
  ntreeseq = c(1000)
  mn.size = 30
  
  if(bsnum == 4){
    core.num = 8
  }else{
    core.num = 10
  }
  
  mse.seq.trend = c() 
  mse.seq.mean = c()
  rsqontrain.seq.trend = c()
  rsqontrain.seq.mean = c()
  
  for (numoftree in ntreeseq) {
    
    print(paste0("Number of tree = ", numoftree,"   ", Sys.time()))
    
    sampsizerate = 0.6
    
    xrf1 = ranger(
      formula = formu1,
      data = data.train,
      num.trees = numoftree,
      min.node.size = mn.size,
      mtry = floor(length(alldrivernames)/3),
      sample.fraction = sampsizerate,
      replace = F,
      num.threads = core.num,
      importance = "none",  
      probability = FALSE,  
      verbose = TRUE,
      seed = 200+bsnum              
    )
    
    mse.trend = mse.RF(xrf1,data.test)
    mse.seq.trend = c(mse.seq.trend,mse.trend)
    rsqontrain.seq.trend = c(rsqontrain.seq.trend,xrf1$r.squared)
    
    rm(xrf1)
    gc()
    
    xrf2 = ranger(
      formula = formu2,
      data = data.train,
      num.trees = numoftree,
      min.node.size = mn.size,
      mtry = floor(length(meandrivernames)/3),
      sample.fraction = sampsizerate,
      replace = F,
      num.threads = core.num,
      importance = "none",  
      probability = FALSE,  
      verbose = TRUE,
      seed = 100+bsnum               
    )
    
    mse.mean = mse.RF(xrf2,data.test)
    mse.seq.mean = c(mse.seq.mean,mse.mean)
    rsqontrain.seq.mean = c(rsqontrain.seq.mean,xrf2$r.squared)
    
    
    rm(xrf2)
    gc()
    
  }
  
  #mse.seq.trend = c(mse.seq.trend,mse.RF(rf1,data.test))
  #mse.seq.mean = c(mse.seq.mean,mse.RF(rf2,data.test))

  plotdata.ntree.error1 = data.frame(ntree = c(ntreeseq), mse = mse.seq.trend,
                                     rsqontest.seq= 1-(mse.seq.trend/var(data.test$AR1_trend)),
                                     rsqontrain.seq = rsqontrain.seq.trend,
                                     biome = rep(shortname,length(mse.seq.trend)),
                                     reactor = rep("AR1_trend",length(mse.seq.trend)))
  
  plotdata.ntree.error2 = data.frame(ntree = c(ntreeseq), mse = mse.seq.mean,
                                     rsqontest.seq= 1-(mse.seq.mean/var(data.test$AR1_mean)),
                                     rsqontrain.seq = rsqontrain.seq.mean,
                                     biome = rep(shortname,length(mse.seq.mean)),
                                     reactor = rep("AR1_mean",length(mse.seq.mean)))
  
  
  plotdata.ntree.save =  rbind(plotdata.ntree.error1,plotdata.ntree.error2)
  
  saveRDS(plotdata.ntree.save,
          file = paste0(datapath,"Ntree_vs_MSE_",shortname,"add.rds"),
          compress = F)
  
  plotdata.ntree = rbind(plotdata.ntree,plotdata.ntree.save)
  
  print(paste0(fullname," is finished!"))
  print(Sys.time())
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(rf1)
  rm(rf2)
  gc()
  
}




ntfiles = list.files(pattern = "Ntree_vs_MSE.*\\.rds$")

plotdata.ntree = c()

for(ff in ntfiles){
  dd0 = readRDS(ff)
  plotdata.ntree = rbind(plotdata.ntree,dd0)
}

for (ii in 1:length(biomeshort)) {
  plotdata.ntree$biome[plotdata.ntree$biome==biometags[ii]] = biomeshort[ii]
}


saveRDS(plotdata.ntree,
        file = paste0(datapath,"Ntree_vs_MSE_allBiomes.rds"),
        compress = F)

plotdata.ntree$biome = factor(plotdata.ntree$biome,
                              levels = biomeshort,
                              labels = biometags)

color.ntree = c("#3b0e80","#803a08")

plot.ntree <- ggplot(plotdata.ntree,
                     aes(x=ntree,y=mse,color = reactor,fill = reactor))+  
  geom_line( linewidth =0.6)+
  facet_wrap2(biome ~reactor,
              strip = strip_nested(),
              scales="free_y",
              strip.position = "right",
              nrow = 7,
              ncol = 2)+
  scale_y_continuous(n.breaks = 4, labels = scales::label_scientific())+
  scale_color_manual(values =  color.ntree)+
  xlab("Number of trees")+
  theme(panel.spacing.y = unit(0.4, "cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray30",linewidth = 0.2),
        strip.background = element_rect(fill="gray80"),
        legend.position="none")


ggsave(paste0(outpath.rf,"fig_ntree_mse.jpg"),plot.ntree,
       width = 15,height =18,units = "cm",dpi = 600)



################ 
#forth loop for plotting importance (new)

color.relaimpo = c("#ab3838","#2f3f8a","#c4a117")


# plotdata.relaimpo.all = c()
# 
# for(bsnum in c(1:7)){
#   
#   fullname = biomefull[bsnum]
#   shortname = biomeshort[bsnum]
#   
#   print(paste0(fullname," is under analysis!"))
#   print(Sys.time())
#   
#   anadata = readRDS(paste0(datapath,"anadata_",
#                            shortname,".rds"))
#   
#   anadata1 = na.omit(anadata)
#   
#   alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
#   meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
#   
#   #plotdrvnames = renamedr(originaldrvnames)
#   
#   
#   data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
#   data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
#   rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
#   rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
#   resultframe1 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_trend.rds"))
#   resultframe2 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_mean.rds"))
#   
#   
#   for(mtype in 1:2){
#     
#     if(mtype == 1){
#       rfm = rf1
#       partialrelation = resultframe1
#       yyname = "AR1_trend"
#       drivernames = alldrivernames
#     }else{
#       rfm = rf2
#       partialrelation = resultframe2
#       yyname = "AR1_mean"
#       drivernames = meandrivernames
#     }
#     
#     
#     data.pred.test = predict(rfm, data = data.test)
#     data.rsq = data.frame(actual = data.test[,mtype],predict = data.pred.test$predictions)
#     rsq.test = rsq.ap(data.rsq$actual, data.rsq$predict)
#     
#     
#     
#     statform.partrela = c()
#     
#     for (drv in drivernames) {
#       condition0 = (partialrelation$xname == drv)&
#         (partialrelation$yname == yyname)&
#         (partialrelation$bsname == shortname)
#       data.trend = partialrelation[condition0,]
#       sigtest.result = cor.test(x=data.trend$xvar,y=data.trend$yvar, 
#                                method ="pearson")
#       p.val = sigtest.result$p.value
#       tau = sigtest.result$estimate
#       
#       if(is.na(p.val*tau)){
#         rela = "non-monotone"
#       }else if(p.val<0.05&tau>0){
#         rela = "positive"
#       }else if(p.val<0.05&tau<0){
#         rela = "negative"
#       }else{
#         rela = "non-monotone"
#       }
#       
#       statform.partrela = rbind(statform.partrela,
#                                 data.frame(relation = rela, xvar = drv,
#                                            yvar = yyname, biome = shortname))
#       
#     }
#     
#     
#     
#     plotdata.relaimpo0 = data.frame(xvar = names(rfm$variable.importance),importance = rfm$variable.importance)
#     rownames(plotdata.relaimpo0) = NULL
#     
#     plotdata.relaimpo = merge(plotdata.relaimpo0,statform.partrela,by = "xvar")
#     plotdata.relaimpo = cbind(plotdata.relaimpo,
#                               data.frame(rsq.on.test = rep(rsq.test,nrow(plotdata.relaimpo))))
#     
#     plotdata.relaimpo$relation = factor(plotdata.relaimpo$relation,
#                                         levels = c("positive",
#                                                    "negative",
#                                                    "non-monotone"))
#     
#     
#     plotdata.relaimpo$xvar = renamedr(plotdata.relaimpo$xvar)
#     
#     
#     plotdata.relaimpo = cbind(plotdata.relaimpo,
#                               relative.importance = plotdata.relaimpo$importance/max(plotdata.relaimpo$importance))
#     
#     
#     plotdata.relaimpo.all = rbind(plotdata.relaimpo.all,plotdata.relaimpo)
#     
#   }
#   
#   rm(anadata)
#   rm(anadata1)
#   rm(data.test)
#   rm(data.train)
#   rm(resultframe1)
#   rm(resultframe2)
#   rm(rf1)
#   rm(rf2)
#   gc()
#   
# }
# 
# saveRDS(plotdata.relaimpo.all,file = paste0(datapath,"plotdata_relaimpo_all.rds"),
#         compress=F)


plotdata.relaimpo.all = readRDS(file = paste0(datapath,"plotdata_relaimpo_all.rds"))

plotdata.relaimpo.mean = plotdata.relaimpo.all[plotdata.relaimpo.all$yvar=="AR1_mean",]

driver_order <- plotdata.relaimpo.mean %>%
  group_by(xvar) %>%
  summarise(mean_importance = mean(relative.importance, na.rm = TRUE)) %>%
  arrange(desc(mean_importance)) %>%
  pull(xvar)

plotdata.relaimpo.mean$xvar <- factor(plotdata.relaimpo.mean$xvar, 
                                      levels = (driver_order),
                                      labels = gsub("HUMANFP_m", "HFP_m", driver_order))


abscolors = c("#267300","#a9a800","#ff5500","#00a884","#ffaa01","#004da7","#ff01c4")


pic.relaimport.mean = ggplot(data = plotdata.relaimpo.mean, 
                             mapping = aes( x = xvar,y = relative.importance,
                                            color = biome,
                                            group = biome))+
  geom_point(size = 3.5,alpha = 0.8,
             shape = 21,
             stroke = 2,
             position = position_jitter(height = 0, width = 0.0))+
  # geom_line(linewidth = 1)+
  # coord_flip()+
  xlab("")+
  ylab("Relative Importance")+
  labs(tag = "a")+
  scale_color_manual(values = abscolors, labels = biometags, name = NULL)+
  scale_x_discrete(expand = expansion(mult = c(0.5/9, 0.5/9)))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  guides(color = guide_legend(nrow = 1))+
  theme_classic()+
  theme(legend.position="top",
        # plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt"),
        # axis.ticks.x = element_blank(),       
        axis.text.x = element_text(color = 'black')
  )


#####

plotdata.partial.trend = c()
plotdata.partial.mean = c()

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
  rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
  resultframe1 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_trend.rds"))
  resultframe2 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_mean.rds"))
  
  plotdata.partial.trend = rbind(plotdata.partial.trend,resultframe1)
  plotdata.partial.mean = rbind(plotdata.partial.mean,resultframe2)
  
}



plotdata.partial.mean$xname <- renamedr(plotdata.partial.mean$xname)
plotdata.partial.mean$xname <- factor(plotdata.partial.mean$xname, 
                                      levels = (driver_order),
                                      labels = gsub("HUMANFP_m", "HFP_m", driver_order))
plotdata.partial.mean$bsname <- factor(plotdata.partial.mean$bsname, levels = biomeshort)


cor.relation = c()
for (ir in 1:nrow(plotdata.partial.mean)) {
  loc = which(plotdata.relaimpo.mean$xvar == plotdata.partial.mean$xname[ir]&
                plotdata.relaimpo.mean$yvar == plotdata.partial.mean$yname[ir]&
                plotdata.relaimpo.mean$biome  == plotdata.partial.mean$bsname[ir])
  cor.relation = c(cor.relation,as.character(plotdata.relaimpo.mean$relation[loc]))
}


plotdata.partial.mean = cbind(plotdata.partial.mean,cor.relation)
plotdata.partial.mean$cor.relation = factor(plotdata.partial.mean$cor.relation,
                                            levels = c("positive",
                                                       "negative",
                                                       "non-monotone"),
                                            labels = c("Positive",
                                                       "Negative",
                                                       "Non-significant"))


plotdata.partial.mean$bsname = factor(plotdata.partial.mean$bsname,
                                      levels = biomeshort,
                                      labels = biometags)


pic.partial.mean = ggplot(data = plotdata.partial.mean, 
                          mapping = aes(x = xvar, y = yvar, color = cor.relation))+
  geom_line(linewidth = 1)+
  facet_grid(bsname ~ xname, scales = "free",
             labeller = labeller(bsname = biometags))+
  scale_color_manual(values = color.relaimpo, name = NULL)+
  labs(tag = "b")+
  xlab("")+
  ylab("AR1 mean")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = NULL,
    axis.text.x = element_text(angle = 60,hjust = 1),
    plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt"),
    legend.position = "bottom"
  )



patch.combine.mean = pic.relaimport.mean/pic.partial.mean + plot_layout(heights = c(2,7))




outfile_patch.mean = paste0(outpath.rf,"fig_relaimport_mean.jpg")




jpeg(outfile_patch.mean, units="cm",
     width= 20,
     height= 24,
     res = 800)

plot.new()

plot(patch.combine.mean)

dev.off()


#####



plotdata.relaimpo.trend = plotdata.relaimpo.all[plotdata.relaimpo.all$yvar=="AR1_trend",]

driver_order <- plotdata.relaimpo.trend %>%
  group_by(xvar) %>%
  summarise(trend_importance = mean(relative.importance, na.rm = TRUE)) %>%
  arrange(desc(trend_importance)) %>%
  pull(xvar)

driver_order.grouped = c(driver_order[wherestr(driver_order,"_m")],
                         driver_order[wherestr(driver_order,"_t")])


plotdata.relaimpo.trend$xvar <- factor(plotdata.relaimpo.trend$xvar, 
                                       levels = (driver_order.grouped),
                                       labels = gsub("HUMANFP_m", "HFP_m", driver_order.grouped))


abscolors = c("#267300","#a9a800","#ff5500","#00a884","#ffaa01","#004da7","#ff01c4")

dash.x = length(wherestr(driver_order,"_m"))+0.5

pic.relaimport.trend = ggplot(data = plotdata.relaimpo.trend, 
                              mapping = aes( x = xvar,y = relative.importance,
                                             color = biome,
                                             group = biome))+
  geom_point(size = 3.5,alpha = 0.95,
             shape = 21,
             stroke = 2,
             position = position_jitter(height = 0, width = 0.0))+
  geom_vline(xintercept = dash.x, linetype = 'dashed')+
  # geom_line(linewidth = 1)+
  # coord_flip()+
  xlab("")+
  ylab("Relative Importance")+
  labs(tag = "a")+
  scale_color_manual(values = abscolors, labels = biometags, name = NULL)+
  scale_x_discrete(expand = expansion(mult = c(0.5/14, 0.5/14)))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  guides(color = guide_legend(nrow = 1))+
  theme_classic()+
  theme(legend.position="top",
        # plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt"),
        # axis.ticks.x = element_blank(),       
        axis.text.x = element_text(color = 'black')
  )


#####

plotdata.partial.trend = c()
plotdata.partial.mean = c()

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
  rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
  resultframe1 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_trend.rds"))
  resultframe2 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_mean.rds"))
  
  plotdata.partial.trend = rbind(plotdata.partial.trend,resultframe1)
  plotdata.partial.mean = rbind(plotdata.partial.mean,resultframe2)
  
}



plotdata.partial.trend$xname <- renamedr(plotdata.partial.trend$xname)
plotdata.partial.trend$xname <- factor(plotdata.partial.trend$xname, 
                                       levels = (driver_order.grouped),
                                       labels = gsub("HUMANFP_m", "HFP_m", driver_order.grouped))
plotdata.partial.trend$bsname <- factor(plotdata.partial.trend$bsname, levels = biomeshort)


cor.relation = c()
for (ir in 1:nrow(plotdata.partial.trend)) {
  loc = which(plotdata.relaimpo.trend$xvar == plotdata.partial.trend$xname[ir]&
                plotdata.relaimpo.trend$yvar == plotdata.partial.trend$yname[ir]&
                plotdata.relaimpo.trend$biome  == plotdata.partial.trend$bsname[ir])
  cor.relation = c(cor.relation,as.character(plotdata.relaimpo.trend$relation[loc]))
}


plotdata.partial.trend = cbind(plotdata.partial.trend,cor.relation)
plotdata.partial.trend$cor.relation = factor(plotdata.partial.trend$cor.relation,
                                             levels = c("positive",
                                                        "negative",
                                                        "non-monotone"),
                                             labels = c("Positive",
                                                        "Negative",
                                                        "Non-significant"))


plotdata.partial.trend$bsname = factor(plotdata.partial.trend$bsname,
                                       levels = biomeshort,
                                       labels = biometags)


pic.partial.trend = ggplot(data = plotdata.partial.trend, 
                           mapping = aes(x = xvar, y = yvar, color = cor.relation))+
  geom_line(linewidth = 1)+
  facet_grid(bsname ~ xname, scales = "free",
             labeller = labeller(bsname = biometags))+
  scale_color_manual(values = color.relaimpo, name = NULL)+
  scale_y_continuous(labels = label_scientific())+
  labs(tag = "b")+
  xlab("")+
  ylab("AR1 trend")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = NULL,
    axis.text.x = element_text(angle = 60,hjust = 1),
    plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt"),
    legend.position = "bottom"
  )



patch.combine.trend = pic.relaimport.trend/pic.partial.trend + plot_layout(heights = c(2,7))




outfile_patch.trend = paste0(outpath.rf,"fig_relaimport_trend.jpg")




jpeg(outfile_patch.trend, units="cm",
     width= 30,
     height= 24,
     res = 800)

plot.new()

plot(patch.combine.trend)

dev.off()




################ 
#seventh loop for plotting cross validation 


color.cv = c("#CD7054","#6CA6CD")



list_trend = list()
list_mean = list()

for(bsnum in c(1:7)){
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  data.train = readRDS(file = paste0(datapath,"datatrain_",shortname,".rds"))
  data.test = readRDS(file = paste0(datapath,"datatest_",shortname,".rds"))
  rf1 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_trend.rds"))
  rf2 = readRDS(file = paste0(datapath,"RFmodel_",shortname,"_","AR1_mean.rds"))
  resultframe1 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_trend.rds"))
  resultframe2 = readRDS(file = paste0(datapath,"Partialrela_",shortname,"_","AR1_mean.rds"))


  
  
  for(mtype in 1:2){
    
    if(mtype == 1){
      rfm = rf1
      partialrelation = resultframe1
      yyname = "AR1_trend"
    }else{
      rfm = rf2
      partialrelation = resultframe2
      yyname = "AR1_mean"
    }
    
    
    
    
    data.pred.test = predict(rfm, data = data.test)
    data.pred.train = predict(rfm, data = data.train)
    data.crossverifiy = data.frame(actual = c(data.test[,mtype],data.train[,mtype]),
                                   predict = c(data.pred.test$predictions,data.pred.train$predictions),
                                   subset = c(rep("test",length(data.pred.test$predictions)),
                                              rep("train",length(data.pred.train$predictions))))
    value.max = max(c(data.crossverifiy$actual,data.crossverifiy$predict))
    value.min = min(c(data.crossverifiy$actual,data.crossverifiy$predict))
    
    plot.range = c(value.min + (value.max-value.min)*0,
                   value.max - (value.max-value.min)*0)
    
    rsq.on.test = rsq.ap(data.test[,mtype],data.pred.test$predictions)
    rsq.on.train = rsq.ap(data.train[,mtype],data.pred.train$predictions)
    
    
    
    text.x = plot.range[2]-(plot.range[2]-plot.range[1])*0.21
    text.y1 = plot.range[1]+(plot.range[2]-plot.range[1])*0.09  
    text.y2 = plot.range[1]+(plot.range[2]-plot.range[1])*0.29 
    
    
    data.crossverifiy$subset = factor(data.crossverifiy$subset,
                                      levels = c("test","train"))
    
    
    plot.cv = ggplot(data.crossverifiy, aes(x = predict, y = actual, group = subset)) + 
      stat_density_2d(geom = "polygon",aes(alpha = after_stat(level), fill = subset))+
      geom_abline(slope=1, intercept=0, 
                  linewidth = 0.5, linetype = 2,
                  color = "gray20")+
      xlim(plot.range)+
      ylim(plot.range)+
      xlab("predict")+
      ylab("observation")+
      scale_color_manual(values = color.cv)+
      scale_fill_manual(values =  color.cv)+
      coord_fixed(ratio = 1)+
      labs(tag =  paste0("  ",letters[abs(mtype-2)*7+bsnum]),
           subtitle = paste0(yyname,"\n",fullname))+
      annotate(geom = "label", x = text.x, y = text.y1, 
               fill = alpha(color.cv[1], 0.6), 
               label = sprintf("italic(R)[ test ]^2 == '%0.2f'",rsq.on.test),
               size = 3,parse = T)+
      annotate(geom = "label", x = text.x, y = text.y2, 
               fill = alpha(color.cv[2], 0.6), 
               label = sprintf("italic(R)[train]^2 == '%0.2f'",rsq.on.train),
               size = 3,parse = T)+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.position = "none",
            plot.tag = element_text(face = "bold"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
    
    filename.pic.cv = paste0(outpath.rf,"cv_pic/figcv_",
                             shortname,"_",yyname,".jpg")
    
    
    ggsave(filename.pic.cv,plot.cv,
           width = 6.5,height = 6.5,units = "cm",dpi = 600)
    
    if(mtype == 1){
      list_trend = c(list_trend,list(plot.cv))
    }else{
      list_mean = c(list_mean,list(plot.cv))
    }
    
    
    
  }
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(rf1)
  rm(rf2)
  gc()
  
}


patch.combine.cv = list_mean[[1]]+list_mean[[2]]+list_mean[[3]]+list_mean[[4]]+
  list_mean[[5]]+list_mean[[6]]+list_mean[[7]]+plot_spacer()+
  list_trend[[1]]+list_trend[[2]]+list_trend[[3]]+list_trend[[4]]+
  list_trend[[5]]+list_trend[[6]]+list_trend[[7]]+plot_spacer()+
  plot_layout(ncol = 4)



outfile_patch = paste0(outpath.rf,"figcv_1combined_new.jpg")
jpeg(outfile_patch, units="cm",
     width= 6.25*4,
     height= 6.25*4,
     res = 600)

plot.new()

plot(patch.combine.cv)

dev.off()
















################
#correlation between drivers

data.driver = c()

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  anadata = readRDS(paste0(datapath,"anadata_",
                           shortname,".rds"))
  
  anadata1 = na.omit(anadata)
  
  alldrivernames = colnames(anadata1)[3:ncol(anadata1)]
  meandrivernames = colnames(anadata1)[8:ncol(anadata1)]
  
  
  needcols = wherestrs(colnames(anadata1),alldrivernames)
  
  data00 = anadata1[,needcols] 
  data.driver = rbind(data.driver,data00)
  
  rm(anadata)
  rm(anadata1)
  rm(data00)
  
  gc()
  print(paste0(fullname," is finished!"))
  print(Sys.time())
}

colnames(data.driver) = renamedr(colnames(data.driver))

corresult = cor(data.driver,method = "pearson")

mycolors1 = c("#CD7054","#FDF5E6","#6CA6CD")
mycolorbar <- colorRampPalette(mycolors1,space = "Lab")


outfile = paste0(outpath.rf,"corroplot_drivers.jpg")
jpeg(outfile, units="cm", width=16, height=16, res=600)
plot.new()


corrplot(corresult, 
         method = "square", type = "lower",
         order = 'FPC',
         addCoef.col = "black",
         number.cex = 0.63,
         col = mycolorbar(200), 
         tl.pos = "ld",
         tl.col = "black",
         tl.cex = 0.7, 
         cl.pos = 'b')

dev.off()



