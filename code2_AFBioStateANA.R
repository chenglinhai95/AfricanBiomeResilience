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




####################################################################

setwd(paste0("E:/AfricaResearch/StabilityTrend/data/AlignProcess/Pool_Aligned_data"))



biomefull = c('Closed Forest','Woody Savanna','Savanna',
              'Dense Grassland','Sparse Grassland',
              'Closed Shrubland','Open Shrubland')
biomeshort = c('CF','WS','SA','DG','SG','CS','OS')

biometags = c('Closed\nForest','Woody\nSavanna','Savanna',
              'Dense\nGrassland','Sparse\nGrassland',
              'Closed\nShrubland','Open\nShrubland')



outpath1 = "E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/corrplot/"
outpath2 = "E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/lmgresultforms/"
outpath3 = "E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/mean_trend/"
outpath.rf = "E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/RFana_noNDVI/"


################ 
#first loop for RF model training

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  numoftree = 300
  sampsizerate = 0.1
  featurenum = round(0.5*length(alldrivernames))
  
  formu1 = as.formula(paste("AR1_trend ~", combinestrvec(alldrivernames,sep = " + ")))
  formu2 = as.formula(paste("AR1_mean ~", combinestrvec(alldrivernames,sep = " + ")))
  
  set.seed(200+bsnum)
  trains = createDataPartition(y = anadata1$AR1_trend,
                               p = 0.7,
                               list = F)
  
  data.train = anadata1[trains,]
  data.test = anadata1[-trains,]
  
  
  set.seed(300+bsnum)
  rf1 = randomForest(formu1,
                     data = data.train,
                     ntree = numoftree,
                     mtry = featurenum,
                     sampsize = round(sampsizerate*nrow(data.train)),
                     importance = T)
  
  
  set.seed(300-bsnum)
  rf2 = randomForest(formu2,
                     data = data.train,
                     ntree = numoftree,
                     mtry = featurenum,
                     sampsize = round(sampsizerate*nrow(data.train)),
                     importance = T)
  
  save(data.train,file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  save(data.test,file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  save(rf1,file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  save(rf2,file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  
  
  
  
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
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  
  
  resultframe1 = c()
  resultframe2 = c()
  
  for (xxx in 1:length(alldrivernames)) {
    
    nc = 18
    cl <- makeCluster(nc,type = "SOCK")
    registerDoParallel(cl)
    
    partialrelation1 =  pdp::partial(rf1, pred.var = alldrivernames[xxx], grid.resolution = 30,
                                     type = "regression",parallel = T,
                                     paropts = list(.packages = "randomForest"))
    
    
    partialrelation2 =  pdp::partial(rf2, pred.var = alldrivernames[xxx], grid.resolution = 30,
                                     type = "regression",parallel = T,
                                     paropts = list(.packages = "randomForest"))
    stopCluster(cl)
    
    resultframe01 = data.frame(xvar = partialrelation1[,1], 
                               yvar = partialrelation1[,2],
                               xname = rep(alldrivernames[xxx],length(partialrelation1[,1])),
                               yname = rep("AR1_trend",length(partialrelation1[,1])),
                               bsname = rep(shortname,length(partialrelation1[,1])))
    
    resultframe02 = data.frame(xvar = partialrelation2[,1], 
                               yvar = partialrelation2[,2],
                               xname = rep(alldrivernames[xxx],length(partialrelation2[,1])),
                               yname = rep("AR1_mean",length(partialrelation2[,1])),
                               bsname = rep(shortname,length(partialrelation2[,1])))
    
    resultframe1 = rbind(resultframe1,resultframe01)
    resultframe2 = rbind(resultframe2,resultframe02)
    
    print(paste0(xxx,"/",length(alldrivernames)))
  }
  
  save(resultframe1,file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  save(resultframe2,file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(rf1)
  rm(rf2)
  rm(resultframe1)
  rm(resultframe2)
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
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  
  
  
  plotdata.ntree.error1 = data.frame(ntree = c(1:rf1$ntree), mse = rf1$mse,
                                     biome = rep(shortname,rf1$ntree),
                                     reactor = rep("AR1_trend",rf1$ntree))
  plotdata.ntree.error2 = data.frame(ntree = c(1:rf2$ntree), mse = rf2$mse,
                                     biome = rep(shortname,rf2$ntree),
                                     reactor = rep("AR1_mean",rf2$ntree))
  
  
  plotdata.ntree = rbind(plotdata.ntree,plotdata.ntree.error1,plotdata.ntree.error2)
  
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(rf1)
  rm(rf2)
  gc()
  
}


plotdata.ntree$biome = factor(plotdata.ntree$biome,
                              levels = biomeshort,
                              labels = biometags)


color.ntree = c("#3b0e80","#803a08")

plot.ntree <- ggplot(plotdata.ntree,
                     aes(x=ntree,y=mse,color = reactor,fill = reactor))+  
  geom_line( linewidth =0.6)+
  facet_grid(vars(reactor), vars(biome), scales="free")+
  scale_y_continuous(n.breaks = 4)+
  scale_color_manual(values =  color.ntree)+
  xlab("Number of trees")+
  theme(panel.spacing.y = unit(0.4, "cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray30",linewidth = 0.2),
        strip.background = element_rect(fill="gray80"),
        legend.position="none")


ggsave(paste0(outpath.rf,"fig_ntree_mse.jpg"),plot.ntree,
       width = 25,height = 8,units = "cm",dpi = 600)



################ 
#forth loop for plotting importance


color.relaimpo = c("#ab3838","#2f3f8a","#c4a117")


list_trend = list()
list_mean = list()


for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
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
    
    
    data.pred.test = predict(rfm, newdata = data.test)
    data.rsq = data.frame(actual = data.test[,mtype],predict = data.pred.test)
    rsq.test = rsq.ap(data.rsq$actual, data.rsq$predict)
    
    
    
    statform.partrela = c()
    
    for (drv in alldrivernames) {
      condition0 = (partialrelation$xname == drv)&
        (partialrelation$yname == yyname)&
        (partialrelation$bsname == shortname)
      data.trend = partialrelation[condition0,]
      mktest.result = cor.test(x=data.trend$xvar,y=data.trend$yvar, 
                               method ="pearson")
      p.val = mktest.result$p.value
      tau = mktest.result$estimate
      
      if(p.val<0.05&tau>0){
        rela = "positive"
      }else if(p.val<0.05&tau<0){
        rela = "negative"
      }else{
        rela = "non-monotone"
      }
      
      statform.partrela = rbind(statform.partrela,
                                data.frame(relation = rela, xvar = drv,
                                           yvar = yyname, biome = shortname))
      
    }
    
    
    
    plotdata.relaimpo0 = data.frame(xvar = rownames(rfm$importance),importance = rfm$importance[,1])
    rownames(plotdata.relaimpo0) = NULL
    
    plotdata.relaimpo = merge(plotdata.relaimpo0,statform.partrela,by = "xvar")
    plotdata.relaimpo = cbind(plotdata.relaimpo,
                              data.frame(rsq.on.test = rep(rsq.test,nrow(plotdata.relaimpo))))
    
    plotdata.relaimpo$relation = factor(plotdata.relaimpo$relation,
                                        levels = c("positive",
                                                   "negative",
                                                   "non-monotone"))
    
    
    plotdata.relaimpo$xvar = renamedr(plotdata.relaimpo$xvar)
    
    yrng <- range(plotdata.relaimpo$importance)
    textloc.x = 1.2
    textloc.y = yrng[2]*0.6
    
    
    tick.loc = signif(yrng[2]*0.6,1)
    
    text.title = paste0( fullname,"\n",yyname)
    
    
    plot.relaimpo = ggplot(plotdata.relaimpo,
                           aes(x=reorder(xvar, importance),y=importance,
                               color = relation, fill = relation)) + 
      geom_bar(stat = 'identity', width = 0.5) +
      coord_flip() +
      scale_color_manual(values =  color.relaimpo)+
      scale_fill_manual(values =  color.relaimpo)+
      scale_y_continuous(breaks = c(0,tick.loc))+
      xlab(NULL)+
      labs(tag = paste0("   ",letters[abs(mtype-2)*7+bsnum]))+
      ggtitle(text.title)+
      annotate(geom = "text", x = textloc.x, y = textloc.y, 
               label = sprintf("italic(R)^2 == '%0.2f'",rsq.test),
               size = 5,parse = T)+
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(),
            axis.text.y = element_text(color = "black",size = 14),
            axis.title.x = element_text(color = "black",size = 14),
            axis.text.x = element_text(color = "gray20",size = 12),
            plot.title = element_text(hjust = -0.2,
                                      size = 13),
            plot.tag.position = c(0,1),
            plot.tag = element_text(face = "bold",
                                    hjust = -0.4,
                                    size = 20),
            legend.position="none",
            legend.title = element_blank(),
            legend.text = element_text(color = "black",size = 14))
    
    
    filename.pic.relaimpo = paste0(outpath.rf,"relaimpo_pic/figrelaimpo_",
                                   shortname,"_",yyname,".jpg")
    
    
    ggsave(filename.pic.relaimpo,plot.relaimpo,
           width = 6,height = 15,units = "cm",dpi = 600)
    
    
    if(mtype ==1 & bsnum==1){
      plot.relaimpo = plot.relaimpo + theme(legend.position="bottom",)
    }
    
    
    if(mtype == 1){
      list_trend = c(list_trend,list(plot.relaimpo))
    }else{
      list_mean = c(list_mean,list(plot.relaimpo))
    }
    
    
    
  }
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(resultframe1)
  rm(resultframe2)
  rm(rf1)
  rm(rf2)
  gc()
  
}

library(patchwork)
patch.trend = c()
patch.mean = c()

for (xx in 1:7) {
  patch.trend <- patch.trend | list_trend[[xx]]
  patch.mean <- patch.mean | list_mean[[xx]]
}

patch.combine <-  (patch.mean + 
                     plot_layout(guides = 'collect') & 
                     theme(legend.position = 'none'))/(patch.trend +
                                                         plot_layout(guides = 'collect')& theme(legend.position = 'bottom'))

outfile_patch = paste0(outpath.rf,"relaimpo_pic/figrelaimpo_1combined.jpg")


jpeg(outfile_patch, units="cm",
     width= 6.5*7,
     height= 14*2.2,
     res = 600)

plot.new()

plot(patch.combine)

dev.off()



################ 
#fifth loop for plotting partial relation


colpalettes<-unique(c(pal_simpsons("springfield")(16),pal_gsea("default")(12),
                      pal_startrek("uniform")(7),pal_lancet("lanonc")(9)))

set.seed(222)
color.xvars = sample(colpalettes,length(alldrivernames))

color.table = data.frame(xvar = alldrivernames,
                         plotname = renamedr(alldrivernames),
                         color = color.xvars,
                         od = 1:length(alldrivernames))

color.table$color[5] = "#ABA80EFF"
color.table$color[9] = "#5F08A1FF"


cplotdata.t5 = c()
list_fig = vector(mode='list', 14)

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
  for(mtype in 1:2){
    
    if(mtype == 1){
      rfm = rf1
      partialrelation = resultframe1
      yyname = "AR1_trend"
      plotorder = 7+bsnum
      tagnote = letters[plotorder]

    }else{
      rfm = rf2
      partialrelation = resultframe2
      yyname = "AR1_mean"
      plotorder = bsnum
      tagnote = letters[plotorder]
    }
    
    
    plotdata.relaimpo0 = data.frame(xvar = rownames(rfm$importance),importance = rfm$importance[,1])
    rownames(plotdata.relaimpo0) = NULL
    imporank = rank(-plotdata.relaimpo0$importance)
    plotdata.relaimpo = cbind(plotdata.relaimpo0,imporank)
    
    
    plotdata.partrela = c()
    
    for (drv in alldrivernames) {
      imporank = plotdata.relaimpo$imporank[plotdata.relaimpo$xvar==drv]
      plotdata.relaimpo0$xvar==drv
      condition0 = (partialrelation$xname == drv)&
        (partialrelation$yname == yyname)&
        (partialrelation$bsname == shortname)
      subdata = partialrelation[condition0,]
      
      xvar.st = subdata$xvar
      
      subdata.st = standarizedf(xvar.st)
      
      subdata[,1] = subdata.st
      imporank = rep(imporank,nrow(subdata))
      subdata = cbind(subdata, imporank)
      
      
      plotdata.partrela = rbind(plotdata.partrela,subdata)
      
    }
    
    
    plotdata.partrela = as.data.frame(plotdata.partrela)
    
    
    plotdata.partrela.t5 = plotdata.partrela[(plotdata.partrela$imporank<=5),]
    cplotdata.t5 = rbind(cplotdata.t5,plotdata.partrela.t5)
    
    xvar.temp = unique(plotdata.partrela.t5$xname)
    color.table.loc = wherestrs(color.table$xvar,xvar.temp)
    color.table.temp = color.table[color.table.loc,]
    
    
    plotdata.partrela.t5$xname = renamedr(plotdata.partrela.t5$xname)
    plotdata.partrela.t5$xname = factor(plotdata.partrela.t5$xname,
                                        levels = color.table.temp$plotname)
    
    
    
    plot.partrela.t5 = ggplot(plotdata.partrela.t5,
                              aes(x=xvar,y=yvar,
                                  color = xname, fill = xname)) + 
      geom_line(linewidth=0.8)+
      labs(tag = tagnote,
           title = fullname)+
      xlab("Influencing Factor")+
      ylab(yyname)+
      scale_color_manual(values = color.table.temp$color)+
      theme(panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),
            axis.line = element_line(),
            axis.text.y = element_text(color = "black"),
            plot.tag = element_text(face = "bold"),
            plot.title = element_text(face = "plain",size=11),
            legend.position='none',
            legend.title = element_blank())
    
    list_fig[[plotorder]] = plot.partrela.t5
    
    filename.pic.partial = paste0(outpath.rf,"partial_pic/figpartial_t5_",
                                  shortname,"_",yyname,".jpg")
    
    
    ggsave(filename.pic.partial,plot.partrela.t5,
           width = 9,height = 6,units = "cm",dpi = 600)
    
    
  }
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(resultframe1)
  rm(resultframe2)
  rm(rf1)
  rm(rf2)
  gc()
  
}


xvar.temp = unique(cplotdata.t5$xname)
color.table.loc = wherestrs(color.table$xvar,xvar.temp)
color.table.temp = color.table[sort(color.table.loc),]

color.table.temp$plotname = factor(color.table.temp$plotname,
                             levels = color.table.temp$plotname)

x = c(1,2,3,1,2,3,1,2,3)
ykey = c(6.2,6.2,6.2,4.2,4.2,4.2,2.2,2.2,2.2)
ylab = c(5.5,5.5,5.5,3.5,3.5,3.5,1.5,1.5,1.5)

color.table.legend = cbind(color.table.temp,x,ykey,ylab)

layout <- "
ABCD
EFGH
IJKL
MNOO
"


plot.leg = ggplot()+
  geom_point(data = color.table.legend, aes(x = x, y = ykey, color = plotname),
             shape = '-',size = 13)+
  scale_color_manual(values = color.table.legend$color)+
  xlim(0.5,3.5)+
  ylim(1,6.25)+
  geom_text(data = color.table.legend, aes(x = x, y = ylab, label = plotname))+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        legend.position='none')
  
  


plotgroup = list_fig[[1]]+list_fig[[2]]+list_fig[[3]]+
  list_fig[[4]]+list_fig[[5]]+list_fig[[6]]+
  list_fig[[7]]+list_fig[[8]]+list_fig[[9]]+
  list_fig[[10]]+list_fig[[11]]+list_fig[[12]]+
  list_fig[[13]]+list_fig[[14]]+plot.leg+
  plot_layout(design = layout)



outfile_pn = paste0(outpath.rf,"partial_pic/figpartial_1_t5combinednew.jpg")
jpeg(outfile_pn, units="cm",
     width= 22,
     height= 20,
     res = 600)

plot.new()

plot(plotgroup)

dev.off()





cplotdata.t5$xname = renamedr(cplotdata.t5$xname)
cplotdata.t5$xname = factor(cplotdata.t5$xname,
                            levels = color.table.temp$plotname)
cplotdata.t5$bsname = factor(cplotdata.t5$bsname,
                             levels = biomeshort,
                             labels = biometags)




cplot.partrela.t5 = ggplot(cplotdata.t5,
                           aes(x=xvar,y=yvar,
                               color = xname, fill = xname)) + 
  facet_grid(rows = vars(yname),cols = vars(bsname), scales = "free")+
  geom_line(size=0.9)+  
  scale_color_manual(values = color.table.temp$color)+
  scale_x_continuous(breaks = c(0.2,0.5,0.8))+
  xlab("Influencing Factor")+
  ylab(NULL)+
  guides(fill = guide_legend(linetype = 1, 
                             linewidth = 2,
                             ncol=2))+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(size = 11),
        legend.position="right",
        legend.title = element_blank(),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.8,"cm"))



tag_facet2 <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.3, vjust = 1.2, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = "bold", family = family, inherit.aes = FALSE)
}

tags = letters[1:14] #c(paste0("a",1:7),paste0("b",1:7))
cplot.partrela.t5.withtag = tag_facet2(cplot.partrela.t5, tag_pool = tags)



filename.pic.partial = paste0(outpath.rf,"partial_pic/figpartial_1_t5combined.jpg")


ggsave(filename.pic.partial,cplot.partrela.t5.withtag,
       width = 22,height =9,units = "cm",dpi = 600)




################ 
#sixth loop for plotting partial relation seperate

colpalettes<-unique(c(pal_simpsons("springfield")(16),pal_gsea("default")(12),
                      pal_startrek("uniform")(7),pal_lancet("lanonc")(9)))




set.seed(222)
color.xvars = sample(colpalettes,length(alldrivernames))

color.table = data.frame(xvar = alldrivernames,
                         plotname = renamedr(alldrivernames),
                         color = color.xvars,
                         od = 1:length(alldrivernames))

color.table$color[5] = "#ABA80EFF"
color.table$color[9] = "#5F08A1FF"


cplotdata.t5 = c()


for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
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
    
    
    ndviloc = wherestrs(partialrelation$xname,c("ndvi"))
    
    partialrelation$xvar[ndviloc] = partialrelation$xvar[ndviloc]/1000
    
    plotdata.partrela = partialrelation
    
    
    
    
    xvar.temp = unique(plotdata.partrela$xname)
    color.table.loc = wherestrs(color.table$xvar,xvar.temp)
    color.table.temp = color.table[color.table.loc,]
    
    
    plotdata.partrela$xname = renamedr(plotdata.partrela$xname)
    plotdata.partrela$xname = addunit(plotdata.partrela$xname)
    
    plotdata.partrela$xname = factor(plotdata.partrela$xname,
                                     levels = addunit(color.table.temp$plotname))
    
    
    #####rescale PPT
    
    whereppt = wherestr(plotdata.partrela$xname,"PPT")
    plotdata.partrela$xvar[whereppt] = plotdata.partrela$xvar[whereppt]*12
    
    
    tagletter = letters[abs(mtype-2)*7+bsnum]
    
    plot.partrela = ggplot(plotdata.partrela,
                           aes(x=xvar,y=yvar,
                               color = xname, fill = xname)) + 
      facet_wrap(vars(xname),nrow = 3,scales = "free_x")+
      geom_line(linewidth=0.8)+
      xlab("")+
      ylab(paste0(yyname,"\n",fullname))+
      labs(tag = paste0(" ",tagletter))+
      scale_y_continuous(n.breaks = 4,
                         labels = scales::number_format(accuracy = 0.001,
                                                        style_positive = "plus",
                                                        style_negative = "minus"))+
      scale_x_continuous(n.breaks = 4)+
      scale_color_manual(values = color.table.temp$color)+
      theme(panel.background = element_rect(color = "black",fill = "white"),
            panel.spacing.x = unit(4, "mm"),
            panel.grid = element_blank(),
            axis.text.y = element_text(color = "black",
                                       margin = ggplot2::margin(t = 0, r = 0.1,
                                                                b = 0, l = 0.1,
                                                                unit = "cm")),
            axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.2,
                                                                 b = 0, l = 0.2,
                                                                 unit = "cm")),
            plot.title = element_text(size = 11),
            plot.tag = element_text(face = "bold",
                                    margin = ggplot2::margin(t = 0.1, r = -0.5,
                                                             b = 0, l = 0.1,
                                                             unit = "cm")),
            plot.margin = ggplot2::margin(r = 0.2,
                                          unit = "cm"),
            legend.position="none",
            legend.title = element_blank())
    
    
    # plot.partrela.fix <- set_panel_size(plot.partrela,
    #                           width  = unit(2.7, "cm"),
    #                           height = unit(1.5, "cm"))
    # 
    
    
    
    
    filename.pic.partial = paste0(outpath.rf,"partial_pic/figpartial_",
                                  tagletter,"_",
                                  shortname,"_",yyname,".jpg")
    
    
    ggsave(filename.pic.partial,plot.partrela,
           width = 25,height = 11,units = "cm",dpi = 600)
    
    # jpeg(filename.pic.partial,width = 18,
    #      height = 12,units = "cm",
    #      res = 600)
    # grid.newpage()
    # grid.draw(plot.partrela.fix)
    # pushViewport(viewport(just = "top"))
    # dev.off()
    
  }
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(resultframe1)
  rm(resultframe2)
  rm(rf1)
  rm(rf2)
  gc()
  
}




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
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
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
    
    
    
    
    data.pred.test = predict(rfm, newdata = data.test)
    data.pred.train = predict(rfm, newdata = data.train)
    data.crossverifiy = data.frame(actual = c(data.test[,mtype],data.train[,mtype]),
                                   predict = c(data.pred.test,data.pred.train),
                                   subset = c(rep("test",length(data.pred.test)),
                                              rep("train",length(data.pred.train))))
    value.max = max(c(data.crossverifiy$actual,data.crossverifiy$predict))
    value.min = min(c(data.crossverifiy$actual,data.crossverifiy$predict))
    
    plot.range = c(value.min + (value.max-value.min)*0,
                   value.max - (value.max-value.min)*0)
    
    rsq.on.test = rsq.ap(data.test[,mtype],data.pred.test)
    rsq.on.train = rsq.ap(data.train[,mtype],data.pred.train)
    
    text.x = plot.range[2]-(plot.range[2]-plot.range[1])*0.21
    text.y1 = plot.range[1]+(plot.range[2]-plot.range[1])*0.09  
    text.y2 = plot.range[1]+(plot.range[2]-plot.range[1])*0.29 
    
    
    data.crossverifiy$subset = factor(data.crossverifiy$subset,
                                      levels = c("test","train"))
    
    
    plot.cv = ggplot(data.crossverifiy, aes(x = predict, y = actual, group = subset)) + 
      stat_density_2d(geom = "polygon",aes(alpha = ..level.., fill = subset))+
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
  rm(resultframe1)
  rm(resultframe2)
  rm(rf1)
  rm(rf2)
  gc()
  
}



library(patchwork)



patch.combine.cv = list_mean[[1]]+list_mean[[2]]+list_mean[[3]]+list_mean[[4]]+
  list_mean[[5]]+list_mean[[6]]+list_mean[[7]]+plot_spacer()+
  list_trend[[1]]+list_trend[[2]]+list_trend[[3]]+list_trend[[4]]+
  list_trend[[5]]+list_trend[[6]]+list_trend[[7]]+plot_spacer()+
  plot_layout(ncol = 4)



outfile_patch = paste0(outpath.rf,"cv_pic/figcv_1combined_new.jpg")
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
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  
  
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



###################
#Additional loop for plotting importance in horizontal distribution of factors



color.relaimpo = c("#ab3838","#2f3f8a","#c4a117")


list_trend = list()
list_mean = list()


for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  
  
  #plotdrvnames = renamedr(originaldrvnames)
  
  
  load(file = paste0(outpath.rf,"datatrain_",shortname,".rda"))
  load(file = paste0(outpath.rf,"datatest_",shortname,".rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"RFmodel_",shortname,"_","AR1_mean.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_trend.rda"))
  load(file = paste0(outpath.rf,"Partialrela_",shortname,"_","AR1_mean.rda"))
  
  
  
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
    
    
    data.pred.test = predict(rfm, newdata = data.test)
    data.rsq = data.frame(actual = data.test[,mtype],predict = data.pred.test)
    rsq.test = rsq.ap(data.rsq$actual, data.rsq$predict)
    
    
    
    statform.partrela = c()
    
    for (drv in alldrivernames) {
      condition0 = (partialrelation$xname == drv)&
        (partialrelation$yname == yyname)&
        (partialrelation$bsname == shortname)
      data.trend = partialrelation[condition0,]
      mktest.result = cor.test(x=data.trend$xvar,y=data.trend$yvar, 
                               method ="pearson")
      p.val = mktest.result$p.value
      tau = mktest.result$estimate
      
      if(p.val<0.05&tau>0){
        rela = "positive"
      }else if(p.val<0.05&tau<0){
        rela = "negative"
      }else{
        rela = "non-monotone"
      }
      
      statform.partrela = rbind(statform.partrela,
                                data.frame(relation = rela, xvar = drv,
                                           yvar = yyname, biome = shortname))
      
    }
    
    
    
    plotdata.relaimpo0 = data.frame(xvar = rownames(rfm$importance),importance = rfm$importance[,1])
    rownames(plotdata.relaimpo0) = NULL
    
    plotdata.relaimpo = merge(plotdata.relaimpo0,statform.partrela,by = "xvar")
    plotdata.relaimpo = cbind(plotdata.relaimpo,
                              data.frame(rsq.on.test = rep(rsq.test,nrow(plotdata.relaimpo))))
    
    plotdata.relaimpo$relation = factor(plotdata.relaimpo$relation,
                                        levels = c("positive",
                                                   "negative",
                                                   "non-monotone"))
    
    
    plotdata.relaimpo$xvar = renamedr(plotdata.relaimpo$xvar)
    
    yrng <- range(plotdata.relaimpo$importance)
    textloc.x = length(alldrivernames) - 2
    textloc.y = yrng[2]*0.6
    
    
    tick.loc = signif(yrng[2]*0.6,1)
    
    text.title = paste0( fullname,"  ",yyname)
    
    
    plot.relaimpo = ggplot(plotdata.relaimpo,
                           aes(x=reorder(xvar, importance,decreasing =T),y=importance,
                               color = relation, fill = relation)) + 
      geom_bar(stat = 'identity', width = 0.3) +
      scale_color_manual(values =  color.relaimpo)+
      scale_fill_manual(values =  color.relaimpo)+
      scale_y_continuous(breaks = c(0,tick.loc),
                         limits = c(0,yrng[2]*1.2),
                         expand = c(0,NA))+
      xlab(NULL)+
      labs(tag = paste0("   ",letters[abs(mtype-2)*7+bsnum]))+
      annotate(geom = "text", x = 9, y = yrng[2], 
               label = text.title,
               size = 2.6)+
      annotate(geom = "text", x = textloc.x, y = textloc.y, 
               label = sprintf("italic(R)^2 == '%0.2f'",rsq.test),
               size = 2.6,parse = T)+
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(),
            axis.text.y = element_text(color = "gray20",size = 6),
            axis.title.y = element_text(color = "black",size = 6),
            axis.text.x = element_text(color = "black",size = 7,
                                       angle = 60, vjust = 1,hjust = 1),
            plot.tag.position = c(0,1),
            plot.tag = element_text(face = "bold",
                                    hjust = 0,
                                    vjust = 0.2,
                                    size = 10),
            legend.position="none")
    
    
    filename.pic.relaimpo = paste0(outpath.rf,"relaimpo_pic/figrelaimpo_",
                                   shortname,"_",yyname,".jpg")
    
    
    ggsave(filename.pic.relaimpo,plot.relaimpo,
           width = 8,height = 3.5,units = "cm",dpi = 600)
    
    
    if(mtype ==1 & bsnum==1){
      plot.relaimpo = plot.relaimpo + theme(legend.position="bottom",
                                            legend.title = element_blank(),
                                            legend.text = element_text(color = "black",size = 10),
                                            legend.key.height = unit(0.3,"cm"),
                                            legend.key.width = unit(0.2,"cm"))
      
      legendplot = as_ggplot(get_legend(plot.relaimpo))
    }
    
    
    if(mtype == 1){
      list_trend = c(list_trend,list(plot.relaimpo))
    }else{
      list_mean = c(list_mean,list(plot.relaimpo))
    }
    
    
    
  }
  
  rm(anadata)
  rm(anadata1)
  rm(data.test)
  rm(data.train)
  rm(resultframe1)
  rm(resultframe2)
  rm(rf1)
  rm(rf2)
  gc()
  
}

library(patchwork)
library(ggpubr)

patch.combine1 = c()
for (xx in 1:7) {
  
  patch.combine1 <- patch.combine1/
    ((list_mean[[xx]]+ theme(legend.position = "none"))|
       (list_trend[[xx]]+ theme(legend.position = "none")))
}


patch.combine1 <- (patch.combine1/legendplot) 

outfile_patch = paste0(outpath.rf,"relaimpo_pic/figrelaimpo_1combined.jpg")


jpeg(outfile_patch, units="cm",
     width= 9.2*2,
     height= 3.1*7,
     res = 600)

plot.new()

plot(patch.combine1)

dev.off()



################ 
#Additional loop for AR1 mean vs trend

library(mgcv)
library(patchwork)

abscolors = c("#267300","#a9a800","#ff5500","#00a884","#ffaa01","#004da7","#ff01c4")

list_fig = list()

for(bsnum in c(1:7)){
  
  fullname = biomefull[bsnum]
  shortname = biomeshort[bsnum]
  
  print(paste0(fullname," is under analysis!"))
  print(Sys.time())
  
  load(paste0("E:/AfricaResearch/StabilityTrend/data/ANA_result/ANA3_biostate/anadata_",
              shortname,".rda"))
  
  anadata$tgdp = anadata$tgdp/anadata$mgdp
  anadata$tpop = anadata$tpop/anadata$mpop
  
  anadata1 = anadata
  
  drivers = anadata1[,3:ncol(anadata1)]
  alldrivernames = colnames(drivers)
  
  alldrivernames = delstr(alldrivernames,c("tclimaet","tclimai","tclimdef",
                                           "tclimpet","tclimsoil","mclimaet",
                                           "tclimtmin","tclimtmax",
                                           "tpop","tgdp",
                                           "mclimai","mclimdef",
                                           "mclimpet","mclimsoil",
                                           "mclimtmin","mclimtmax",
                                           "mpop","mgdp",
                                           "mndvi","tndvi"))
  
  f0 = ggplot(data = anadata1, mapping = aes(x = AR1_mean, y = AR1_trend))+
    geom_smooth(color = abscolors[bsnum],
                method = "gam", 
                formula = y ~ s(x, bs = "cs"),
                level = 0.95)+
    geom_hline(yintercept = 0,linetype = "dashed")+
    xlab("AR1 mean")+
    ylab("AR1 trend")+
    labs(tag = letters[bsnum],
         title = fullname)+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.position = "none",
          plot.tag = element_text(face = "bold"),
          plot.title = element_text(face = "plain",size=11))
  
  dataxx = anadata1[,1:2]
  
  mmm =  gam(AR1_trend~s(AR1_mean),data = dataxx)
  datanewx = data.frame(seq(range(dataxx$AR1_mean)[1],range(dataxx$AR1_mean)[2],
                        length.out = 100))
  colnames(datanewx) = "AR1_mean"

  datapi = predict.gam(mmm,newdata = datanewx,interval = "prediction",level = 0.95)
  
  
  
  
  list_fig = c(list_fig,list(f0))
  
  rm(anadata)
  rm(anadata1)
  
  
  gc()
  print(paste0(fullname," is finished!"))
  
}

pf = list_fig[[1]]+list_fig[[2]]+list_fig[[3]]+
  list_fig[[4]]+list_fig[[5]]+list_fig[[6]]+list_fig[[7]]+
  plot_layout(nrow = 2)

outfile_mt = paste0(outpath3,"mean_trend.jpg")
jpeg(outfile_mt, units="cm",
     width= 22,
     height= 11,
     res = 600)

plot.new()

plot(pf)

dev.off()

######################################
##test function for data y range by x intervals 


breakquantile <- function(data0, cx = 2, cy = 1, bs = 100,level = 0.95){
  
  rangex = range(data0[,cx])
  breakx = seq(rangex[1],rangex[2],length.out = bs+1)
  
  vups = c()
  vdowns = c()
  vmeans = c()
  for (i in 1:bs) {
    itvmin = breakx[i]
    itvmax = breakx[i+1]
    
    datai = data0[data0[,cx]>=itvmin&data0[,cx]<=itvmax,]
    
    itvy = datai[,cy]
    vup = quantile(itvy,level+(1-level)/2)
    vdown = quantile(itvy,(1-level)/2)
    
    vups = c(vups,vup)
    vdowns = c(vdowns,vdown)
    vmeans = c(vmeans,mean(c(itvmin,itvmax)))
  }
  df = data.frame(vmeans,vups,vdowns)
  
  return(df)
}




