

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
    df0[,c] = 1000*(df0[,c]-mins[c])/(maxs[c]-mins[c])
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


####################################################################

setwd(paste0("F:/AfricaResearch/StabilityTrend/data/AlignProcess/Pool_Aligned_data"))

outpath = "F:/AfricaResearch/StabilityTrend/paper_fig1_EF/"


map_afbs = raster("AFABSs.tif")
map_afbs_ar1mean = raster("AF_istsg_ar1m.tif")
map_afbs_ar1trend = raster("AF_istsg_AR1t_senslope_6yr.tif")
mapar1trendpvalue = raster("AF_istsg_AR1t_mkpvalue_6yr.tif")



map_afbs_ar1trendsig = sign(map_afbs_ar1trend)
map_afbs_ar1trendsig[mapar1trendpvalue>=0.05] = 0
map_afbs_ar1trendsig[(map_afbs<1)] = NA
writeRaster(map_afbs_ar1trendsig,filename = "AF_istsg_ar1trend_sig.tif")

################################################# AR1 trend

sta_data =  data.frame(bstype = as.vector(map_afbs),
                       ar1mean = as.vector(map_afbs_ar1mean),
                       ar1trend = as.vector(map_afbs_ar1trend),
                       ar1trendsig = as.vector(map_afbs_ar1trendsig))

sta_data = sta_data[(sta_data$bstype>=1),]

sta_data = na.omit(sta_data)

sta_data0 = sta_data

#sta_data = cbind.data.frame(sta_data,note = rep("",nrow(sta_data)))


pic_temp = ggplot()+
  geom_density(data = sta_data,mapping = aes(x = ar1trend,after_stat(count),
                                             group = ar1trendsig,
                                             fill = ar1trendsig,
                                             color = ar1trendsig),position = "stack")+
  facet_grid(bstype ~ .,switch="both",scales="free_y")

yranges = ggplot_build(pic_temp)$layout$panel_scales_y



staform = c()
incnotes = c()
decnotes = c()

for (bsnum in c(1:7)) {
  temdata = sta_data[sta_data$bstype==bsnum,]
  total = nrow(temdata)
  siginc = sum(temdata$ar1trendsig==1)
  sigdec = sum(temdata$ar1trendsig==-1)
  noteinc = paste0(sprintf("%.1f",siginc*100/total),"%")
  notedec = paste0(sprintf("%.1f",sigdec*100/total),"%")
  
  ym = mean(yranges[[bsnum]]$range$range)
  
  df0 = data.frame(x = c(-0.04,0.04),y = c(ym,ym),
                   note = c(notedec,noteinc),bstype = c(bsnum,bsnum))
  
  staform = rbind(staform,df0)
  incnotes = c(incnotes,noteinc)
  decnotes = c(decnotes,notedec)
}





biomelabs = c('Closed\nForest','Woody\nSavanna','Savanna',
              'Dense\nGrassland','Sparse\nGrassland',
              'Closed\nShrubland','Open\nShrubland')

sta_data$bstype = factor(sta_data$bstype,levels = c(1:7),labels = biomelabs)
sta_data$ar1trendsig = factor(sta_data$ar1trendsig,levels = c(1,-1,0),
                              labels = c("Sig. Increase",
                                         "Sig. Decrease",
                                         "Non-sig"))


staform$bstype = factor(staform$bstype,levels = c(1:7),labels = biomelabs)


sigcolors = c("#ca5862","#11217e","#ccbf87")

pic_trend = ggplot()+
  geom_density(data = sta_data,mapping = aes(x = ar1trend,after_stat(count),
                                group = ar1trendsig,
                                fill = ar1trendsig,
                                color = ar1trendsig),position = "stack")+
  geom_text(data = staform, mapping = aes(x = x,y = y, label = note), size = 4)+
  facet_grid(bstype ~ .,switch="both",scales="free_y")+
  geom_hline(yintercept = 0)+
  scale_x_continuous(breaks = c(-0.05,0,0.05))+
  scale_fill_manual(values = sigcolors)+
  scale_color_manual(values = sigcolors)+
  xlab("AR1 trend")+
  theme(strip.placement = "inside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(color="black",angle = 0,size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="gray20",size = 10),
        axis.line.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

ggsave(filename = paste0(outpath,"AR1trend_statistic.jpg"),pic_trend,
       width = 3.8,height = 3.8,units = "in",dpi = 600)


##################################################### AR1 mean


abscolors = c("#267300","#a9a800","#ff5500","#00a884","#ffaa01","#004da7","#ff01c4")


pic_mean = ggplot(data = sta_data,mapping = aes(x = ar1mean,
                                                fill = bstype,
                                                color = bstype))+
  geom_density()+
  facet_grid(bstype ~ .,switch="both",scales="free_y")+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = abscolors)+
  scale_color_manual(values = abscolors)+
  xlab("AR1 mean")+
  theme(strip.placement = "inside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(color="black",angle = 0,size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="gray20",size = 10),
        axis.line.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")


ggsave(filename = paste0(outpath,"AR1mean_statistic.jpg"),pic_mean,
       width = 3.8,height = 3.8,units = "in",dpi = 600)

#########################################################  AR1 sig


sigcolors1 = c("#11217e","#ccbf87","#ca5862")
abscolors = c("#267300","#a9a800","#ff5500","#00a884","#ffaa01","#004da7","#ff01c4")

sta_data1 =sta_data0

sta_data1$ar1trendsig = factor(sta_data1$ar1trendsig,levels = c(-1,0,1),
                              labels = c("Sig. Decrease",
                                         "Non-sig",
                                         "Sig. Increase"))
biomelabs = c('Closed\nForest','Woody\nSavanna','Savanna',
              'Dense\nGrassland','Sparse\nGrassland',
              'Closed\nShrubland','Open\nShrubland')

sta_data1$bstype = factor(sta_data1$bstype,levels = c(1:7),labels = biomelabs)


sta_data1inc = sta_data1[sta_data1$ar1trendsig =="Sig. Increase", ]
sta_data1dec = sta_data1[sta_data1$ar1trendsig =="Sig. Decrease", ]





pic_pie = ggplot(sta_data1, aes(x = "", fill = ar1trendsig)) + 
  geom_bar(position = "fill") +
  # geom_text(aes(x = 1.7,label = scales::percent(after_stat(count)/sum(after_stat(count)))),
  #           stat = 'count',
  #           position = position_fill(vjust = 0.5),
  #           size = 8, color = "black") +
  scale_fill_manual(values = sigcolors1)+
  scale_color_manual(values = sigcolors1)+
  coord_polar(theta = 'y')+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())



percent_data = c()

for (bs0 in biomelabs) {
  bstype = bs0
  yinc = sum(sta_data1inc$bstype == bs0)/nrow(sta_data1inc)
  ydec = -sum(sta_data1dec$bstype == bs0)/nrow(sta_data1dec)
  yinclab = paste0(sprintf("%.1f",abs(yinc)*100),"%")
  ydeclab = paste0(sprintf("%.1f",abs(ydec)*100),"%")
  df0 = data.frame(bstype,yinc,ydec,yinclab,ydeclab)
  
  percent_data = rbind(percent_data,df0)
}

percent_data$bstype = factor(percent_data$bstype,levels = rev(biomelabs))
abscolors1 = rev(abscolors) 


pic_bar = ggplot(data = percent_data, mapping = aes(fill = bstype))+
  geom_bar(mapping = aes(x = bstype, y = yinc),stat = 'identity')+
  geom_bar(mapping = aes(x = bstype, y = ydec),stat = 'identity')+
  geom_text(mapping = aes(x = bstype, y = yinc+0.2, label = yinclab),
            size = 8, color = "black")+
  geom_text(mapping = aes(x = bstype, y = ydec-0.2, label = ydeclab),
            size = 8, color = "black")+
  geom_hline(yintercept = 0,linewidth = 1.5)+
  scale_fill_manual(values = abscolors1)+
  scale_y_continuous(limits = c(-0.6,0.6))+
  coord_flip()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())


ggsave(filename = paste0(outpath,"AR1sig_contribution.jpg"),pic_bar,
       width = 6,height = 3.5,units = "in",dpi = 600)


ggsave(filename = paste0(outpath,"AR1sig_percent.jpg"),pic_pie,
       width = 3.5,height = 3.5,units = "in",dpi = 600)




