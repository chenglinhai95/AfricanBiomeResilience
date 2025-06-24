
library("raster")
library("terra")
library("foreach")
library("doParallel")
library("doSNOW")
library("signal")
library("zoo") 
library("trend") 
library("Kendall")

###functions
##############################################

AR1mean <- function(vec0, min.dlr=0.2, min.sample = 10){
  vec1 = as.vector(vec0)
  
  x = vec1[1:(length(vec1)-1)]
  y = vec1[2:length(vec1)]
  
  df0 = data.frame(x,y)
  
  if (nrow(na.omit(df0))/nrow(df0)< min.dlr){
    return(NA)
  }
  
  if (nrow(na.omit(df0))< min.sample){
    return(NA)
  }

  lm.mod = lm(y~x, df0, na.action = na.omit)
  ar1 = lm.mod$coefficients[['x']]
  return(ar1)
}



compute_sen_slope <- function(time_series) {
  n <- length(time_series)
  slopes <- numeric(n * (n - 1) / 2)
  idx <- 1
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      slopes[idx] <- (time_series[j] - time_series[i]) / (j - i)
      idx <- idx + 1
    }
  }
  
  return(median(slopes))
}



sliding_window_AR1_complete <- function(ts, window_size = 23 * 6, step_size = 23) {
  n <- length(ts)
  ar1_values <- c() 
  for (start_idx in seq(1, n - window_size+1, by = step_size)) {
    locend = start_idx + window_size - 1
    if(locend>n){
      locend = n
    }
    window <- ts[start_idx:locend]
    ar1_values <- c(ar1_values, AR1mean(window))
  }
  return(ar1_values)
}



AR1trend <- function(vec0, winsize, step = 23){
  ar1.seq = sliding_window_AR1_complete(ts = vec0, window_size = winsize, step_size = step)
  original_slope <- compute_sen_slope(ar1.seq)
  mk_result <- MannKendall(ar1.seq)
  mkpvalue <- mk_result$sl[1]
  df.return = data.frame(senslope = original_slope,
                         mkpvalue = mkpvalue)
  return(df.return)
}



AR1meanNtrend <- function(vec0,trend.window = c(23*4,23*6,23*8), freq.ts = 23){
  ts.length = length(vec0)
  vec1 = vec0/10000
  vec1 = tanh(vec1^2) ##transform to kNDVI
  
  xtt = ts(vec1,start = 1,frequency = 23)
  stt = stl(xtt,s.window = "periodic",s.degree = 0,t.window = 23*10,
            t.degree = 1,robust = T, na.action = na.omit)
  
  res.seq = stt$time.series[,3]

  
  ar1m = AR1mean(res.seq)
  ar1t4yr = AR1trend(res.seq,winsize = trend.window[1])
  ar1t6yr = AR1trend(res.seq,winsize = trend.window[2])
  ar1t8yr = AR1trend(res.seq,winsize = trend.window[3])
  
  cnames = c("AR1m",
  "AR1t_senslope_4yr","AR1t_mkpvalue_4yr",
  "AR1t_senslope_6yr","AR1t_mkpvalue_6yr",
  "AR1t_senslope_8yr","AR1t_mkpvalue_8yr")
  result.return = cbind(ar1m,ar1t4yr,ar1t6yr,ar1t8yr)
  colnames(result.return) = cnames
  return(result.return)
}


AR1meanNtrend_batch <- function(mdata){
  if(is.vector(mdata)){
    mdata = matrix(mdata,nrow = 1)
  }
  
  irs = nrow(mdata)
  result.matrix = c()
  for (i in 1:irs){
    result.matrix = rbind(result.matrix, AR1meanNtrend(mdata[i,]))
  }
  return(result.matrix)
}

setbatch <- function(ntotal,boxsize = 20){
  result.m = c()
  for (i in seq(1,ntotal,boxsize)) {
    s1= i
    s2 = i+boxsize-1
    if(s2>ntotal){
      s2 = ntotal
    }
    result.m = rbind(result.m,cbind(s1,s2))
  }
  return(result.m)
}

################################

setwd("D:/Earthfuture_reanalysis/istsg_AF")
outpath <- "D:/Earthfuture_reanalysis/istsg_result/"

refrast = rast("D:/Earthfuture_reanalysis/AFABSs.tif")

#####align raster data
# for (yr in 2001:2020) {
#   rr0 = rast(paste0("Istsg_real",yr,"_AF.tif"))
#   rr1 = project(rr0,refrast,method = "bilinear")
#   rr2 = crop(rr1,refrast)
#   writeRaster(rr2,paste0("AF_Istsg_Aligned_",yr,"_AF.tif"))
#   print(paste0("Year ",yr, " contains " ,nlyr(rr2)," layers!")) 
# }



###################

ABSmap <- raster("D:/Earthfuture_reanalysis/AFABSs.tif")
ext_ref = extent(ABSmap)

m_ABSmap = as.matrix(ABSmap)
matrix_loc = which(m_ABSmap>=1&m_ABSmap<=7)

##get file list in order
rnames <- paste0("AF_Istsg_Aligned_",2001:2020,"_AF.tif")


#example for clip
# r.stack.test = stack(rnames[20])
# r.test = r.stack.test[[23]]
# nrow(r.test)
# ncol(r.test)
# nrow(ABSmap)
# ncol(ABSmap)

# ##save data table
# T1 = Sys.time()
# dataforana = matrix(NA, nrow = length(matrix_loc), ncol = length(rnames)*23)
# colnumber = 0
# for (yr in 2001:2020) {
#   r.stack0 = stack(paste0("AF_Istsg_Aligned_",yr,"_AF.tif"))
#   for (lyr in 1:23) {
#     r0 = r.stack0[[lyr]]
#     m0 = as.matrix(r0)
#     vec0 = m0[matrix_loc]
#     colnumber = colnumber+1
#     dataforana[,colnumber] = vec0
#     print(paste0("AF_Istsg_Aligned_",yr,"_AF.tif",' at layer ', lyr))
#   }
# }
# T2 = Sys.time()
# print(T2-T1)
# saveRDS(dataforana,
#         file = paste0(outpath,'dataforana.rds'),
#         compress = F)


# testanadata = dataforana[420:500,]
# testresult = AR1meanNtrend_batch(testanadata)



##extract subset of data

total.rows = length(matrix_loc)
breaknum = 30
breakrows = round(seq(from = 0, to = total.rows,length.out = breaknum+1))


# print(Sys.time())
# for(datapart in 1:breaknum){
#   loc.start = breakrows[datapart]+1
#   loc.end = breakrows[datapart+1]
#   dataforana.subset = dataforana[loc.start:loc.end,]
#   saveRDS(dataforana.subset,
#           file = paste0(outpath,'dataforana_sub',datapart,'.rds'),
#           compress = F)
#   print(paste0('Data subset ',datapart,' is saved!'))
#   print(Sys.time())
# }
# rm(dataforana)
# gc()


##AR1 calculation

result.all = c()

for(datapart in 1:breaknum){
  
  print(paste0('Data subset ',datapart,' start at ',Sys.time()))
  dataforana.subset = readRDS(file = paste0(outpath,'dataforana_sub',datapart,'.rds'))

  batchtable = setbatch(nrow(dataforana.subset),boxsize = 200)

  nseries = nrow(batchtable)
  
  nc = 16
  cl <- makeCluster(nc,type = "SOCK")
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = nseries,style = 3)
  progress = function(n) setTxtProgressBar(pb,n)
  opts = list(progress = progress)
  
  
  result.subset <- foreach(x = 1:nseries,
                           .combine='rbind',
                           .inorder = TRUE,
                           .packages=c('signal','zoo',
                                       'Kendall','trend'),
                           .options.snow = opts,
                           .errorhandling = "pass") %dopar% 
    {
      s1 = batchtable[x,1]
      s2 = batchtable[x,2]
      m0 = dataforana.subset[s1:s2,]
      rr = AR1meanNtrend_batch(m0)
      return(rr)
    }
  
  stopCluster(cl)
  saveRDS(result.subset, file = paste0(outpath,'result_sub',datapart,'.rds'))
  result.all = rbind(result.all,result.subset)
  print(paste0('Data subset ',datapart,' finished at ',Sys.time()))
}

saveRDS(result.all, file = paste0(outpath,'result_total.rds'))



blankmatrix = matrix(NA, nrow = nrow(m_ABSmap), ncol = ncol(m_ABSmap))

for (bandname in colnames(result.all)) {

  mm00 = blankmatrix
  mm00[matrix_loc] = as.vector(result.all[[bandname]])
  rr00 = raster(mm00)
  crs(rr00) = crs(ABSmap)
  extent(rr00) = extent(ABSmap)  
  if(bandname == "senslope"){
    bandname = "ar1t_senslope"
  }
  if(bandname == "pvalue"){
    bandname = "FFTpvalue"
  }
  writeRaster(rr00, file = paste0(outpath,'AF_istsg_',bandname,'.tif'),
              overwrite=TRUE)
}
