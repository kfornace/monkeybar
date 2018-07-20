library(maptools)
library(raster)
require(fields)
require(parallel)
require(dismo)
require(snowfall)

newproj<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r.km<-c(100,200,500,1000,2000,3000,4000,5000,7500,10000,20000)

# North & South
k1<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[1],".csv"),header=T)
k2<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[2],".csv"),header=T)
k3<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[3],".csv"),header=T)
k4<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[4],".csv"),header=T)
k5<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[5],".csv"),header=T)
k6<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[6],".csv"),header=T)
k7<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[7],".csv"),header=T)
k8<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[8],".csv"),header=T)
k9<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[9],".csv"),header=T)
k10<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept",r.km[10],".csv"),header=T)
k11<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_sept16",r.km[11],".csv"),header=T)

# Banggi
bk1<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[1],".csv"),header=T)
bk2<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[2],".csv"),header=T)
bk3<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[3],".csv"),header=T)
bk4<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[4],".csv"),header=T)
bk5<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[5],".csv"),header=T)
bk6<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[6],".csv"),header=T)
bk7<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[7],".csv"),header=T)
bk8<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[8],".csv"),header=T)
bk9<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[9],".csv"),header=T)
bk10<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[10],".csv"),header=T)
bk11<-read.csv(paste("kds.resol.BANGGI_BIOLdist_sept",r.km[11],".csv"),header=T)

# To match all columns
bk1$in_city<-NA
bk2$in_city<-NA
bk3$in_city<-NA
bk4$in_city<-NA
bk5$in_city<-NA
bk6$in_city<-NA
bk7$in_city<-NA
bk8$in_city<-NA
bk9$in_city<-NA
bk10$in_city<-NA
bk11$in_city<-NA

# rearrange
k1.<-k1[,c(1:13,15:65,14)]
k2.<-k2[,c(1:13,15:65,14)]
k3.<-k3[,c(1:13,15:65,14)]
k4.<-k4[,c(1:13,15:65,14)]
k5.<-k5[,c(1:13,15:65,14)]
k6.<-k6[,c(1:13,15:65,14)]
k7.<-k7[,c(1:13,15:65,14)]
k8.<-k8[,c(1:13,15:65,14)]
k9.<-k9[,c(1:13,15:65,14)]
k10.<-k10[,c(1:13,15:65,14)]
k11.<-k11[,c(1:13,15:65,14)]

# add pop.den to k1. and bk1
popk1<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[1],".csv"),header=T)
popk2<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[2],".csv"),header=T)
popk3<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[3],".csv"),header=T)
popk4<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[4],".csv"),header=T)
popk5<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[5],".csv"),header=T)
popk6<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[6],".csv"),header=T)
popk7<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[7],".csv"),header=T)
popk8<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[8],".csv"),header=T)
popk9<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[9],".csv"),header=T)
popk10<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[10],".csv"),header=T)
popk11<-read.csv(paste("kds.resol.NORTHSOUTH_BIOLdist_popDENnew",r.km[11],".csv"),header=T)

# Banggi
popb1<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[1],".csv"),header=T)
popb2<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[2],".csv"),header=T)
popb3<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[3],".csv"),header=T)
popb4<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[4],".csv"),header=T)
popb5<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[5],".csv"),header=T)
popb6<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[6],".csv"),header=T)
popb7<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[7],".csv"),header=T)
popb8<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[8],".csv"),header=T)
popb9<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[9],".csv"),header=T)
popb10<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[10],".csv"),header=T)
popb11<-read.csv(paste("kds.resol.BANGGI_BIOLdist_popDENnew",r.km[11],".csv"),header=T)

# Include pop den
k1.$pop.den <- popk1$pop.den.NEW
bk1$pop.den <- popb1$pop.den.NEW
k2.$pop.den <- popk2$pop.den.NEW
bk2$pop.den <- popb2$pop.den.NEW
k3.$pop.den <- popk3$pop.den.NEW
bk3$pop.den <- popb3$pop.den.NEW
k4.$pop.den <- popk4$pop.den.NEW
bk4$pop.den <- popb4$pop.den.NEW
k5.$pop.den <- popk5$pop.den.NEW
bk5$pop.den <- popb5$pop.den.NEW
k6.$pop.den <- popk6$pop.den.NEW
bk6$pop.den <- popb6$pop.den.NEW
k7.$pop.den <- popk7$pop.den.NEW
bk7$pop.den <- popb7$pop.den.NEW
k8.$pop.den <- popk8$pop.den.NEW
bk8$pop.den <- popb8$pop.den.NEW
k9.$pop.den <- popk9$pop.den.NEW
bk9$pop.den <- popb9$pop.den.NEW
k10.$pop.den <- popk10$pop.den.NEW
bk10$pop.den <- popb10$pop.den.NEW
k11.$pop.den <- popk11$pop.den.NEW
bk11$pop.den <- popb11$pop.den.NEW

# Combine
a1<-rbind(k1.,bk1)
a2<-rbind(k2.,bk2)
a3<-rbind(k3.,bk3)
a4<-rbind(k4.,bk4)
a5<-rbind(k5.,bk5)
a6<-rbind(k6.,bk6)
a7<-rbind(k7.,bk7)
a8<-rbind(k8.,bk8)
a9<-rbind(k9.,bk9)
a10<-rbind(k10.,bk10)
a11<-rbind(k11.,bk11)

a1t.<-a1[-which(is.na(a1$species)==F & a1$species!="Pk"),]
a2t.<-a2[-which(is.na(a2$species)==F & a2$species!="Pk"),]
a3t.<-a3[-which(is.na(a3$species)==F & a3$species!="Pk"),]
a4t.<-a4[-which(is.na(a4$species)==F & a4$species!="Pk"),]
a5t.<-a5[-which(is.na(a5$species)==F & a5$species!="Pk"),]
a6t.<-a6[-which(is.na(a6$species)==F & a6$species!="Pk"),]
a7t.<-a7[-which(is.na(a7$species)==F & a7$species!="Pk"),]
a8t.<-a8[-which(is.na(a8$species)==F & a8$species!="Pk"),]
a9t.<-a9[-which(is.na(a9$species)==F & a9$species!="Pk"),]
a10t.<-a10[-which(is.na(a10$species)==F & a10$species!="Pk"),]
a11t.<-a11[-which(is.na(a11$species)==F & a11$species!="Pk"),]

a1t..<-a1t.[-which(a1t.$no_mem==0),]
a2t..<-a2t.[-which(a2t.$no_mem==0),]
a3t..<-a3t.[-which(a3t.$no_mem==0),]
a4t..<-a4t.[-which(a4t.$no_mem==0),]
a5t..<-a5t.[-which(a5t.$no_mem==0),]
a6t..<-a6t.[-which(a6t.$no_mem==0),]
a7t..<-a7t.[-which(a7t.$no_mem==0),]
a8t..<-a8t.[-which(a8t.$no_mem==0),]
a9t..<-a9t.[-which(a9t.$no_mem==0),]
a10t..<-a10t.[-which(a10t.$no_mem==0),]
a11t..<-a11t.[-which(a11t.$no_mem==0),]

# all NS
kt<-list(a1t..,a2t..,a3t..,a4t..,a5t..,a6t..,a7t..,a8t..,a9t..,a10t..,a11t..)
ncs<-sum(kt[[1]]$loc_typ=="case")
enum <-which(kt[[1]]$loc_typ=="enum")
kmp <- which(kt[[1]]$loc_typ=="kmp")
pres <-which(kt[[1]]$loc_typ=="case")

#### Sampling
pab.reps<-100

# create data.frames
for(j in 1:pab.reps){
  abs.samp<-c(sample(kmp,(ncs/2),replace=F),sample(enum,(ncs/2),replace=F))

  sam.ds1<-rbind(kt[[1]][abs.samp,],kt[[1]][pres,])
  sam.ds2<-rbind(kt[[2]][abs.samp,],kt[[2]][pres,])
  sam.ds3<-rbind(kt[[3]][abs.samp,],kt[[3]][pres,])
  sam.ds4<-rbind(kt[[4]][abs.samp,],kt[[4]][pres,])
  sam.ds5<-rbind(kt[[5]][abs.samp,],kt[[5]][pres,])
  sam.ds6<-rbind(kt[[6]][abs.samp,],kt[[6]][pres,])
  sam.ds7<-rbind(kt[[7]][abs.samp,],kt[[7]][pres,])
  sam.ds8<-rbind(kt[[8]][abs.samp,],kt[[8]][pres,])
  sam.ds9<-rbind(kt[[9]][abs.samp,],kt[[9]][pres,])
  sam.ds10<-rbind(kt[[10]][abs.samp,],kt[[10]][pres,])
  sam.ds11<-rbind(kt[[11]][abs.samp,],kt[[11]][pres,])

  nl<-list(sam.ds1,sam.ds2,sam.ds3,sam.ds4,sam.ds5,sam.ds6,sam.ds7,sam.ds8,sam.ds9,sam.ds10,sam.ds11)
  
  # assign each point into one of these 10 zones - do this for each NS, as might want to use separately
  i<-1
    nl[[i]]$yfoldzone<-NA
    nl[[i]]$foldzone<-NA
    nl[[i]]$yfoldzone<-ifelse(nl[[i]]$coords.x2>7.05,1,ifelse(nl[[i]]$coords.x2>6.825,2,ifelse(nl[[i]]$coords.x2>6.64,3,ifelse(nl[[i]]$coords.x2>6.46,4,5))))
    nl[[i]][nl[[i]]$yfoldzone==1,]$foldzone<-ifelse(nl[[i]][nl[[i]]$yfoldzone==1,]$coords.x1<117.14,"a","b")
    nl[[i]][nl[[i]]$yfoldzone==2,]$foldzone<-ifelse(nl[[i]][nl[[i]]$yfoldzone==2,]$coords.x1<116.76,"c","d")
    nl[[i]][nl[[i]]$yfoldzone==3,]$foldzone<-ifelse(nl[[i]][nl[[i]]$yfoldzone==3,]$coords.x1<116.76,"e","f")
    nl[[i]][nl[[i]]$yfoldzone==4,]$foldzone<-ifelse(nl[[i]][nl[[i]]$yfoldzone==4,]$coords.x1<116.84,"g","h")
    nl[[i]][nl[[i]]$yfoldzone==5,]$foldzone<-ifelse(nl[[i]][nl[[i]]$yfoldzone==5,]$coords.x1<116.84,"i","j")
    nl[[i]]$foldzoneN<-nl[[i]]$foldzone
    nl[[i]]$foldzoneN<-ifelse(nl[[i]]$foldzone=="a" | nl[[i]]$foldzone=="b" , "AB",
                              ifelse(nl[[i]]$foldzone=="c" | nl[[i]]$foldzone=="d" , "CD",
                                     ifelse(nl[[i]]$foldzone=="h" | nl[[i]]$foldzone=="i" | nl[[i]]$foldzone=="j" , "HIJ",
                                            nl[[i]]$foldzone )))
    nl[[i]]$occ<-ifelse(nl[[i]]$loc_typ=="case",1,0)
  
  # variables that apply across spatial scales
  vars.mult<-c(
    "prop.fcCYanst",      
    "fcCYanst.F.par",     
    "prop.clCYanst",      
    "clCYanst.F.par",     
    "prop.fl.L5Y",        
    "fl.L5Y.F.par",      
    "prop.flCYanst",      
    "flCYanst.F.par",     
    "prop.fganst",       
    "fganst.F.par",       
    "mean.elev",          
    "mean.asp",          
    "mean.slp",           
    "mean.ndvi",         
    "sd.ndvi",
    "pop.den")
  
  mk<-match(vars.mult,names(nl[[1]]))
  
  nln1<-nl[[1]][,mk]
  names(nln1)<-paste(names(nln1),"_ns1",sep="")
  nln2<-nl[[2]][,mk]
  names(nln2)<-paste(names(nln2),"_ns2",sep="")
  nln3<-nl[[3]][,mk]
  names(nln3)<-paste(names(nln3),"_ns3",sep="")
  nln4<-nl[[4]][,mk]
  names(nln4)<-paste(names(nln4),"_ns4",sep="")
  nln5<-nl[[5]][,mk]
  names(nln5)<-paste(names(nln5),"_ns5",sep="")
  nln6<-nl[[6]][,mk]
  names(nln6)<-paste(names(nln6),"_ns6",sep="")
  nln7<-nl[[7]][,mk]
  names(nln7)<-paste(names(nln7),"_ns7",sep="")
  nln8<-nl[[8]][,mk]
  names(nln8)<-paste(names(nln8),"_ns8",sep="")
  nln9<-nl[[9]][,mk]
  names(nln9)<-paste(names(nln9),"_ns9",sep="")
  nln10<-nl[[10]][,mk]
  names(nln10)<-paste(names(nln10),"_ns10",sep="")
  nln11<-nl[[11]][,mk]
  names(nln11)<-paste(names(nln11),"_ns11",sep="")

  ds<-cbind(nln1,nln2,nln3,nln4,nln5,nln6,nln7,nln8,nln9,nln10,nln11)

  ds$coords.x1<-nl[[1]]$coords.x1
  ds$coords.x2<-nl[[1]]$coords.x2
  ds$occ<-nl[[1]]$occ
  ds$no_mem<-nl[[1]]$no_mem
  ds$age<-nl[[1]]$age
  ds$gend<-nl[[1]]$gend
  ds$name<-nl[[1]]$name
  ds$date<-nl[[1]]$date
  ds$mn_d_cs<-nl[[1]]$mn_d_cs
  ds$mn_d_km<-nl[[1]]$mn_d_km
  ds$job<-nl[[1]]$job
  ds$d.nearest.road<-nl[[1]]$d.nearest.road
  ds$foldzoneN<-nl[[1]]$foldzoneN
  ds$d.nearest.clinic<-nl[[1]]$d.nearest.clinic
  ds$elev<-nl[[1]]$elev
  ds$asp<-nl[[1]]$asp
  ds$slp<-nl[[1]]$slp
  
  dse <- ds
  dse$ofs<-1/dse$no_mem
  dse$fzn<-ifelse(dse$foldzoneN=="AB",1,ifelse(dse$foldzoneN=="CD",2,ifelse(dse$foldzoneN=="e",3,
              ifelse(dse$foldzoneN=="f",4,ifelse(dse$foldzoneN=="g",5,6)))))
  
  assign(paste("dse",j,sep=""),dse)
}

dsl<-list(dse1,dse2,dse3,dse4,dse5,dse6,dse7,dse8,dse9,dse10,
          dse11,dse12,dse13,dse14,dse15,dse16,dse17,dse18,dse19,dse20,
          dse21,dse22,dse23,dse24,dse25,dse26,dse27,dse28,dse29,dse30,
          dse31,dse32,dse33,dse34,dse35,dse36,dse37,dse38,dse39,dse40,
          dse41,dse42,dse43,dse44,dse45,dse46,dse47,dse48,dse49,dse50,
          dse51,dse52,dse53,dse54,dse55,dse56,dse57,dse58,dse59,dse60,
          dse61,dse62,dse63,dse64,dse65,dse66,dse67,dse68,dse69,dse70,
          dse71,dse72,dse73,dse74,dse75,dse76,dse77,dse78,dse79,dse80,
          dse81,dse82,dse83,dse84,dse85,dse86,dse87,dse88,dse89,dse90,
          dse91,dse92,dse93,dse94,dse95,dse96,dse97,dse98,dse99,dse100)

# BRT
sfInit(cpus = 4, parallel = TRUE)
sfLibrary(seegSDM)
model_list <- sfLapply(dsl,gbm.step,gbm.x = c(1:176), gbm.y = 179, family = "bernoulli", tree.complexity = 5, learning.rate = 0.002,bag.fraction = 0.5)
sfStop()

# stats
nvars<-176
nres<-11
list.stats<-matrix(nrow=pab.reps,ncol=8)
for(i in 1:pab.reps){
  list.stats[i,]<-as.numeric(as.character(model_list[[i]]$cv.statistics[c(1:6,9:10)]))
}
out.stats<-as.data.frame(list.stats)
names(out.stats)<-c("deviance.mean","deviance.se","correlation.mean","correlation.se",
                    "discrimination.mean","discrimination.se","cv.threshold","cv.threshold.se")
out.stats
hist(out.stats$discrimination.mean)

# Relative influences
pab.reps<-100
nvars<-176
nres<-11
top.vars<-16

relinf<-matrix(NA,nrow=nvars,ncol=1+pab.reps)
relinf.ds<-as.data.frame(relinf)
relinf.ds[,1]<-names(dse1[,c(1:nvars)])
names(relinf.ds)<-c("Var",1:pab.reps)
for(k in 1:pab.reps){
  tmp.sum<-summary(model_list[[k]],plotit=F)
  for(i in 1:length(relinf.ds[,1])){
    relinf.ds[i,k+1]<-tmp.sum$rel.inf[tmp.sum$var==relinf.ds$Var[i]]
  }
}

# plot
boxplot(t(relinf.ds[,2:101]),names=relinf.ds$Var,las=2,cex.axis=0.8)

# collate results and save to file
#write.table(out.stats,paste("NS_ALL_stats.out.mean.JSVwout50.txt"))
#write.table(relinf.ds,paste("NS_ALL_var.ri.out.JSVwout50.txt"))

outstats<-read.table("NS_ALL_stats.out.mean.JSVwout50.txt",header=T)
relinf.ds<-read.table("NS_ALL_var.ri.out.JSVwout50.txt",header=T)

nridf<-matrix(nrow=pab.reps,ncol=nvars)
for(i in 1:pab.reps){
  for(k in 1:nvars){
    nridf[i,k]<-relinf.ds[k,i+1]
  }
}

relinf.ds$Var
relinf.ds$varname<-relinf.ds$Var
relinf.ds$median<-NA
for(i in 1:length(relinf.ds[,1])){
  relinf.ds$median[i]<-median(as.numeric(relinf.ds[i,2:(pab.reps+1)]))
}
ord.ri<-relinf.ds[order(-relinf.ds$median),]
rim<-t(ord.ri[1:top.vars,2:101])

par(mar=c(15,5,2,1))
par(mfrow=c(1,1))
require(RColorBrewer)
cv<-brewer.pal(5,'Set1')
col.vec<-ifelse(ord.ri$type=='Distance',cv[1],
                ifelse(ord.ri$type=='Population',cv[2],
                       ifelse(ord.ri$type=='cover',cv[3],
                              ifelse(ord.ri$type=='frag',cv[4],cv[5]))))
boxplot(rim,names=ord.ri$varname[1:top.vars],las=2,cex.axis=1,ylim=c(0,6),cex.lab=1.25,col=col.vec,ylab='Relative importance',main='Figure 3')
legend(13,15,c('Distance','Population','Proportional cover','Fragmentation','Environmental'),cv,bty='n')
box(col="black",lwd=1)

# re-plot
relinf.ds$Var
varname<-c("Clearing (previous year) 1 km",
"Mean aspect 1 km",
"Clearing P:A (previous year) 0.1 km",
"Gain (all years) 0.5 km",
"Clearing P:A (previous year) 5 km",
"Loss P:A (previous 5 years) 5 km",
"Mean aspect 2 km",
"Loss (previous year) 0.5 km",
"SD NDVI 0.1 km", 
"Population density 2 km",
"Clearing P:A (previous year) 0.2 km",
"Clearing P:A (previous year) 0.5 km",
"Cover P:A (previous year) 0.1 km",
"Loss P:A (previous 5 years) 4 km",
"Loss P:A (previous year) 5 km", 
"Mean slope 0.5 km")

type<-c('cov',
        'env',
        'frag',
        'cov',
        'frag',
        'frag',
        'env',
        'cov',
        'env',
        'pop',
        'frag',
        'frag',
        'frag',
        'frag',
        'frag',
        'env')

relinf.ds$median<-NA
for(i in 1:length(relinf.ds[,1])){
  relinf.ds$median[i]<-median(as.numeric(relinf.ds[i,2:(pab.reps+1)]))
}
ord.ri<-relinf.ds[order(-relinf.ds$median),]
rim<-t(ord.ri[1:top.vars,2:101])

par(mar=c(20,5,2,1))
par(mfrow=c(1,1))
require(RColorBrewer)
cv<-brewer.pal(4,'Set1')
col.vec<- ifelse(type=='pop',cv[1],
                       ifelse(type=='cov',cv[2],
                              ifelse(type=='frag',cv[4],cv[3])))
boxplot(rim,names=varname,las=2,cex.axis=1,ylim=c(0,4.1),cex.lab=1.25,col=col.vec,ylab='Relative importance',main='',
        staplewex=0,cex=0)
legend(11,4.2,c('Proportional cover','Environmental','Fragmentation','Population'),c(cv[2],cv[3],cv[4],cv[1]),bty='n')
box(col="black",lwd=2)

# Maringal effect curves
par(mfrow=c(4,4))
par(mar=c(4,3,2,1))
for(k in 1:16){
  
  xvals<-matrix(nrow=100,ncol=pab.reps)
  for(i in 1:pab.reps){
    xvals[,i]<-plot(model_list[[i]],which(relinf.ds$Var==ord.ri$Var[k]),return.grid=T)[,1]
  }
  yvals<-matrix(nrow=100,ncol=pab.reps)
  for(i in 1:pab.reps){
    yvals[,i]<-plot(model_list[[i]],which(relinf.ds$Var==ord.ri$Var[k]),return.grid=T)[,2]
  }
  
  # plot all marginal effect curves for var with highest mean rel inf
  median.y<-rep(NA,length(yvals[,1]))
  mean.y<-rep(NA,length(yvals[,1]))
  upp.y<-rep(NA,length(yvals[,1]))
  lwr.y<-rep(NA,length(yvals[,1]))
  mean.x<-rep(NA,length(yvals[,1]))
  
  for(i in 1:length(yvals[,1])){
    median.y[i]<-median(yvals[i,])
    mean.y[i]<-mean(yvals[i,])
    mean.x[i]<-mean(xvals[i,])
    lwr.y[i]<-quantile(yvals[i,], probs = 0.025)
    upp.y[i]<-quantile(yvals[i,], probs = 0.975)
  }
  
  plot(mean.x,mean.y,lwd=4,col="white",type="l",ylim=c(min(lwr.y),max(upp.y)),
       ylab="",las=1,cex.axis=0.75,cex.lab=0.75,xlab='',main='')
  mtext(varname[k],side=1,line=2.5,cex=0.65)
  box(col="black",lwd=2)
  polygon(c(mean.x,rev(mean.x)),c(lwr.y,rev(upp.y)),col=col.vec[k],border=NA)
  lines(mean.x,mean.y,lwd=3)
}

mpred<-matrix(NA,length(dsl[[1]][,1]),pab.reps)
xcoords<-matrix(NA,length(dsl[[1]][,1]),pab.reps)
ycoords<-matrix(NA,length(dsl[[1]][,1]),pab.reps)
obs<-matrix(NA,length(dsl[[1]][,1]),pab.reps)

for(i in 1:pab.reps){
  mpred[,i] <- predict(model_list[[i]], dsl[[i]], n.trees=model_list[[i]]$gbm.call$best.trees, type="response")
  xcoords[,i] <- dsl[[i]]$coords.x1
  ycoords[,i] <- dsl[[i]]$coords.x2
  obs[,i] <- dsl[[i]]$occ
}

hist(mpred)
hist(obs)
resi<-sqrt((obs-mpred)^2)
hist(resi)
range(resi)

require(rgdal)
admin2<-readOGR(".","MYS_adm2")
kud<-subset(admin2,admin2$NAME_2=="Kudat" | admin2$NAME_2=="Kota Marudu" | admin2$NAME_2=="Kota Belud" | admin2$NAME_2=="Pitas")

par(mfrow=c(1,1))
par(mar=c(5,5,3,2))

# Plot predictions
col.vec<-rev(heat.colors(max(mpred)*100))
plot(xcoords,ycoords,col="white",xlim=c(116.5,117.3),ylim=c(6.25,7.4),xlab="Longitude",ylab="Latitude",las=1,
     main='a)')
plot(admin2,xlim=c(116.5,117.3),ylim=c(6.25,7.4),las=1,col="grey90",add=T,lwd=0.5)
points(xcoords,ycoords,col=col.vec[round(mpred*100)],pch=16,cex=0.6)
box(lwd=1.5)

arrows(117.1,6.3,117.2815,6.3,length=0,lwd=5,col="grey50")
text(117.195,6.27,label="20 km")
text(116.7,7.375,"Probability of occurence")
legend(116.5,7.35,"0.10",pch=16,cex=1,col=col.vec[10],bty="n")
legend(116.5,7.3,"0.30",pch=16,cex=1,col=col.vec[30],bty="n")
legend(116.5,7.25,"0.50",pch=16,cex=1,col=col.vec[50],bty="n")
legend(116.5,7.2,"0.70",pch=16,cex=1,col=col.vec[70],bty="n")
legend(116.5,7.15,"0.90",pch=16,cex=1,col=col.vec[90],bty="n")
length(newd$coords.x1[newd$loc_typ=="case"])
length(newd$coords.x1[newd$loc_typ=="enum"])

# Plot prediction accuracy
col.vec<-(topo.colors(max(resi)*100))
plot(xcoords,ycoords,col="white",xlim=c(116.5,117.3),ylim=c(6.25,7.4),xlab="Longitude",ylab="Latitude",las=1,
     main='b)')
plot(admin2,xlim=c(116.5,117.3),ylim=c(6.25,7.4),las=1,col="grey90",add=T,lwd=0.5)
points(xcoords,ycoords,col=col.vec[round(resi*100)],pch=16,cex=0.6)
box(lwd=1.5)

arrows(117.1,6.3,117.2815,6.3,length=0,lwd=5,col="grey50")
text(117.195,6.27,label="20 km")
text(116.6,7.375,"Prediction error")
legend(116.5,7.35,"0.01",pch=16,cex=1,col=col.vec[1],bty="n")
legend(116.5,7.3,"0.10",pch=16,cex=1,col=col.vec[10],bty="n")
legend(116.5,7.25,"0.20",pch=16,cex=1,col=col.vec[20],bty="n")
legend(116.5,7.2,"0.30",pch=16,cex=1,col=col.vec[30],bty="n")
legend(116.5,7.15,"0.40",pch=16,cex=1,col=col.vec[40],bty="n")
legend(116.5,7.1,"0.50",pch=16,cex=1,col=col.vec[50],bty="n")
length(newd$coords.x1[newd$loc_typ=="case"])
length(newd$coords.x1[newd$loc_typ=="enum"])