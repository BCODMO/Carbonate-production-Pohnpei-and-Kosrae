# R van Woesik and C Cacciapaglia 2018
# 
# Keeping up with sea-level rise: carbonate production rates in Pohnpei, western Pacific Ocean
# 
# Code to Krige reef growth throughout Pohnpeian reefs 


library(raster)
library(rgdal)
library(readr)
setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Pohnpei")

Ppts <- read_csv("pohnpei site model3.csv",col_types = cols(X1 = col_skip()))
coordinates(Ppts)<- ~palsites.real.lon+palsites.real.lat

setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Pohnpei")
setwd(paste(getwd(),"/shapefiles",sep=''))
#get shapefiles
r<-raster('newras_out.gri')
r<-extend(r,c(158.05,158.5,6.7,7.3))
ri<-boundaries(r,'outer')
ro<-boundaries(r,'inner')
r<-sum(ri,ro,na.rm=T)
r[r==0]<-NA
m<-readOGR('.','pohn inout')
ri<-mask(r,m)
ro<-mask(r,m,inverse=T) #### OUTER REEF BOUNDARY



reefs<-readOGR(".",'Pohnpei reefs')

#preefs<-reefs
# plot(preefs[preefs$RB_ATTRIB=='patch island',])
# plot(preefs[preefs$RB_ATTRIB=='barrier atoll-bank',])
# plot(preefs[preefs$RB_ATTRIB=='barrier island',])
# plot(preefs[preefs$RB_ATTRIB=='fringing island',])
# plot(preefs[preefs$RB_ATTRIB=='patch atoll-bank',])

preefs<-crop(reefs,c(158.05,158.5,6.7,7.3))
plot(preefs)
outer<-preefs[preefs$RB_ATTRIB=='barrier island',]
plot(outer,col='blue',border='blue')
patch<-preefs[preefs$RB_ATTRIB=='patch island',]
plot(patch,col='red',border='red',add=T)
inner<-preefs[preefs$RB_ATTRIB=='fringing island',] 
plot(inner,col='goldenrod',border='goldenrod',add=T)




#PATCH REEF SHAPE #converge the backreef and the patch reefs
backreef<-rasterize(outer,r)
bi<-boundaries(backreef,'outer')
#bo<-boundaries(backreef,'inner')
#b<-sum(bi,bo,na.rm=T)
b<-bi
b[b==0]<-NA
b<-mask(b,r,inverse=T)
patchr<-rasterize(patch,r)
patch_back<-sum(b,patchr,na.rm=T);patch_back[patch_back==0]<-NA;patch_back[patch_back>0]<-1 ### PATCH REEF SHAPE


### INNER REEF SHAPE
innerr<-rasterize(inner,r)
ii<-boundaries(innerr,'inner')
#bo<-boundaries(backreef,'inner')
#b<-sum(bi,bo,na.rm=T)
i<-ii
i[i==0]<-NA
i<-mask(i,r,inverse=T)
inneredge<-i### INNER REEF SHAPE

ro
inneredge
patch_back



# prjlonlat<-CRS('+proj=longlat +datum=WGS84')
# palaufore<-spTransform(palaufore,prjlonlat)
# palinner<-spTransform(palinner,prjlonlat)
# palpatch<-spTransform(palpatch,prjlonlat)
#crop extent and re-plot
#e<-raster::`extent`(c(134.162336,134.652938,7.18335,7.790720))
# palaufore<-crop(palaufore,e)
# palinner<-crop(palinner,e)
# palpatch<-crop(palpatch,e)

e<-c(158.058,158.3752,6.75,7.1)

# land<-readOGR("C:/Users/Chris/Desktop/van Woesik/Palau2017/ACCRETION OF REEFS/PALAU 2017/Shapefiles","land_no babeldaob_1")
# land<-spTransform(land,prjlonlat)
# land<-crop(land,e)
plot(ro)
plot(inneredge, add=T)
plot(patch_back,add=T)

ro<-crop(ro,e)
inneredge<-crop(inneredge,e)
patch_back<-crop(patch_back,e)
image(ro)
plot(inneredge, add=T)
plot(patch_back,add=T)


# plot(land,col='burlywood',border='burlywood',add=T)
compassRose<-function(x,y,rot=0,cex=1,cex.dir=1,llwd=1,col='black') { 
  oldcex<-par(cex=cex) 
  mheight<-strheight("M") 
  xylim<-par("usr") 
  plotdim<-par("pin") 
  xmult<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])*plotdim[2]/plotdim[1] 
  point.angles<-seq(0,7*pi/4,by=pi/4)+pi*rot/180 
  crspans<-rep(c(mheight*3,mheight/2),4) 
  xpoints<-cos(point.angles)*crspans*xmult+x 
  ypoints<-sin(point.angles)*crspans+y 
  polygon(xpoints,ypoints,lwd=llwd,border=col) 
  txtxpoints<-cos(point.angles[c(1,3,5,7)])*1.33*crspans[1]*xmult+x 
  txtypoints<-sin(point.angles[c(1,3,5,7)])*1.33*crspans[1]+y 
  text(txtxpoints,txtypoints,c("E","N","W","S"),cex=cex.dir,col=col) 
  par(oldcex) 
}
points(Ppts)
proj4string(Ppts)<-CRS('+proj=longlat +datum=WGS84')

library(rgdal)
library(automap)
library(gstat)
PptsO<-spTransform(Ppts, CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")) #make sure zone fits!!!!

PptsO@data$eros<-PptsO@data$BFj+PptsO@data$BUj #add erosion.


raslistfull<-list()
maskedras<-list()
for(i in 1:3){
  PptsO
  if(i==1){Ppts<-subset(PptsO,P.locat=='outer')}
  if(i==2){Ppts<-subset(PptsO,P.locat=='inner')}
  if(i==3){Ppts<-subset(PptsO,P.locat=='patch')}
  
  #set up grid
  resg<-100
  extra<-10000
  grd <- expand.grid(x = seq(from = bbox(PptsO)[1]-extra, to = bbox(PptsO)[3]+extra, by = resg), y = seq(from = bbox(PptsO)[2]-extra, to = bbox(PptsO)[4]+extra, by = resg))
  #plot(grd,pch='.', cex = .4)
  #points(Ppts, pch = 16, col = "red", cex = 1)
  gridded(grd) = ~x+y
  
  
  proj4string(grd)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  
  
  if(i==1){Panis<- vgm(psill=10,Nug=.5 ,"Sph", range=16000, anis= c(35, .6))} #anisotrophic krige, 40% for outer reefs, at 35 degrees
  if(i!=1){Panis<- vgm(psill=10,Nug=.5 ,"Sph", range=10000, anis= c(35, 1))} #non-anisotrophic, shorter range
  Kpannis <- krige(Ppts@data$BFj ~ 1, locations=Ppts, newdata=grd, Panis) 
  
  Kpannis<-raster(Kpannis)
  Rcoords<-coordinates(Kpannis)
  Rcoords<-SpatialPoints(Rcoords)
  proj4string(Rcoords)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  Rcoords<-(spTransform(Rcoords,CRS=CRS('+proj=longlat +datum=WGS84')))
  
  
  #########3 
  exten<-t(bbox(grd))
  exten<-SpatialPoints(exten)
  proj4string(exten)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  extenD<-spTransform(exten,CRS=CRS('+proj=longlat +datum=WGS84'))
  
  #########
  Kpannis2<-raster(nrows=Kpannis@nrows,ncols=Kpannis@ncols)
  rasterize(Rcoords,Kpannis2)
  extent(Kpannis2)<-extenD@bbox
  values(Kpannis2)<-values(Kpannis)
  
  # plot(Kpannis2)
  # plot(palaufore, border="black", axes=TRUE, las=1,lwd=3,add=T)
  # points(spTransform(Ppts,CRS('+proj=longlat +datum=WGS84')),pch=16,cex=1,col='blue')
  
  raslistfull[i]<-Kpannis
  #masking the rasters based on the shapefiels of palau
  
  if(i==1){
    Kpannis2<-resample(Kpannis2,ro)
    Pmsk<-mask(Kpannis2,ro)
    newmsk<-buffer(Pmsk,175,doEdge=TRUE)
    Pmsk1<-mask(Kpannis2,newmsk)
    maskedras[i]<-Pmsk1
    # plot(Pmsk1,maxpixels=500000,col=colorRampPalette(rev(c('blue','turquoise2', 'green','yellow','orange','red')))(200),interpolate=F,colNA='black')
    # plot(Pmsk1,maxpixels=500000,col=(rainbow(200, start=.8, end=.3)),interpolate=F)
    
    par(bg='black')
    image(Pmsk1,maxpixels=500000,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2),na.rm=T)),xlab='',ylab='') #, xaxs='r',yaxs = 'r',col.axis='white'
    compassRose(134.2453,7.72,col='white')
    box(col='white')
    axis(1,col='white',col.axis='white')
    axis(2,col='white',col.axis='white')
  }
  
  if(i==2){
    Kpannis2<-resample(Kpannis2,inneredge)
    Pmsk2<-mask(Kpannis2,inneredge)
    maskedras[i]<-Pmsk2
    
    image(Pmsk2,add=T,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2),na.rm=T)),colorbar=F)
  }
  
  if(i==3){
    Kpannis2<-resample(Kpannis2,patch_back)
    Pmsk<-mask(Kpannis2,patch_back)
    newmsk<-buffer(Pmsk,100,doEdge=TRUE)
    Pmsk3<-mask(Kpannis2,newmsk)
    maskedras[i]<-Pmsk3
    
    image(Pmsk3,add=T,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2),na.rm=T)),colorbar=F)
  }  
  
}



dev.off()
par(bg='black')
library(plotrix)
ar=range(values(Pmsk1),na.rm=T)
br=range(values(Pmsk2),na.rm=T)
cr=range(values(Pmsk3),na.rm=T);allr<-c(ar,br,cr)
##################################################################
image(Pmsk1,maxpixels=500000,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,20),xlab='',ylab='') #, xaxs='r',yaxs = 'r',col.axis='white'
compassRose(158.083,6.815095,col='white',cex=.8)
color.legend(158.0858,7.04661,158.15,7.05661,legend=seq(0,20,length.out=5), rect.col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),cex=1,col='white')
box(col='white')
axis(1,col='white',col.axis='white')
axis(2,col='white',col.axis='white')
image(Pmsk2,add=T,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,20),colorbar=F)
image(Pmsk3,add=T,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,20),colorbar=F)
text(158.1179,7.03861,'Net carbonate accretion',col='white')
text(158.1179,7.028,expression(paste('(kg CaCO'[3],' m'^'-2',' yr'^'-1',')',sep='')),col='white')

# 
# library(plotKML)
# setwd("E:/ACCRETION OF REEFS/Accretion Model/Palau/Earth plot/outer")
# Pmsk1a<-Pmsk1
# plotKML(Pmsk1a,z.lim=c(0,20))
# 
# setwd("E:/ACCRETION OF REEFS/Accretion Model/Palau/Earth plot/inner")
# Pmsk2a<-Pmsk2
# plotKML(Pmsk2a,z.lim=c(0,20))
# 
# setwd("E:/ACCRETION OF REEFS/Accretion Model/Palau/Earth plot/patch")
# Pmsk3a<-Pmsk3
# plotKML(Pmsk3a,z.lim=c(0,20))


#Save the files ofr easy access after Krige
setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Figures/Pohnpei")
writeRaster(Pmsk1,'Pohn_outer erosion Bfj.nc',overwrite=TRUE)
writeRaster(Pmsk2,'Pohn_inner erosion Bfj.nc',overwrite=TRUE)
writeRaster(Pmsk3,'Pohn_patch erosion Bfj.nc',overwrite=TRUE)
# 





######################################### Figure
#########################################
setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Figures/Pohnpei")
pouter<-raster('Pohn_outer erosion Bfj.nc')
pinner<-raster('Pohn_inner erosion Bfj.nc')
ppatch<-raster('Pohn_patch erosion Bfj.nc')


library(OpenStreetMap)
e<-c(158.050,158.3855,6.74,7.07)
map <- openmap(c(e[4],e[1]), c(e[3],e[2]),type="bing")
map$tiles[[1]]$projection@projargs<-proj4string(pouter)
map$bbox$p1<-c(e[1],e[4])
map$bbox$p2<-c(e[2],e[3])
map$tiles[[1]]$bbox$p1<-c(e[1],e[4])
map$tiles[[1]]$bbox$p2<-c(e[2],e[3])

cols=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200)

#par(bg='white')
sz<-3



##############################
##############################
##############################



setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Figures")
tiff('Pohnpei parrot fish erosion krige.png',width=541*sz,height=530*sz, res = 300)


plot(map)
par(mar=c(0,0,0,0))
library(plotrix)
ar=range(values(pinner),na.rm=T)
br=range(values(pouter),na.rm=T)
cr=range(values(ppatch),na.rm=T);allr<-c(ar,br,cr)
##################################################################
image(pinner,maxpixels=500000,col=cols,zlim=c(0,.5),xlab='',ylab='',add=T) #, xaxs='r',yaxs = 'r',col.axis='white'
compassRose(158.083,6.815095,col='white',cex=.8)
color.legend(158.0858,7.04561,158.15,7.052,legend=(seq(0,.5,length.out=3)), rect.col=cols,cex=1,col='white')
image(ppatch,add=T,col=cols,zlim=c(0,.5),colorbar=F)
image(pouter,add=T,col=cols,zlim=c(0,.5),colorbar=F)
text(158.1179,7.03961,'Net carbonate accretion',col='white')
text(158.1179,7.0291,expression(paste('(kg CaCO'[3],' m'^'-2',' yr'^'-1',')',sep='')),col='white')

box(lwd=3,col='white')
axis(1,tck = .015,col='white',col.axis='white',lwd=2)
axis(2,tck = .015,lwd=2,col='white',col.axis='white')
axis(3,tck = .015,lwd=2,col='white',col.axis='white')
axis(4,tck = .015,lwd=2,col='white',col.axis='white')
text(158.1,6.7515,"158.1°E",cex=.975,col='white')
text(158.2,6.7515,'158.2°E',cex=.975,col='white')
text(158.3,6.7515,'158.3°E',cex=.975,col='white')

text(158.07,7.0,'7.0°N',cex=.975,col='white')
text(158.07,6.9,'6.9°N',cex=.975,col='white')

dev.off()



