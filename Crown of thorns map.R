#Crown of thorns map

library(maptools)
library(sp)
library(rgdal)
library(plotGoogleMaps)
library(spatstat)
#library(GISTools)
library(PBSmapping)
library(raster)
library(rgeos)
####################################


setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae/Shapefiles") #where you get files
reefs<-readOGR(".",'Kosrae reefs')
savef<-"E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae/" #where you want to save the stuff
#unique(reefs$RB_ATTRIB)








################## Finding outer reef edge
#get shapefiles
r<-raster('newras_out.gri')
r<-extend(r,c(162.88,163.04,5.25,5.38))
ri<-boundaries(r,'outer')
ro<-boundaries(r,'inner')
r<-sum(ri,ro,na.rm=T)
r[r==0]<-NA

#m<-readOGR('.','Kos inout')
#ri<-mask(r,m)
#ro<-mask(r,m,inverse=T) #### OUTER REEF BOUNDARY
outerreef<-r
plot(outerreef)

preefs<-reefs
plot(preefs[preefs$RB_ATTRIB=='fringing island'&preefs$RB_DEPTH_C==1,],col='lightblue',add=T)
plot(preefs[preefs$RB_ATTRIB=='patch island'&preefs$RB_DEPTH_C==1,],col='red',add=T)
plot(preefs[preefs$RB_ATTRIB=='barrier island'&preefs$RB_DEPTH_C==1,],col='green',add=T)

plot(preefs[preefs$RB_ATTRIB=='fringing island'&preefs$RB_DEPTH_C==2,],col='blue',add=T)
#plot(preefs[preefs$RB_ATTRIB=='patch island'&preefs$RB_DEPTH_C==2,],col='orange',add=T)
plot(preefs[preefs$RB_ATTRIB=='barrier island'&preefs$RB_DEPTH_C==2,],col='forestgreen',add=T)





#ri<-readOGR(dsn='inner reefs.kml',layer='inner1')

tkml <- getKMLcoordinates(kmlfile="inner reefs.kml", ignoreAltitude=T)
#make polygon

for(i in 1:5){
  p1 = Polygon(tkml[i])
  p2 = Polygons(list(p1), ID = "inner")
  p3= SpatialPolygons(list(p2),proj4string=CRS("+init=epsg:4326"))
  assign(paste('inner',i,sep=''),p3)
}


ir<-union(inner1,inner2)
for(i in 3:5){
  ir<-union(ir,get(paste('inner',i,sep='')))
}

plot(ro)
plot(ir,add=T)



r<-rasterize(ir,outerreef)
ri<-boundaries(r,'outer')
ro<-boundaries(r,'inner')
r<-sum(ri,ro,na.rm=T)
r[r==0]<-NA
innerreef<-r



set.seed(18)
p.randomouter<- sampleRandom(outerreef,18, xy = TRUE, sp=TRUE, na.rm = TRUE)
set.seed(13)
p.randominner<- sampleRandom(innerreef, 6, xy = TRUE, sp=TRUE, na.rm = TRUE)

image(outerreef,col='lightblue')
image(innerreef,col='orange',add=T)
points(p.randomouter, pch=16, col='blue', cex=1,lwd=2)
points(p.randominner, pch=16, col='red', cex=1,lwd=2)




#Make a table:
p.randomouter@coords;p.randominner@coords;p.randompatch@coords
stratlist<-c(rep('outer',length(p.randomouter)),rep('inner',length(p.randominner)))
Psampling<-rbind(p.randomouter@coords,p.randominner@coords)
rownames(Psampling)<-stratlist


setwd(savef)  #SAVE CSV OF SITES
#write.csv(Psampling,'Kosrae sites.csv')

setwd(savef)  #READ IN CSV of SITES
sites<-read.csv('Kosrae sites.csv',header=T)
head(sites)
inner<-sites[sites$X=='inner',]
outer<-sites[sites$X=='outer',]
patch<-sites[sites$X=='patch',]




image(outerreef,col='lightblue')
image(innerreef,col='orange',add=T)

points(outer$x,outer$y,pch=LETTERS[c(as.numeric(rownames(outer)))], col='blue', cex=2,lwd=2)
points(inner$x,inner$y,pch=LETTERS[c(as.numeric(rownames(inner)))], col='red', cex=2,lwd=2)
points(outer$x,outer$y,pch=16, col='black', cex=1)
points(inner$x,inner$y,pch=16, col='black', cex=1)






################################################# 



library(sp)
library(readr)
#summarize transect data
#setwd("C:/Users/Chris/Desktop/van Woesik/Palau2017/ACCRETION OF REEFS/Accretion Model")
setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Pohnpei")

rtd<-reef.transectdata<-read_csv("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae/Kosrae corals 2018 July 11.csv")
rtd$COT<-as.numeric(rtd$COT)
rtd$Site[rtd$Site=='k']<-'K'

library(plyr)
sumdata <- ddply(rtd, c("Site"), summarise,
                 N    = sum(COT,na.rm=T),
                 COT.per.100km2   = sum(COT,na.rm=T)/6)



sit<-sites
coordinates(sit)<- ~x+y

image(outerreef,col='lightblue')
image(innerreef,col='orange',add=T)

sit@data$COT<-c(sumdata$COT.per.100km2)
points(sit,cex=sit$COT,pch=16)



##############################################
#Krige

Ppts<-sit
proj4string(Ppts)<-CRS('+proj=longlat +datum=WGS84')
names(Ppts)[1]<-'P.locat'


library(rgdal)
library(automap)
library(gstat)

PptsO<-spTransform(Ppts, CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")) #make sure zone fits!!!!
PptsO
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



raslistfull<-list()
maskedras<-list()
i<-1
  Ppts<-PptsO
 # if(i==1){Ppts<-subset(PptsO,P.locat=='outer')}
 # if(i==2){Ppts<-subset(PptsO,P.locat=='inner')}
#  if(i==3){Ppts<-subset(PptsO,P.locat=='patch')}
  
  #set up grid
  resg<-100
  extra<-5000
  grd <- expand.grid(x = seq(from = bbox(Ppts)[1]-extra, to = bbox(Ppts)[3]+extra, by = resg), y = seq(from = bbox(Ppts)[2]-extra, to = bbox(Ppts)[4]+extra, by = resg))
  #plot(grd,pch='.', cex = .4)
  #points(Ppts, pch = 16, col = "red", cex = 1)
  gridded(grd) = ~x+y
  
  exten<-t(bbox(grd))
  exten<-SpatialPoints(exten)
  proj4string(exten)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  extenD<-spTransform(exten,CRS=CRS('+proj=longlat +datum=WGS84'))
 # e<-c(162.8801,163.0401,5.249895,5.380058)
  proj4string(grd)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  
  
Panis<- vgm(psill=.3,Nug=.13 ,"Ste", range=3200, anis= c(35, .7)) #anisotrophic krige, 40% for outer reefs, at 35 degrees
#Panis<- vgm(psill=10,Nug=.5 ,"Sph", range=10000, anis= c(35, 1))} #non-anisotrophic, shorter range
  Kpannis <- krige(Ppts@data$COT ~ 1, locations=Ppts, newdata=grd, Panis) 
  
  Kpannis<-raster(Kpannis)
  Rcoords<-coordinates(Kpannis)
  Rcoords<-SpatialPoints(Rcoords)
  proj4string(Rcoords)<-CRS("+proj=merc +zone=53s +ellps=WGS84 +datum=WGS84")
  Rcoords<-(spTransform(Rcoords,CRS=CRS('+proj=longlat +datum=WGS84')))
  
  Kpannis2<-raster(nrows=Kpannis@nrows,ncols=Kpannis@ncols)
  rasterize(Rcoords,Kpannis2)
  extent(Kpannis2)<-extenD@bbox
  values(Kpannis2)<-values(Kpannis)
  
  # plot(Kpannis2)
  # plot(palaufore, border="black", axes=TRUE, las=1,lwd=3,add=T)
  # points(spTransform(Ppts,CRS('+proj=longlat +datum=WGS84')),pch=16,cex=1,col='blue')
  
  raslistfull[i]<-Kpannis
  #masking the rasters based on the shapefiels of palau
  
  arc<-colorRampPalette((c('blue','turquoise2', 'green','yellow','orange','red')))(200)
  arc2<-colorRampPalette((c('black','white')))(200)
  arc3<-rev(terrain.colors(200))
    
  Kpannis2<-resample(Kpannis2,outerreef)
    Pmsk<-mask(Kpannis2,reefs)
    #newmsk<-buffer(Pmsk,250,doEdge=TRUE)
    Pmsk1<-mask(Kpannis2,newmsk)
    maskedras[i]<-Pmsk1
    Pmsk[Pmsk<0]<-0
    #newmsk2<-buffer(innerreef,75,doEdge=TRUE)
    
    
    plot(Pmsk,col=arc3,main=expression(paste('COT per 100m'^'2',sep='')))
    Pmsk2<-mask(Kpannis2,newmsk2)
    image(Pmsk2,col=arc3,add=T,zlim=c(0,2.6))
    points(sit,cex=.5+sit$COT,pch=1)
    
    # plot(Pmsk1,maxpixels=500000,col=colorRampPalette(rev(c('blue','turquoise2', 'green','yellow','orange','red')))(200),interpolate=F,colNA='black')
    # plot(Pmsk1,maxpixels=500000,col=(rainbow(200, start=.8, end=.3)),interpolate=F)
    
    # par(bg='black')
    # image(Pmsk1,maxpixels=500000,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2),na.rm=T)),xlab='',ylab='') #, xaxs='r',yaxs = 'r',col.axis='white'
    # compassRose(134.2453,7.72,col='white')
    # box(col='white')
    # axis(1,col='white',col.axis='white')
    # axis(2,col='white',col.axis='white')
  }
  
  # if(i==2){
  #   Kpannis2<-resample(Kpannis2,innerreef)
  #   Pmsk2<-mask(Kpannis2,innerreef)
  #   
  #   maskedras[i]<-Pmsk2
  #   
  #   image(Pmsk2,add=T,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2))),colorbar=F)
  # }
  # 
  # if(i==3){
  #   Kpannis2<-resample(Kpannis2,patch_back)
  #   Pmsk<-mask(Kpannis2,patch_back)
  #   newmsk<-buffer(Pmsk,200,doEdge=TRUE)
  #   Pmsk3<-mask(Kpannis2,newmsk)
  #   maskedras[i]<-Pmsk3
  #   
  #   image(Pmsk3,add=T,col=colorRampPalette(rev(c('blue','turquoise2','green','yellow','red')))(200),zlim=c(range(values(Kpannis2))),colorbar=F)
  # }  
  
}



dev.off()
par(bg='black')
library(plotrix)
Pmsk1[Pmsk1<0]<-0
Pmsk2[Pmsk2<0]<-0
ar=range(values(Pmsk1),na.rm=T)
br=range(values(Pmsk2),na.rm=T);allr<-c(ar,br)


#cr=range(values(Pmsk3),na.rm=T);allr<-c(ar,br,cr)
##################################################################
image(Pmsk1,maxpixels=500000,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,2),xlab='',ylab='') #, xaxs='r',yaxs = 'r',col.axis='white'
compassRose(158.083,6.815095,col='white',cex=.8)
color.legend(162.9009,5.366731,162.9409,5.371731,legend=seq(0,2,length.out=5), rect.col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),cex=1,col='white')
box(col='white')
axis(1,col='white',col.axis='white')
axis(2,col='white',col.axis='white')
image(Pmsk2,add=T,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,2),colorbar=F)
#image(Pmsk3,add=T,col=colorRampPalette((c('blue','turquoise2','green','yellow','red')))(200),zlim=c(0,12),colorbar=F)
text(158.1179,7.03861,'Net carbonate accretion',col='white')
text(158.1179,7.028,expression(paste('(kg CaCO'[3],' m'^'-2',' yr'^'-1',')',sep='')),col='white')


