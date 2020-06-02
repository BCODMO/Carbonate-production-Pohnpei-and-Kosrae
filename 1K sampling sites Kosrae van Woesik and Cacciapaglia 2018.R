
# R van Woesik and C Cacciapaglia
# 
#
# 
# Code to randomly select stratified study sites on on Kosrae's reefs 

#Libraries that you will need
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


getf<- "E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae" #where you get files
reefs<-readOGR(".",'Kosrae reefs')
savef<-getf #where you want to save the stuff
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
plot(preefs[preefs$RB_ATTRIB=='patch island'&preefs$RB_DEPTH_C==2,],col='orange',add=T)
plot(preefs[preefs$RB_ATTRIB=='barrier island'&preefs$RB_DEPTH_C==2,],col='forestgreen',add=T)





ri<-readOGR(dsn='inner reefs.kml',layer='inner1')

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





#Google earth figure
library(OpenStreetMap)
e<-c(162.88,163.04,5.25,5.38)
map <- openmap(c(e[4],e[1]), c(e[3],e[2]),type="bing")
map$tiles[[1]]$projection@projargs<-proj4string(youter)
map$bbox$p1<-c(e[1],e[4])
map$bbox$p2<-c(e[2],e[3])
map$tiles[[1]]$bbox$p1<-c(e[1],e[4])
map$tiles[[1]]$bbox$p2<-c(e[2],e[3])

plot(map)
points(outer$x,outer$y,pch=LETTERS[c(as.numeric(rownames(outer)))], col='black', cex=2.2,lwd=2)
points(inner$x,inner$y,pch=LETTERS[c(as.numeric(rownames(inner)))], col='yellow', cex=2.2,lwd=2)
points(outer$x,outer$y,pch=LETTERS[c(as.numeric(rownames(outer)))], col='white', cex=2,lwd=2)
points(outer$x,outer$y,pch=16, col='white', cex=1)
points(inner$x,inner$y,pch=16, col='white', cex=1)



