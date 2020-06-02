
# R van Woesik and C Cacciapaglia
# 
#
# 
# Code to randomly select stratified study sites on on Pohnpeis reefs 


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


getf<-"E:/ACCRETION OF REEFS/POHNPEI 2018/Pohnpei" #where you get files
savef<- "E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Pohnpei" #where you want to save files (csv/figs)

setwd(paste(getf,"/Chris Pohnpei",sep=''))
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


setwd(getf)
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

#area
cellStats(mask(area(ro),ro),sum) #outer
cellStats(mask(area(inneredge),inneredge),sum) #inner
cellStats(mask(area(patch_back),patch_back),sum) #patch

#sample, and add points 
set.seed(8)
p.randomouter<- sampleRandom(ro,8, xy = TRUE, sp=TRUE, na.rm = TRUE)
set.seed(8)
p.randominner<- sampleRandom(inneredge, 6, xy = TRUE, sp=TRUE, na.rm = TRUE)
set.seed(10)
p.randompatch<- sampleRandom(patch_back,10, xy = TRUE, sp=TRUE, na.rm = TRUE)

newout<-crop(ro,bbox(outer)) #just for visuals
image(newout,col='lightblue')
image(inneredge,col='yellow',add=T)
image(patch_back,col='pink',add=T)
points(p.randomouter, pch=16, col='blue', cex=1,lwd=2)
points(p.randominner, pch=16, col='goldenrod', cex=1,lwd=2)
points(p.randompatch, pch=16, col='red', cex=1,lwd=2)


#Make a table:
p.randomouter@coords;p.randominner@coords;p.randompatch@coords
stratlist<-c(rep('outer',length(p.randomouter)),rep('inner',length(p.randominner)),rep('patch',length(p.randompatch)))
Psampling<-rbind(p.randomouter@coords,p.randominner@coords,p.randompatch@coords)
rownames(Psampling)<-stratlist


setwd(savef)  #SAVE CSV OF SITES
#write.csv(Psampling,'Pohnpei sites.csv')

setwd(savef)  #READ IN CSV of SITES
sites<-read.csv('Pohnpei sitesCC.csv',header=T)
head(sites)
inner<-sites[sites$X=='inner',]
outer<-sites[sites$X=='outer',]
patch<-sites[sites$X=='patch',]




#figure
image(newout,col='lightblue')
image(inneredge,col='yellow',add=T)
image(patch_back,col='pink',add=T)
image(ri,add=T,col='black')
points(outer$x,outer$y,pch=LETTERS[c(as.numeric(rownames(outer)))], col='blue', cex=2,lwd=2)
points(inner$x,inner$y,pch=LETTERS[c(as.numeric(rownames(inner)))], col='goldenrod', cex=2,lwd=2)
points(patch$x,patch$y,pch=LETTERS[c(as.numeric(rownames(patch)))], col='red', cex=2,lwd=2)
points(outer$x,outer$y,pch=16, col='black', cex=1)
points(inner$x,inner$y,pch=16, col='black', cex=1)
points(patch$x,patch$y,pch=16, col='black', cex=1)



#Google earth figure
library(OpenStreetMap)
e<-c(158.058,158.3749,6.75,7.06)
map <- openmap(c(e[4],e[1]), c(e[3],e[2]),type="bing")
map$tiles[[1]]$projection@projargs<-proj4string(pouter)
map$bbox$p1<-c(e[1],e[4])
map$bbox$p2<-c(e[2],e[3])
map$tiles[[1]]$bbox$p1<-c(e[1],e[4])
map$tiles[[1]]$bbox$p2<-c(e[2],e[3])
plot(map)
points(outer$x,outer$y,pch=LETTERS[c(as.numeric(rownames(outer)))], col='lightblue', cex=2,lwd=2)
points(inner$x,inner$y,pch=LETTERS[c(as.numeric(rownames(inner)))], col='yellow', cex=2,lwd=2)
points(patch$x,patch$y,pch=LETTERS[c(as.numeric(rownames(patch)))], col='orangered', cex=2,lwd=2)
points(outer$x,outer$y,pch=16, col='white', cex=1)
points(inner$x,inner$y,pch=16, col='white', cex=1)
points(patch$x,patch$y,pch=16, col='white', cex=1)

