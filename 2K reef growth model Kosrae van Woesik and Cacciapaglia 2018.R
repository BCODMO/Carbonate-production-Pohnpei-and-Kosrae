
# R van Woesik and C Cacciapaglia 2018
# 
# Keeping up with sea-level rise: carbonate production rates in Palau and Yap, western Pacific Ocean
# 
# Code to model reef growth on Palauan reefs 

#need to evalueate galaxea fascicularis


library(sp)
library(readr)
#summarize transect data
#setwd("C:/Users/Chris/Desktop/van Woesik/Palau2017/ACCRETION OF REEFS/Accretion Model")
setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae")

rtd<-reef.transectdata<-read_csv("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Kosrae/Kosrae corals 2018 July 30 Finalcc.csv")



library(plyr)
sumdata <- ddply(rtd, c("Site","Species"), summarise,
                 N    = length(Diff),
                 sum = sum(Diff),
                 mean   = mean(Diff),
                 sd   = sd(Diff),
                 planar.proportion=sum(Diff)/6000)

coraldata<-read.csv("E:/ACCRETION OF REEFS/Accretion Model/coraldensitypalau Nov 2018.csv")#coral names, density measured from palau samples, weight/displacement method g/ml or g/cm^3
morphologydata2<-read.csv("E:/ACCRETION OF REEFS/Accretion Model/species model data RvW 2.csv")
morphologydata<-morphologydata2[,-1]
#convert rugostiy-class to rogosity correction values

Morphs<-unique(as.character(coraldata$Morphology))[-11]
rugosfac<-c(.1,0.138,.4,.348,1,.333,.2,1,.5,1)
Morphs.rugfac<-data.frame(Morphs,rugosfac)

extensionrates<- read_csv("E:/ACCRETION OF REEFS/Accretion Model/edited extension rates Dec 2018.csv")
extensionrates<-as.data.frame(extensionrates)

sumdata$Species[which(sumdata$Species=='Goniopora')]<-'Goniopora sp.'
sumdata$Species[which(sumdata$Species=='Psammocora')]<-'Psammocora sp.'
sumdata$Species[which(sumdata$Species=='Euphyllia')]<-'Euphyllia sp.'
sumdata$Species[which(sumdata$Species=='Fungia')]<-'Fungia sp.'
sumdata$Species<-gsub('Dipsastraea','Dipsastrea',sumdata$Species)
#sumdata$Species<-gsub('Dipsastrea','Favia',sumdata$Species)
#sumdata$Species<-gsub('Astrea','Montastrea',sumdata$Species,ignore.case = FALSE)

rtds<-sumdata

#setup unique list to find corals
tlist<-unique(rtds$Species) #list of unique objects in transects, remove non-corals
slist<-tlist[-(grep('algae',tlist))]
slist<-slist[(grep(" ",slist))]
slist<-slist[-(grep('sponge',slist,ignore.case = TRUE))]
slist<-slist[-(grep('Soft',slist,ignore.case = TRUE))]
slist<-slist[-(grep('Sand',slist,ignore.case = TRUE))]
slist<-slist[-(grep('Dead',slist,ignore.case = TRUE))]

tlist[which(tlist%in%slist==F)] #removed

rtds<-sumdata

spplist<-sort(slist)

sites<-unique(sumdata$Site)
sites

# Set up to run through transects i to N (
# for (i in 1:length(Transect)){
#   }
GP<-numeric()
BFj<-numeric()
BUj<-numeric()
pfspeceros<-matrix(ncol=21,nrow=24)
for (j in 1:(length(sites))){
  
  Dlist<-NA;extensionratematch<-NA;Densitymatch<-NA;morphologies<-NA;Lklist<-NA
  #k is the coral species
  #m is the parrotfish species
  #n is the macroeroder
  #o is the urchin species
  
  rtdsj<-subset(rtds,Site==sites[j])#subset to focus on site level
  
  
  #Equation 1 coral Carbonate Production
  #CP[k]= (Rk*(xk)) * (Dk*Lk) *10
  # CPk is the rate of carbonate production by calcifying species k (kg m-2 y-1)
  #xk is the planar proportion cover of coral species k, per site
  #Rk is the rugosity correction factor based on coral morphology
  #Dk is the skeletal density of coral species k (g cm-3)
  #Lk is the linear-extension rate of coral species k (cm y-1)
  
  ######################## START CORAL SPECIES LOOP k FOR SITE j
  CP<-numeric()
  for(k in 1:length(spplist)){
    Rk<-NA;xk<-NA;Dk<-NA;Lk<-NA;prin<-NA;erm<-NA;dm<-NA;types<-NA #reset values
    
    #match species list with subset, if species doesn't match, apply zero
    if(sum(rtdsj$Species==spplist[k])==0){
      CP[k]=0
      G_s<-unlist(strsplit(spplist[k],' '))
    }
    if(sum(rtdsj$Species==spplist[k])>=1){ #if the species exists in the transect...
      
      G_s<-unlist(strsplit(spplist[k],' ')) #separate genus and species from the spplist[i]
      #print(G_s,justify='right')
      
      
      ############## RUGOSTIY FACTOR, Rk
      #unique(extensionrates$`growth form`) #check the rugosity correction factor for the different morphs
      #Rk<-coraldata[which(coraldata[,1]==spplist[k]),4];colnames(coraldata)[4] #rugosity
      #Rk<-as.numeric(Rk)
      Rk<-1
      
      types<-NA
      if(sum(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2],na.rm=T)>0){
        types<-unique((extensionrates[which(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2]),6]))
        types<-types[!is.na(types)]
        types<-types[1]
      }
      
      # if(sum(coraldata$genus==G_s[1]&coraldata$species==G_s[2],na.rm=T)>0){
      #   types<-as.character(coraldata[which(coraldata$genus==G_s[1]&coraldata$species==G_s[2]),7])
      # }
      
      # if(sum(morphologydata$spplist==spplist[k])>0){
      #   types<-as.character(morphologydata[which(morphologydata$spplist==spplist[k]),8])
      # }
      #print(paste('morphology =',types))
      if(sum(which(Morphs.rugfac==types))>0){
        Rk<-Morphs.rugfac[which(Morphs.rugfac==types),2]}
      
      #################PLANAR COVER, xk  
      #planar proportion is the 7th column. will be perfect match
      xk<-rtdsj[which(rtdsj$Species==spplist[k]),7]
      
      
      ################# DENSITY, Dk  (g cm-3)
      #density of coral species, must match genus+species, or genus average| mean taken between two datasets, Pratchett et al., and measured density from Palau corals
      Dktally<-c() #reset blank
      #crosscheck with both datasets
      ####### species specific data
      # if(sum(coraldata$genus==G_s[1]&coraldata$species==G_s[2],na.rm=T)>0){ #match in Genus&species
      #   Dktally[1]<-coraldata[which(coraldata$genus==G_s[1]&coraldata$species==G_s[2]),6]
      #   dm<-'species match'
      # }
      if(sum(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2],na.rm=T)>0){
        Dktally<-c(Dktally,mean(extensionrates[which(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2]),7],na.rm=T))
        dm<-'species match'
      }
      # if(sum(morphologydata$spplist==spplist[k])>0){ #match in Genus&species
      #   Dktally<-c(Dktally,morphologydata[which(morphologydata$spplist==spplist[k]),9])
      #   dm<-'species match'
      # }
      
      ####### if no species specific data, try genus for each set:
      if(sum(Dktally,na.rm=T)<.1){
        #Dktally<-(coraldata[which(coraldata$genus==G_s[1]),6])
        Dktally<-c(Dktally,extensionrates[which(extensionrates$genus==G_s[1]),7])
        dm<-'genus match'
      }
      
      
      Dktally<-unlist(Dktally)
      Dk<-mean(Dktally,na.rm=T)
      
      
      ###################### linear-extension rate, Lk
      ####### species specific data
      if(sum(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2],na.rm=T)>0){
        Lk<-mean(extensionrates[which(extensionrates$genus==G_s[1]&extensionrates$species==G_s[2]),5],na.rm=T)
        Lk<-Lk/10 #convert from mm to cm
        erm<-'species match'
      }
      
      # if(sum(morphologydata$spplist==spplist[k])>0){ #match in Genus&species
      #   Lk<-c(Lk,morphologydata[which(morphologydata$spplist==spplist[k]),10])
      #   erm<-'species match'
      # }
      
      ####### if no species specific data, try genus for each set:
      if(sum(Lk,na.rm=T)==0){
        Lk<-c(Lk,mean(extensionrates[which(extensionrates$genus==G_s[1]),5],na.rm=T))
        Lk<-Lk/10 #convert from mm to cm
        erm<-'genus match'
      } 
      
      
      Lk<-mean(Lk,na.rm=T)
      
      
    }
    
    
    ######## TOTALLING CORAL CARBONATE PRODUCTION    
    #print('   ',quote = FALSE)
    CP[k]= (Rk*(xk)) * (Dk*Lk) *10 ### run the equation 
    Dlist[k]<-Dk
    Lklist[k]<-mean(Lk,na.rm=T)
    extensionratematch[k]<-erm
    Densitymatch[k]<-dm
    morphologies[k]<-types[1] #need to combine more and find out the most common morph? changing morph?
  } 
  ###########################END CORAL SPECIES LOOP k for site j
  
  cbinefull<-cbind(spplist,CP,Dlist,Lklist,extensionratematch,Densitymatch,morphologies)
  cbine<-cbinefull[!is.na(cbinefull[,2]),]
  
  
  #barplot(as.numeric(cbine[,2]),names.arg = cbine[,1],las=2) #view coral carbonate production
  
  #CP[k]= (Rk*(xk)) * (Dk*Lk) *10
  # CPk is the rate of carbonate production by calcifying species k (kg m-2 y-1)
  #xk is the planar proportion cover of coral species k, per site
  #Rk is the rugosity correction factor based on coral morphology
  #Dk is the skeletal density of coral species k (g cm-3)
  #Lk is the linear-extension rate of coral species k (cm y-1)
  
  # Mass flux, 1 gram/second/square centimeter  =  10 kilogram/second/square meter,
  # x10 converts g cm y-1 to kg m-2 y-1 
  
  
  
  ######################## CORALINE ALGAE PRODUCTION + halimeda
  CPt<-sum(CP,na.rm=T)#print(c(CPt,'coral carbonate production kg m-2 y-1'))
  #Coralline algae
  CA<-rtdsj[which(rtdsj$Species=='CA'),7];print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=='Coralline algae'),7]);print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=="Coralline algae Porolithon" ),7]);print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=="Coralline algae Lithothamnion" ),7]);print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=="Macroalgae Amphiroa" ),7]);print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=="Halimeda" ),7]);print(CA)
  CA<-c(CA,rtdsj[which(rtdsj$Species=="Macroalgae Halimeda" ),7]);print(CA)
  
  CA<-sum(CA,na.rm=T)
  CAP= 0.018*(CA)*10 #2
  #0.018 is the mean rate of carbonate production by coralline algae j for the
  #Caribbean in g cm y-1 (Perry et al. 2012) #Check out equation??
  #check Chisholm 2000 for extra data. states between 1.5 to 10.3 kg caco3 yr for 100%cover... australia. we are 10x lower in caribbean
  #CA = coralline algae planar cover
  #CAP = carbonate production of coralline algae
  cbineCA<-rbind(cbine,c('CA',CAP,rep(NA,5)))
  pcbine<-cbineCA
  pcbine[,3]<-round(as.numeric(cbineCA[,3]),2)
  pcbine[,2]<-round(as.numeric(cbineCA[,2]),2)
  par(mar=c(10,3,1,1))
  
  
  #save cbinCA to plot all
  if(j==1){cbinesave<-cbineCA}
  if(j!=1){cbinesave<-rbind(cbinesave,cbineCA)}
  
  
  try(barplot(as.numeric(cbineCA[,2]),names.arg = cbineCA[,1],las=2))
  
  
  print('#######################################################################')
  print(c('site',j,sites[j]))
  print(pcbine,quote=F)
  
  #GROSS CARBONATE PRODUCTION (GPi) for all organisms per site (kg m-2 y-1)
  GP[j]= sum(CPt, CAP)
  print(paste(GP[j],'gross carbonate production'))
  
  
  
  #########################################
  library(readxl)
  #Bioerosion - parrotfish, sea urchins, endolithic macroborers, endolithic microborers
  ftd<- read_excel("fishes Kosrae 2018 july 12 Final.xlsx") #will change directory name based on island later
  ftd<-as.data.frame(ftd)
  ftd$size<-as.numeric(ftd$size)
  pfishdata<-read.csv("E:/ACCRETION OF REEFS/Accretion Model/pfishdata aug29 2017.csv") #.csv with fish species, volume of bite, sp, br
  #ftd<-fish.transectdata<-read.csv("pfish transect data.csv") #fish type and abundance at each site
  #Parrotfishes
  fspplist<-unique(ftd[,'species'])
  ftd$size<-as.numeric(ftd$size)
  
  fsumdata <- ddply(ftd, c('site',"species"), summarise,
                    N    = length(size),
                    sum = sum(size),
                    mean   = mean(size),
                    sd   = sd(size))
  fsumdataj<-subset(fsumdata,site==j)
  
  ftdj<-subset(ftd,sitenum==j)
  Bf<-numeric()
  for (n in 1:length(fspplist)){
    G_sF<-fspplist[n] #G_sF<-unlist(strsplit(fspplist[n],' '))
    fspplist[n]
    
    fsizes<-ftdj[which(ftdj$species==fspplist[n]),'size'] #sizes of species n
    sp<-1/(1+exp(-(-2.46142+0.08864*fsizes)))#scar proportion - based on size of parrotfish (Bonaldo and Bellwood, 2008)+ong and holland 2010
    
    
    #biterate
    reeftime<-9 #assume parrotfish spend 9hrs a day biting coral (Bonaldo and Bellwood, 2008) #hist of biterates through time. it actually is variable throughout day peaking in afternoon... could make a distribution...
    
    brc<-pfishdata[which(pfishdata$Species==G_sF),5] #bite rate constant  (Species specific from Pete Mumby)
    br<-60*abs(((4.31+brc-0.355)-(0.045*reeftime*fsizes))) #biterate based on brc and time, and size of fishes
    
    v<-exp(1.3172+.0624*fsizes)/1000 #Ong and Holland 2010
    
    ###find site specific coral density #RELATE TO FLOW RATE FOR SPECIES DENSITY, WILL INCORPORATE INTO CORALDATASHEET, or convert here, take base rate, and multiply by flow factor. using insitu density measurements for each species!!!
    
    corals.j<-rtdsj[match(cbine[,1],rtdsj[,2]),]
    totcor<-sum(corals.j$sum)
    corals.j<-cbind(corals.j,as.numeric(cbine[,3]));colnames(corals.j)[8]<-'density'
    D<-sum((corals.j$sum*corals.j$density/totcor),na.rm=T) #average density of entire site weighted by cover
    
    Bfn.all= v * sp * br * D * 365 * 0.001 #take the erosion of each individual fish in entire site
    Bf[n]= sum(Bfn.all) #find average erosion of species in all six transects
  }
  #Bfn rate of bioerosion by an individual species and size class n (g ind-1 y-1)
  # v is volume of single bite (cm3) --- size dependent
  #sp is proportion of bites leaving scars
  #br is the daily bite rate (bites d-1)
  #D is the average density of coral skeleton for the region (~ 1.67 g cm-3) #better
  #0.001 is a conversion factor converting g ind-1 yr-1 to kg ind-1 yr-1.
  
  # Calculate the bioerosion according to abundance of each fish species (Af1)
  
  #BF[n]=Bf[n]*subset(fsumdataj,species==fspplist[n])$sum# erosion *species 
  
  
  # and sum up all fishes per site
  hhh2<-cbind(fspplist,Bf)
  hhh2<-hhh2[match(fsumdataj$species,hhh2[,1]),]
  hhh2<-cbind(hhh2,fsumdataj$N)
  hhh2<-hhh2[hhh2[,1]!='NA',]
  colnames(hhh2)<-c('species','bioerosion kg yr-1 site-1','abundance')
  print(hhh2,quote=F)
  
  pfspeceros[j,]<-Bf
  colnames(pfspeceros)<-fspplist
  
  print(c(sum(Bf,na.rm=T)/(120*6*4),'site bioerosion parrotfish kg m-2 yr-1')) #total site erosion kg m-2 yr-1
  
  BFj[j]<-(sum(Bf,na.rm=T)/(120*6*4)) #120m transects, 4m wide ,6 transects per site
  
  
  
  
  ####################### URCHIN BIOEROSION
  utd<-urchin.transectdata<-read.csv("urchins Kosrae 2018 july 11 Final.csv")
  utdj<-subset(utd,site==sites[j])
  uspplist<-as.character(unique(utd$species))
  
  #Bioerosion by sea urchins
  #Urchin erosion rates, with x measured in mm, which is the size of the urchin tests
  #Calculate the erosion rate (kg urchin-1 year-1) for each
  #individual urchin using one of the following equations:
  #where x is the urchins test size
  
  utdjd<-subset(utdj,species=='Diadema')
  xd <- utdjd$diameter.cm
  BD = (0.000001*xd^3.4192)*0.365 *0.57
  #plot(x,BD)
  BDt<-sum(BD)
  #From Januchowski-Hartley et al 2017
  #Diadematidae bioerosion (BD) = (0.000001*x3.4192)*0.365
  ######
  utdje<-subset(utdj,species=='Echinometra')
  xe<- utdje$diameter.cm
  BE = (0.0004*xe^1.9786)*0.365 *0.57
  #lines(x,BE,col='blue')
  BEt<-sum(BE)
  #From Januchowski-Hartley et al 2017
  #Echinometra mathaei bioerosion (BE ) = (0.0004*x1.9786)*0.365
  
  #General equation for all other bioeroding species (BG) = (0.0001*x2.323)*0.365
  utdjg<-subset(utdj,species=='other')
  xg<-utdjg$diameter.cm
  BG= (0.0001*xg^2.323)*0.365 * 0.57
  #lines(x,BG,lwd=2)
  BGt<-sum(BG)
  #From Januchowski-Hartley et al 2017
  
  #To calculate bioerosion by urchins in kg m-2 year-1, we summed the erosion rates of all individual urchins within each
  #transect (UE), and divided by the surface areas within each transect.
  #0.57 is a correction factor to account for the proporion of sediment re-ingested during grazing (Hunter 1977)
  # conversion 0.001 (g ind y-1 to kg ind y-1) 
  
  # Calculate the bioerosion according to abundance and size of each urchin
  
  uspplist2<-c("Echinometra","Diadema","other")
  uerosion<-c(BEt, BDt, BGt)
  abundu<-c(nrow(utdje),nrow(utdjd),nrow(utdjg))
  hhh3<-(cbind(uspplist2,uerosion,abundu))
  colnames(hhh3)<-c("species",'erosion kg yr-1','abundance')
  print(hhh3,quote=F)
  #Sum up all urchins
  BUj[j]= sum(BDt, BEt, BGt)  # 9
  print(c(BUj[j],'total urchin site bioerosion kg yr-1'))
  
  NPj<-GP[j]-(BFj[j]+BUj[j])
  print("   ",quote=F)
  print(paste('Net= ',round(NPj,3),' = Gross ',round(GP[j],3),' - (urchin ',round(BUj[j],3),' + pfish ',round(BFj[j],3),")",sep=''))
  
  rugv<-ddply(rtd, c("Site"), summarise,
              sum = sum(rugosity,na.rm=T),
              mean   = mean(rugosity,na.rm=T),
              sd   = sd(rugosity,na.rm=T))
  
} #END j SITE LOOP 


#Transect Rugosity, site averaged
rugv<-ddply(rtd, c("Site"), summarise,
            sum = sum(rugosity,na.rm=T),
            mean   = mean(rugosity,na.rm=T),
            sd   = sd(rugosity,na.rm=T))
r<-10/rugv$mean
r[is.na(r)]<-1

#Sedimentation 
sed<- rep(0.4,24) #kgcaco3m2 (Glynn 1997)

## negative sedimentation for southern Kosrae, Utwe Bay
sed[23]<-(-0.4)


#Bioerosion of macroborers (clinoid sponge, molluscs, polychaete worms)
macb.erosion.constant<-10 #Glynn 1997
macblist<-c('Cliona sponge',"Sponge encrusting","encrusting Sponge",'Sponge turpois','Turpios sponge',"Cliona boring sponge","Sponge")
macb<-rtds[which((rtds[,2]%in%macblist)),] #island wide sum of macro borers.
BM<-sum(macb[,'planar.proportion'])/length(unique(rtds[,'Site']))*macb.erosion.constant #planar sum /24 sites.
#within J #BM is the regional rate of macrbioerosion (for example Whitcher 2011, St John was 0.4 kg m-2 y
#regional rate based on island wide macro prevolance



#FINAL SUM OF ALL COMPONENTS
NP<-(GP*r)+sed-(BFj+BUj+BM)
#net production = Gross production*transect rugosity + sedimentation - parrotfish - urchin - macroborer erosion






palsites <- read.csv("Kosrae sites.csv" )
#palsites<-palsites[order(palsites$letter),]

P.locat<-palsites$X #c(rep('outer',8),rep('inner',6),rep('patch',10))
islandDF<-data.frame(P.locat,NP,GP*r,BFj,BUj,palsites$real.lon,palsites$real.lat)
# hist(islandDF[islandDF$P.locat=='outer',2])
# hist(islandDF[islandDF$P.locat=='inner',2])
# hist(islandDF[islandDF$P.locat=='patch',2])










#write.csv(islandDF,'Kosrae site model3.csv')##################################################################################











sumplot <- ddply(rtds, c("Species"), summarise,
                 sum = sum(planar.proportion*10,na.rm=T))
sumplot2<-sumplot[rev(order(sumplot$sum)),]
dev.off()
par(mar=c(12,4,2,0))
barplot(sumplot2$sum, ylab='Kosrae planar proportion',names.arg=sumplot2[,1],cex.names=.9,cex.lab=1.2,las=2,yaxt="n"); axis(2, las=3)
#sum(sumplot2$sum)# barplot(GP-BFj-BUj,add=T,names.arg=sites)



cbinesave<-data.frame(cbinesave)
cbinesave$CP<-as.numeric(as.character(cbinesave$CP))

cbinesum<-ddply(cbinesave, c("spplist"), summarise,
                sumcp  = sum(CP)
)

cbinesum<-cbinesum[order(cbinesum$sumcp,decreasing=T),]
namearg<-c(as.character(cbinesum[,1]))
namearg[which(cbinesum[,1]=='CA')]<-'Coralline algae'


setwd("E:/ACCRETION OF REEFS/Accretion model 2018/Accretion Pohnpei Kosrae 2018/Figures/Kosrae")
tiff('Kosrae carb prod.tif',width=4000,height=1600,res=300)
par(mar=c(10,5,2,1))
barplot(cbinesum$sumcp,names.arg=namearg,cex.names=.9,cex.lab=1.3,las=2,ylab=expression(paste('Gross carbonate production (kg',m^2,yr^-1,')')),ylim=c(0,50))

dev.off()


#top 5
cbinesum[1:5,]
sum(cbinesum[1:5,'sumcp'])
sum(cbinesum[,'sumcp'])
sum(cbinesum[1:5,'sumcp'])/sum(cbinesum[,'sumcp'])




# 
# 
# 
# 
# #Bioerosion of macroborers (clinoid sponge, molluscs, polychaete worms)
# #Depends on substrate avaiable
# #Regional rate for Palau? assume 0.4 for now
# macb.erosion.constant<-10 #Glynn 1997
# macblist<-c('Cliona sponge',"Sponge encrusting",'Sponge turpois','Turpios sponge')
# macb<-rtds[which((rtds[,2]%in%macblist)),] #island wide sum of macro borers.
# BM<-sum(macb[,'planar.proportion'])/length(unique(rtds[,'Site']))*macb.erosion.constant #planar sum /24 sites.
# #within J #BM is the regional rate of macrbioerosion (for example Whitcher 2011, St John was 0.4 kg m-2 y
# #regional rate based on island wide macro prevolance?
# 
# Y<-rtdsj[which(rtdsj$Species=='Carbonate'),'planar.proportion'] #planar prop of carbonate at site
# 
# BI[j]=BM*Y #10
# #Yi is the per cent cover of exposed substrate in zone i
# #where BIi is the rate of bioerosion by macroborers at site i (kg m-2 y-1) #dependent on the amount of 
# 
# 
# Bi[i]=BM[i]*Y[i] #11
# #where BIi is the rate of bioerosion by microborers at site i (kg m-2 y-1)
# #BM is the regional rate of microbioerosion (for example Vogel et al. 2000, in the Bahamas was set at 0.27 kg m-2 y), 0.635 kg m-2 y-1(Tribollet et al 2002 (GBR)
# #Yi is the per cent cover of exposed substrate in zone i
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Scarids
# #Bioerosion rate BR (kg.ind-1 yr-1) 
# #BR = v*sprop*br*d*365
# #Where v is bite volume (cm3), 
# #sprop is the proportion of bites leaving scars, 
# #br is bite rate (bites day-1) 
# #and d is substratum density (kg cm-3), here taken to be 1.49*10-3 kg cm-3.
# 
# # Parrotfishes, macroeroders, and microeroders
# Bi[i]= sum (BF[i], BU[i], BI[i], Bi[i]) #12
# 
# 
# #############
# #Sediment generation
# 
# Dpi= Bi*dp[i]     #13
# #Dpi is the rate at which sediment is generated through direct production
# #Bi is rtae of sediment generation by bioerosion at site i
# #dp is the ratio of directly produced to bioeroded grains determined by thin section
# 
# ##############################
# #Net carbonate production
# Bj<-BUj+BFj
# NP[i]= GP[i] - Bj 
# # Net carbonate production = gross carbonate production minus rate of bioerosion
# #(kg m-2 y-1)
# ##########################
# #Vertical erosion and accretion
# 
# E[i]=NP[i]/D[i] * 1000    
# #E is the rate of vertical erosion (mm yr-1)
# #D is the average density of coral skeleton (in Whitcher 1670 kg m-3) (after Perry et al. 2012) 
# # 1000 converts m y-1 to mm y-1
# 
# #Mean density of coral framework for positively accreting sites
# D[i]= sum(xk/Xi * dk) *1000   
# #Di is the weighted mean density of coral framework (kg m-3)
# #xk is the percenatge cover of coral species k at site i
# #Xi is the total percentage cover of corals at site i
# #dk is the density of coral species k 
# #1000 converts g cm-3 to kg m-3
# 
# # Net rate of accretion by reef framework
# F[i]= NP[i]/D[i] *1000          
# 
# #F is the rate of framework accretion at site i (mm y-1)
# 
# 
# ################################
# #Sediment in void space (18%)
# 
# #################################
# #Vertical accretion
# A[i] = F[i]+ (F[i]* sediment/frame) +(Fi[i]* void/frame)     
# #A is the rate of reef accretion at site i (mm y-1)
# #sediment is the mean per cent of linear depth within the cores that is composed of sediment
# #void is mean percentage of linear depth within cores composed of void space
# #frame is the mean percentage of linear depth within the cores composed of in-place coral framework and rubble
# # For Caribbean 
# # sediment was set at 33%
# # void was set at 18%
# #frame was set in 49%
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########################## old code, in case needed
# # difference<-rtd$cm[1]
# # for(i in 2:nrow(rtd)){
# #   di<-rtd$cm[i]-rtd$cm[i-1]
# #   #print(di)
# #   difference<-c(difference,di)
# # }
# # difference[which(difference<=0)]<-rtd$cm[which(difference<0)]
# # rtd<-cbind(rtd,difference)
