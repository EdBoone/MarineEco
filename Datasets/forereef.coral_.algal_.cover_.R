##########################27 October 2014#######################
####Creating aggregated data for coral and algal cover##########
###########Focusing on 1 habitat: the forereef##################
#########################Coral cover############################

library(car)
coral<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Coral cover.csv", head=T)
head(coral)

coral$Date<-as.Date(coral$Date, "%Y-%M")
coral$Date<-as.character(coral$Date)
coral$Year<-substr(coral$Date, 1,4)

###Taking all of the substrates that are not "coral" out of the dataset###
coral<-subset(coral, Substrate!="Turf"&Substrate!="Soft Coral"&
                Substrate!="Millepora"&Substrate!="Sand"&
                Substrate!="Macroalgae"&Substrate!="Crustose Coralline Algae / Bare Space")

###Adding all of the different genera of corals together within each quadrat###
coral.one<-aggregate(coral$Percent.Cover, by=list(coral$Year, coral$Site, 
                                                  coral$Habitat, coral$Transect, coral$Quadrat), FUN=sum)
names(coral.one)<-c("Year", "Site", "Habitat", "Transect", "Quadrat", "Cover")

###Taking the average of all quadrats within each transect###
coral.two<-aggregate(coral.one$Cover, by=list(coral.one$Year, coral.one$Site, 
                                                   coral.one$Habitat, coral.one$Transect), FUN=mean)
names(coral.two)<-c("Year", "Site", "Habitat", "Transect", "Cover")

###Taking the average of all transects within each site by habitat combination###
coral.three<-aggregate(coral.two$Cover, by=list(coral.two$Year,  
                                                     coral.two$Site, coral.two$Habitat), FUN=mean)
names(coral.three)<-c("Year", "Site", "Habitat", "Cover")

###Taking the average of all sites for each habitat###
coral.four<-aggregate(coral.three$Cover, by=list(coral.three$Year, coral.three$Habitat), FUN=mean)
names(coral.four)<-c("Year", "Habitat", "Cover")

###Taking the subset for just the forereef data###
FO.coral.cover<-subset(coral.four, Habitat=="Outer 10")

###Save as csv file###
write.csv(FO.coral.cover, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef coral cover 2005 to 2013.csv")

#########################Turf algal cover############################
turf<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Coral cover.csv", head=T)
head(turf)

turf$Date<-as.Date(turf$Date, "%Y-%M")
turf$Date<-as.character(turf$Date)
turf$Year<-substr(turf$Date, 1,4)

###Taking only the subset of "turf" substrates###
turf<-subset(turf, Substrate=="Turf")

###Adding all of the different genera of turfs together within each quadrat###
turf.one<-aggregate(turf$Percent.Cover, by=list(turf$Year, turf$Site, 
                                                turf$Habitat, turf$Transect, turf$Quadrat), FUN=sum)
names(turf.one)<-c("Year", "Site", "Habitat", "Transect", "Quadrat", "Cover")

###Taking the average of all quadrats within each transect###
turf.two<-aggregate(turf.one$Cover, by=list(turf.one$Year, turf.one$Site, 
                                            turf.one$Habitat, turf.one$Transect), FUN=mean)
names(turf.two)<-c("Year", "Site", "Habitat", "Transect", "Cover")

###Taking the average of all transects within each site by habitat combination###
turf.three<-aggregate(turf.two$Cover, by=list(turf.two$Year,  
                                              turf.two$Site, turf.two$Habitat), FUN=mean)
names(turf.three)<-c("Year", "Site", "Habitat", "Cover")

###Taking the average of all sites for each habitat###
turf.four<-aggregate(turf.three$Cover, by=list(turf.three$Year, turf.three$Habitat), FUN=mean)
names(turf.four)<-c("Year", "Habitat", "Cover")

###Taking the subset for just the forereef data###
FO.turf.cover<-subset(turf.four, Habitat=="Outer 10")

###Save as csv file###
write.csv(FO.turf.cover, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef turf cover 2005 to 2013.csv")

#########################Macroalgal cover############################
macro<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Coral cover.csv", head=T)
head(macro)

macro$Date<-as.Date(macro$Date, "%Y-%M")
macro$Date<-as.character(macro$Date)
macro$Year<-substr(macro$Date, 1,4)

###Taking all of the substrates that are not "macro" out of the dataset###
macro<-subset(macro, Substrate=="Macroalgae")

###Adding all of the different genera of macros together within each quadrat###
macro.one<-aggregate(macro$Percent.Cover, by=list(macro$Year, macro$Site, 
                                                  macro$Habitat, macro$Transect, macro$Quadrat), FUN=sum)
names(macro.one)<-c("Year", "Site", "Habitat", "Transect", "Quadrat", "Cover")

###Taking the average of all quadrats within each transect###
macro.two<-aggregate(macro.one$Cover, by=list(macro.one$Year, macro.one$Site, 
                                              macro.one$Habitat, macro.one$Transect), FUN=mean)
names(macro.two)<-c("Year", "Site", "Habitat", "Transect", "Cover")

###Taking the average of all transects within each site by habitat combination###
macro.three<-aggregate(macro.two$Cover, by=list(macro.two$Year,  
                                                macro.two$Site, macro.two$Habitat), FUN=mean)
names(macro.three)<-c("Year", "Site", "Habitat", "Cover")

###Taking the average of all sites for each habitat###
macro.four<-aggregate(macro.three$Cover, by=list(macro.three$Year, macro.three$Habitat), FUN=mean)
names(macro.four)<-c("Year", "Habitat", "Cover")

###Taking the subset for just the forereef data###
FO.macro.cover<-subset(macro.four, Habitat=="Outer 10")

###Save as csv file###
write.csv(FO.macro.cover, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef macro cover 2005 to 2013.csv")

#############Combining coral, turf, and macroalgal cover into 1 dataframe###
all<-rbind(FO.coral.cover, FO.turf.cover, FO.macro.cover)
all$Substrate<-c(rep("Coral", length(FO.coral.cover$Year)),
                 rep("Turf", length(FO.turf.cover$Year)),
                 rep("Macroalgae", length(FO.macro.cover$Year)))
write.csv(all, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef coral and algal cover 2005 to 2013.csv")

