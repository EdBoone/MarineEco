####################################27 October 2014##########################################
######Aggregated fish biomass and abundance by functional group for the forereef habitat#####
#############################Biomass#########################################################
library(car)
library(reshape)

fish<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Fish data complete.csv", head=T)
functional.group<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Functional groups.csv", 
                           head=T)

###Merge the two datasets together
fish.func.group<-merge(fish, functional.group, by.x=c("Taxonomy"),
                       by.y=c("Species"))

##################Remove Transect 3 from all further analysis################################
###################Remove "unknown" functional group and convict tangs#######################
fish.func.group<-subset(fish.func.group, Transect<3|Transect>3)
fish.func.group<-na.omit(fish.func.group)
fish.func.group<-subset(fish.func.group, 
                               fish.func.group$Functional.group!="Unknown")
fish.func.group<-subset(fish.func.group, 
                               fish.func.group$Taxonomy!="Acanthurus triostegus")

##################Summing biomass across transects for each site#############################
biom.agg<-aggregate(fish.func.group$Biomass, by=list(fish.func.group$Taxonomy,
          fish.func.group$Functional.group, fish.func.group$Year, fish.func.group$Site, 
          fish.func.group$Habitat), FUN=sum)
names(biom.agg)<-c("Species", "Functional.group", "Year", "Site", "Habitat", "Biomass")

##################Summing biomass across species for each functional group###################
biom.agg2<-aggregate(biom.agg$Biomass, by=list(biom.agg$Functional.group, 
                                biom.agg$Year, biom.agg$Site, biom.agg$Habitat), FUN=sum)
names(biom.agg2)<-c("Functional.group", "Year", "Site", "Habitat", "Biomass")

##################Averaging biomass across sites#############################################
##########Unit for biomass is "g/750 meters squared"#########################################
biom<-aggregate(biom.agg2$Biomass, by=list(biom.agg2$Year, biom.agg2$Functional.group, 
                biom.agg2$Habitat), FUN=mean)
names(biom)<-c("Year", "Functional.group", "Habitat", "Biomass")

################Subsetting the data for just the forereef habitat###########################
FO.biom<-subset(biom, Habitat=="FO")
write.csv(FO.biom, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef herbivore biomass by functional group and year.csv")

#############################Count/abundance###############################################
#                                                                                         #
#                                                                                         #
#############################Count/abundance###############################################
library(car)
library(reshape)

fish<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Fish data complete.csv", head=T)
functional.group<-read.csv("X:/Dropbox/SAMSI/Marine Ecology/For R/Functional groups.csv", 
                           head=T)

###Merge the two datasets together
fish.func.group<-merge(fish, functional.group, by.x=c("Taxonomy"),
                       by.y=c("Species"))

##################Remove Transect 3 from all further analysis################################
###################Remove "unknown" functional group and convict tangs#######################
fish.func.group<-subset(fish.func.group, Transect<3|Transect>3)
fish.func.group<-na.omit(fish.func.group)
fish.func.group<-subset(fish.func.group, 
                        fish.func.group$Functional.group!="Unknown")
fish.func.group<-subset(fish.func.group, 
                        fish.func.group$Taxonomy!="Acanthurus triostegus")

##################Summing count across transects for each site#############################
abund.agg<-aggregate(fish.func.group$Count, by=list(fish.func.group$Taxonomy,
            fish.func.group$Functional.group, fish.func.group$Year, fish.func.group$Site, 
            fish.func.group$Habitat), FUN=sum)
names(abund.agg)<-c("Species", "Functional.group", "Year", "Site", "Habitat", "Count")

##################Summing count across species for each functional group###################
abund.agg2<-aggregate(abund.agg$Count, by=list(abund.agg$Functional.group, 
                                abund.agg$Year, abund.agg$Site, abund.agg$Habitat), FUN=sum)
names(abund.agg2)<-c("Functional.group", "Year", "Site", "Habitat", "Count")

##################Averaging count across sites#############################################
##########Unit for count is "number of individuals/750 meters squared"#####################
abund<-aggregate(abund.agg2$Count, by=list(abund.agg2$Year, abund.agg2$Functional.group, 
                                           abund.agg2$Habitat), FUN=mean)
names(abund)<-c("Year", "Functional.group", "Habitat", "Count")

################Subsetting the data for just the forereef habitat###########################
FO.abund<-subset(abund, Habitat=="FO")
write.csv(FO.abund, "X:/Dropbox/SAMSI/Marine Ecology/From R/Forereef herbivore abundance by functional group and year.csv")
