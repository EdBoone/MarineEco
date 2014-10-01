###################################################################
##
##
##
###################################################################
setwd("~/VCU/Google Drive/MarineEco")  # Ed's Mac Directory

x1 <- read.csv( "MCR_LTER_Fish_Coral_Abundance_Change_20140605.csv",
                header = TRUE)
head(x1)

LTER1 <- x1[x1$Site == "LTER 1", ]
LTER2 <- x1[x1$Site == "LTER 2", ]
LTER3 <- x1[x1$Site == "LTER 3", ]
LTER4 <- x1[x1$Site == "LTER 4", ]
LTER5 <- x1[x1$Site == "LTER 5", ]
LTER6 <- x1[x1$Site == "LTER 6", ]

pdf("MeanLiveCoral.pdf")
plot( LTER1$Year, LTER1$Mean_live_coral, type = "l", col = "blue")
lines( LTER2$Year, LTER2$Mean_live_coral, lty = 2 , col = "red" )
lines( LTER3$Year, LTER3$Mean_live_coral, lty = 3 , col = "red" )
lines( LTER4$Year, LTER4$Mean_live_coral, lty = 4 , col = "red" )
lines( LTER5$Year, LTER5$Mean_live_coral, lty = 5 , col = "red" )
lines( LTER6$Year, LTER6$Mean_live_coral, lty = 6 , col = "red" )
dev.off()

pdf("MeanTotalCoral.pdf")
plot( LTER1$Year, LTER1$Mean_total_coral, type = "l", col = "blue")
lines( LTER2$Year, LTER2$Mean_total_coral, lty = 2 , col = "red" )
lines( LTER3$Year, LTER3$Mean_total_coral, lty = 3 , col = "red" )
lines( LTER4$Year, LTER4$Mean_total_coral, lty = 4 , col = "red" )
lines( LTER5$Year, LTER5$Mean_total_coral, lty = 5 , col = "red" )
lines( LTER6$Year, LTER6$Mean_total_coral, lty = 6 , col = "red" )
dev.off()

plot( LTER1$Year, LTER1$Mean_live_coral/LTER1$Mean_total_coral, 
      type = "l", 
      col = "blue",
      ylim = c(0,1))
lines( LTER2$Year, LTER2$Mean_live_coral/LTER2$Mean_total_coral, lty = 2 , col = "red" )
lines( LTER3$Year, LTER3$Mean_live_coral/LTER3$Mean_total_coral, lty = 3 , col = "red" )
lines( LTER4$Year, LTER4$Mean_live_coral/LTER4$Mean_total_coral, lty = 4 , col = "red" )
lines( LTER5$Year, LTER5$Mean_live_coral/LTER5$Mean_total_coral, lty = 5 , col = "red" )
lines( LTER6$Year, LTER6$Mean_live_coral/LTER6$Mean_total_coral, lty = 6 , col = "red" )

pdf("CoralDweller.pdf")
plot(x1$Mean_live_coral, x1$Coral_Dweller,
     xlab="Live Coral",
     ylab="Change in Coral Dweller",
     col = x1$Site)
dev.off()

pdf("Corrallivore.pdf")
plot(x1$Mean_live_coral, x1$Corallivore,
     xlab="Live Coral",
     ylab="Change in Corallivore",
     col = x1$Site)
dev.off()
