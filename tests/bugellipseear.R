#sfiles <- list.files(full.names=TRUE,"C:\\JonathanSwinton\\Vennerable\\pkg\\Vennerable\\R")
#for (sfile in sfiles) source(sfile)
library(Vennerable)

	phi <- 0.8; dex <- 1.7;dey <- 2.5; a<- 7.6; e<- 0.9
	x0 <- c( -0.9, -5.0)
	VE <- list()
	dx <- 0.2
	VE[[1]] <- Vennerable:::newTissueFromEllipse (x0+c(0,0),-phi ,e,-a,Set=1,dx=dx)
	VE[[2]] <- Vennerable:::newTissueFromEllipse (x0+c(dex,0),phi ,e,a,Set=2,dx=dx)
#undebug(Vennerable:::newTissueFromEllipse )
	VE[[3]] <- Vennerable:::newTissueFromEllipse (x0+c(-dey,dey),-phi ,e,-a,Set=3,dx=dx)
	VE[[4]] <- Vennerable:::newTissueFromEllipse (x0+c(dex+dey,dey),phi ,e,a,Set=4,dx=dx)

	TM <- VE[[1]]
	TM2 <- Vennerable:::addSetToDrawing(TM,VE[[2]],set2Name=paste("Set",2,sep=""))
#undebug(.find.triangle.within.face)
#undebug(.find.point.within.face)
	TM3 <- Vennerable:::addSetToDrawing(TM2,VE[[3]],set2Name=paste("Set",3,sep=""))
 