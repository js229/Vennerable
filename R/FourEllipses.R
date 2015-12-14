#warning("Entering FourEllipses")

#######################
# Four squares



compute.E4 <- function(V,doWeights=FALSE,s=.25,dx=0.2) {
	if (doWeights) { warning("Cant do a weighted E4") }
	if (NumberOfSets(V) != 4) { stop("fournotfour")}
	env <- new.env()
	loaded <- data(VennDiagrams,package="Vennerable",envir=env)
	stopifnot("VennDiagrams" %in% loaded)
	VennDiagrams <- get("VennDiagrams",envir=env)
	type <- "ellipses"
	if (! type  %in% names(VennDiagrams)) {
		stop(sprintf("No diagram of type %s cached\n",type))
	}
	n <- 4
	if (length(VennDiagrams[[type]])< 4 ) {
		stop(sprintf("No diagram of type %s and size %d cached (max %d)\n",type,n,length(VennDiagrams[[type]])))
	}
	TM <- VennDiagrams[[type]][[n]]

	VD <- new("VennDrawing",TM,V)
	SetLabels <- .default.SetLabelPositions(VD)
	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
	VD <- VennSetFaceLabels(VD,FaceLabels)
	VD <- .square.universe(VD,doWeights)
	VD

}	

make.E4 <- function(dx=0.02) {

	phi <- 0.8; dex <- 1.5;dey <- 2.5; a<- 7.6; e<- 0.9
	x0 <- c( -0.9, -5.0)
	VE <- list()
	VE[[1]] <- newTissueFromEllipse (x0+c(0,0),-phi ,e,-a,Set=1,dx=dx)
	VE[[2]] <- newTissueFromEllipse (x0+c(dex,0),phi ,e,a,Set=2,dx=dx)
	VE[[3]] <- newTissueFromEllipse (x0+c(-dey,dey),-phi ,e,-a,Set=3,dx=dx)
	VE[[4]] <- newTissueFromEllipse (x0+c(dex+dey,dey),phi ,e,a,Set=4,dx=dx)

	TM <- VE[[1]]
	for (ix in 2:4) {
		TM <- addSetToDrawing(TM,VE[[ix]],set2Name=paste("Set",ix,sep=""))
	}
	TM
}	
