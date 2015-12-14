#warning("Entering 03VennDrawing")

# need to integrate DrawnVenn and AWFE classes
#**************************************************************************
#
# # $Id: 01DrawnVenn.R,v 1.1 2007/10/16 10:59:52 js229 Exp $

# $Log: 01DrawnVenn.R,v $
# Revision 1.1  2007/10/16 10:59:52  js229
# Compiles properly again
#

################################################
setGeneric("PlotUniverse",function(object,gp){standardGeneric("PlotUniverse")})
#setGeneric("IntersectionMidpoints",function(object){standardGeneric("IntersectionMidpoints")}) 
#setGeneric("SetLabelPositions",function(object){standardGeneric("SetLabelPositions")}) 
setGeneric("Areas",function(object){standardGeneric("Areas")}) 
setGeneric("ComputeAreas",function(object,nintervals){standardGeneric("ComputeAreas")}) 
setGeneric("VisibleRange",function(object){standardGeneric("VisibleRange")}) 
setGeneric("VennGetUniverseRange",function(object){standardGeneric("VennGetUniverseRange")}) 
setGeneric("VennSetUniverseRange",function(object,universe){standardGeneric("VennSetUniverseRange")}) 
setGeneric("VennSetSetLabels",function(object,SetLabels){standardGeneric("VennSetSetLabels")}) 
setGeneric("VennGetSetLabels",function(object){standardGeneric("VennGetSetLabels")}) 
setGeneric("VennSetFaceLabels",function(object,FaceLabels){standardGeneric("VennSetFaceLabels")}) 
setGeneric("VennGetFaceLabels",function(object){standardGeneric("VennGetFaceLabels")}) 



setClass("VennDrawing",representation("TissueDrawing","Venn",
	universe="matrix",SetLabels="data.frame",FaceLabels="data.frame"))

setMethod("show","VennDrawing",function(object) {
	show(as(object,"Venn"))
	show(as(object,"TissueDrawing"))
	invisible(	object)
})

setMethod("VennGetUniverseRange","VennDrawing",function(object)object@universe)
setMethod("VennSetUniverseRange","VennDrawing",function(object,universe){object@universe<-universe;object})
setMethod("VisibleRange","TissueDrawing",function(object){
	dxxy <- do.call("rbind",lapply(names(object@setList),.face.toxy,type="set",drawing=object))
	apply(dxxy,2,range)
})

.default.SetLabelPositions <- function(object){
	yscale <- diff(VisibleRange(object)[,2]); smidge <- 0.01*yscale
	sxxy <- lapply(names(object@setList),.face.toxy,type="set",drawing=object)
	setNames <- VennSetNames(as(object,"Venn"))
	smaxy <- do.call(rbind,lapply(sxxy,function(x){x[which(x[,2]==max(x[,2]))[1],] }))
	VLabels <- data.frame(Label=rep("unset",nrow(smaxy)),x=NA,y=NA,hjust=I("left"),vjust=I("bottom"))
	VLabels[,2:3] <- smaxy
	VLabels$Label <- setNames 
	VLabels
}

setMethod("VennSetSetLabels","VennDrawing",function(object,SetLabels) {object@SetLabels <- SetLabels; object})
setMethod("VennGetSetLabels","VennDrawing",function(object) {object@SetLabels})
setMethod("VennSetFaceLabels","VennDrawing",function(object,FaceLabels) {object@FaceLabels <- FaceLabels; object})
setMethod("VennGetFaceLabels","VennDrawing",function(object) {object@FaceLabels})



setMethod("PlotUniverse","VennDrawing", function(object,gp) {
	if(missing(gp)) { gp <- NULL }
	uv <- VennGetUniverseRange(object)
	grid.rect(x=mean(uv[,1]),y=mean(uv[,2]),
		width=diff(uv[,1]),height=diff(uv[,2]),default.units="native",gp=gp)
	}
)

.square.universe <- function(object,doWeights=FALSE,smudge=0.05) {
	if (FALSE & doWeights) { # never attempt to weight the Dark Matter for now
		# minimal square box 
		minimal.square.universe.area <- diff(VisibleRange(object)[,1])*diff(VisibleRange(object)[,2])
		V <- as(object,"Venn")
		visible.area <- .WeightVisible(V)
		dark.matter.area <- .WeightUniverse(V) - .WeightVisible(V)
		dark.matter.scale.squared <- (dark.matter.area + visible.area)/minimal.square.universe.area
		if (is.na(dark.matter.scale.squared) | dark.matter.scale.squared  < 1 + smudge) {
	#		warning("Square box is too large for the weight")
			dark.matter.scale.squared <- 1 + smudge
		}
	} else {
		dark.matter.scale.squared <- 1.2
	}
	dark.matter.scale <- sqrt(dark.matter.scale.squared)
	delta <- apply(VisibleRange(object),2,diff) * (dark.matter.scale-1)
	universe <- VisibleRange(object)
	universe <- universe+matrix(c(-delta,delta),ncol=2,byrow=TRUE)
		
	object@universe <- universe
	object
}


CreateViewport <- function(object) {
	xData <- VennGetUniverseRange(object)[,1]
	yData <- VennGetUniverseRange(object)[,2]
	makevp.eqsc(xData,yData)
}

UpViewports <- function() {
	upViewport()
	upViewport()
	upViewport()
}

VennThemes<- function(drawing,colourAlgorithm,increasingLineWidth) {
	gpList <- list()
	if (is.null(gpList[["Face"]])) {
		gpList[["Face"]]<- FaceColours(drawing=drawing,colourAlgorithm=colourAlgorithm)
	}
	if (is.null(gpList[["FaceText"]])) {
		gpList[["FaceText"]] <- FaceTextColours(drawing=drawing,colourAlgorithm=colourAlgorithm)
	}
	if (is.null(gpList[["Set"]])) {
		gpList[["Set"]] <- SetColours(drawing=drawing,colourAlgorithm=colourAlgorithm,increasingLineWidth)
	}
	if (is.null(gpList[["SetText"]])) {
		gpList[["SetText"]] <- SetTextColours(drawing=drawing)
	}
	gpList
}

FaceColours <- function(drawing,faceNames,colourAlgorithm) {
	if(missing(faceNames)) {	
		faceNames <- .faceNames(drawing)
	}
	faceSignatures <- .faceSignatures(drawing)[faceNames]
	sigs <- setdiff(faceSignatures,"DarkMatter"); nSets <- max(nchar(sigs))
	faceSignatures[faceSignatures == "DarkMatter"] <-  paste(rep("0",nSets),collapse="")
	#nSets <-NumberOfSets(drawing)
	nFaces <- length(faceNames)
	nSignatures <- length(unique(faceSignatures))
	DarkMatterColour <- "pink"
	setcounts <- sort(sapply(faceSignatures,function(sig){
			sigs <- strsplit(sig,"")[[1]]
			setcount <- sum(as.numeric(sigs))
			setcount
	}))
	if (missing(colourAlgorithm)) {
		if ( nSignatures >12) {
			colourAlgorithm <- "signature"
		} else {
			colourAlgorithm <- "sequential"
		}
	}
	if (colourAlgorithm=="signature") {
		countmax <- max(setcounts)
		fillcols <- c(DarkMatterColour ,brewer.pal(countmax,'YlOrRd'))
		setcolours <-fillcols[1+setcounts]; names(setcolours) <- names(setcounts)
	} else if (colourAlgorithm=="sequential"){
		fillcols <- brewer.pal(12,"Set3")
		if (nSignatures > length(fillcols)) {
			fillcols <- rep(fillcols,times=1+nSignatures/length(fillcols))
		}  
		setcolours <- c(DarkMatterColour,fillcols)[1:length(setcounts)]
		names(setcolours ) <- names(setcounts)
	} else if (colourAlgorithm=="binary"){
		setcolours <- ifelse(setcounts %%2 == 0 , "white", "blue")
	} 

	gp <- lapply(names(faceSignatures),function(x)gpar(col=setcolours [x],fill=setcolours [x],lty=0)); 	
	names(gp) <- names(faceSignatures)
	gp
}



FaceTextColours <- function(drawing,faceNames,colourAlgorithm) {
	gp <- FaceColours(drawing=drawing,faceNames=faceNames,colourAlgorithm=colourAlgorithm)
	if (!missing(colourAlgorithm)) {
		if ( colourAlgorithm=="binary") {
			bcols <- unique(sapply(gp,function(x)x$col))
			stopifnot(length(bcols)==2)
			gp <- lapply(gp,function(agp){
				res<-agp;
				res$col<- if (res$col==bcols[1]){bcols[2]}else{bcols[1]};
				res$fill<-res$col;res})
		} else {
		gp <- lapply(gp,function(agp){res<-agp;res$col<-"black";res$fill<-res$col;res})
		}
	} else {
		gp <- lapply(gp,function(agp){res<-agp;res$col<-"black";res$fill<-res$col;res})
	}
	Nsets <- NumberOfSets(drawing)
	fontsize <- if (Nsets<=3) { 20 } else {10 }
	gp <- lapply(gp,function(agp){res<-agp;res$fontsize<-fontsize;res})
	gp
}

SetTextColours <- function(drawing) {
	gp <- SetColours(drawing=drawing)
#	gp <- lapply(gp,function(agp){res<-agp;res$col<-"black";res$fill<-res$col;res})
	Nsets <- NumberOfSets(drawing)
	fontsize <- if (Nsets<=3) { 20 } else {10 }
	gp <- lapply(gp,function(agp){res<-agp;res$fontsize<-fontsize;res})
	gp
}



SetColours <- function(drawing,colourAlgorithm,increasingLineWidth) {
	if (missing(colourAlgorithm)) { colourAlgorithm <- "sequential"}
	if (missing(increasingLineWidth)) { increasingLineWidth <- FALSE}
	nSets <-length(drawing@setList)
	if (colourAlgorithm=="binary") {
		setcolours <-rep("blue",nSets)
		names(setcolours) <- names(drawing@setList)
	} else {
		fillcols <- brewer.pal(9,'Set1')
		if (nSets > length(fillcols)) {
			fillcols <- rep(fillcols,times=1+nSets/length(fillcols))
		}
 		setcolours <-fillcols[1:nSets]; 
		names(setcolours) <- names(drawing@setList)
	}
	gpList <- lapply(names(setcolours ),function(x)gpar(col=setcolours [[x]],fill=NA,lty=1,lwd=3)); 	
	if (increasingLineWidth) {
		gpList <- lapply(1:nSets,function(x){gp <- gpList[[x]];gp$lwd <- nSets - x + 1; gp})
	}
	names(gpList ) <- names(setcolours)
	gpList 
}


PlotVennGeometry <- function(C3,gpList,show=list(FaceText="weight")) {
	show.default <- list(Universe=TRUE,Sets=TRUE,SetLabels=TRUE,
		DarkMatter=FALSE,
		Faces=TRUE,
		FaceText="weight")
	unshown <- names(show)[! names(show) %in% names(show.default)]
	if (length(unshown)>0) {
		warning(sprintf("Unknown show parameters %s",paste(unshown,collapse=",")))
	}
	dmw <- Weights(C3)[dark.matter.signature(C3)]
	if (!is.na(dmw) & dmw > 0 ) {
		show.default$DarkMatter <- TRUE
	}
	show <- show[names(show) %in% names(show.default)]

	show.default[names(show)] <- show
	gp <- VennThemes(drawing=C3)
	if (!missing(gpList)) {
		for(sname in names(gpList)) {
			gp[[sname]] <- gpList[[sname]]
		}
	}
	
 	
	CreateViewport(C3)
	
	if(show.default$Universe) {
		PlotUniverse(C3)
	}
#	if(show.default$DarkMatter) {
#		PlotDarkMatter(C3)
#	}
	if (show.default$Faces) {
		PlotFaces(C3,gp=gp[["Face"]])
	}
	if(show.default$Sets) {
		PlotSetBoundaries(C3,gp=gp[["Set"]])
	}
	if(show.default$SetLabels) {
		PlotSetLabels (C3,gp=gp[["SetText"]]) 
	}

	if (length(show.default$FaceText)>0) {
		PlotIntersectionText(C3,element.plot=show.default$FaceText,
			gp=gp[["FaceText"]],
			show.dark.matter=show.default$DarkMatter)	
	}

	UpViewports()	
}

setMethod("plot",signature(x="VennDrawing",y="missing"),function(x,y,...)PlotVennGeometry(C3=x,...))

PlotIntersectionText <- function(object,gp,element.plot="weight",show.dark.matter=TRUE) {
	if (missing(gp)) gp <- FaceTextColours(object)
	V <- as(object,"Venn")
	nSets <- NumberOfSets(V)
	VI <- VennGetFaceLabels(object);
	rownames(VI) <- VI$FaceName

	
	VI$Annotation <- ""
	if( "weight" %in% element.plot ) {
		VI$WeightText <- "" 
		weights <- Weights(V)[names(Weights(V)) %in% VI$Signature]
		for (ix in seq_along(weights)) {
			VI[VI$Signature==names(weights)[ix],"WeightText"] <- weights[ix]
		}
		VI$Annotation <- paste(VI$Annotation,VI$WeightText,sep="\n")
	}
	if ("signature" %in% element.plot |"indicator" %in% element.plot) {
		VI$Annotation <- paste(VI$Annotation,VI$Signature,sep="\n")
	}
	if ("sets" %in% element.plot ) {
		sets <- seq_len(nSets)
		setpick <- strsplit(VI$Signature,split="")
		setnums <- sapply(setpick,function(w){paste(as.character(sets[w==1]),collapse="")})
 		VI$Annotation <- paste(VI$Annotation,setnums,sep="\n")
	}
	if ("elements" %in% element.plot) {
		elements <- V@IntersectionSets
		if (is.null(elements)) {
			warning("No intersection sets elements known")
			VINames <- ""
		} else {
			VINames <- sapply(elements,paste,collapse="") # 
		}
		elements <- VINames[VI$Signature]
		VI$Annotation <- paste(VI$Annotation ,elements,sep="\n")
	}
	VI$Annotation <- sub("^\n","",VI$Annotation)
	if (!show.dark.matter) {
		VI <- VI[rownames(VI)!="DarkMatter",]
	}
	if (!"hjust" %in% colnames(VI)) { VI$hjust <- "centre" }
	if (!"vjust" %in% colnames(VI)) { VI$vjust <- "centre" }
	hj <-sapply( VI$hjust,function(EXPR){switch(EXPR,left=0,right=1,center=,centre=0.5)})
	vj <-sapply( VI$vjust,function(EXPR){switch(EXPR,top=1,bottom=0,center=,centre=0.5)})
	for (ij in 1:nrow(VI)) {
		grid.text(x=VI$x[ij],y=VI$y[ij],hjust=hj[ij],
		gp=gp[[VI$FaceName[ij]]],
		vjust=vj[ij],label=VI$Annotation[ij],default.units="native")
	}
}

.default.FaceLabelPositions <- function(object){	
	dm <-  dark.matter.signature(object)
	ilabels <- data.frame(internalPointsofFaces(as(object,"TissueDrawing")))
	colnames(ilabels) <- c("x","y")
	ilabels$FaceName <- rownames(ilabels)
	sigs <- unlist(.faceSignatures(object))
	sigdf <- data.frame(Signature=sigs,FaceName=names(sigs),stringsAsFactors=FALSE)
	df <- merge(sigdf,ilabels)
	df[df$FaceName=="DarkMatter","Signature"] <- dm
	df$hjust <- I("centre")
	df$vjust <- I("centre")
	df[df$Signature==dm,c("hjust")] <- c("right")
	df[df$Signature==dm,c("vjust")] <- c("top")
	df[df$Signature==dm,c("x","y")] <- VisibleRange(object)[2,]
	df
}

setMethod("Areas","VennDrawing",function(object) {
	areas <- faceAreas(as(object,"TissueDrawing"))
	names(areas)[names(areas)=="DarkMatter"] <- dark.matter.signature(object)

	areas
})



PlotSetLabels <- function(object,gp) {
	VLabels <- VennGetSetLabels(object)
	if(missing(gp)) gp <- SetColours(object)
	if(nrow(VLabels)==0){ warning("Can't show Set labels"); return()}
#	print(VLabels)
	# just may not be vectorised...
	hj <-sapply( VLabels$hjust,function(EXPR){switch(EXPR,left=0,right=1,center=,centre=0.5)})
	vj <-sapply( VLabels$vjust,function(EXPR){switch(EXPR,top=1,bottom=0,center=,centre=0.5)})

	for (ij in 1:nrow(VLabels)) {
		grid.text(x=VLabels$x[ij],y=VLabels$y[ij],hjust=hj[ij],
		vjust=vj[ij],gp=gp[[ij]],label=as.character(VLabels$Label[ij]),default.units="native")
	}
}


makevp.eqsc <- function(xrange,yrange) {
	# cf Fig 7.4 of Murrell R Graphics
	pushViewport(plotViewport(name="Vennmar",c(1,1,1,1)))
	pushViewport(viewport(name="Vennlay",layout=grid.layout(1,1,widths=diff(xrange),heights=diff(yrange),respect=TRUE)))
	pushViewport(viewport(name="Vennvp",layout.pos.row=1,layout.pos.col=1,xscale=xrange,yscale=yrange))
	}

PlotDarkMatter <- function(VD) {
	ur <- VennGetUniverseRange(VD)
	grey <- brewer.pal(8,"Greys")[2]
	grid.polygon(x=ur[c(1,1,2,2),1],y=ur[c(1,2,2,1),2],gp=gpar(fill=grey))
	.PlotFace.TissueDrawing(VD,"DarkMatter",gp=gpar(fill="white"),doDarkMatter=TRUE)
	}

