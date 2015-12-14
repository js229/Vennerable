#warning("Entering ThreeSquares")

	


######################
compute.S2 <- function(V,doWeights=TRUE,doEuler=FALSE) {
	stopifnot(NumberOfSets(V)==2)
	Weight <- Weights(V)
	if (!doWeights){
		Weight <- 1+0*Weight
		if (doEuler) {
			Weight[Weights(V)==0] <- 0
		} 
	}
	wab <- Weight["00"]
	wAb <- Weight["10"]
	waB <- Weight["01"]
	wAB <- Weight["11"]

	wA <- wAb+ wAB
	wB <- waB + wAB

	s1 <- sqrt(wA)
	s2 <- sqrt(wB)

	# squares are centred at (-+d/2,0) and have sides s1, s2
	if (wAb==0) {
		#d <- (s2-s1)/4 to have it completely inside
		d <- (s2-s1)/2
		if (!doEuler) {
			d <- (s2-s1)/2+ 0.1*s1
		}
	} else if (waB==0) {
		#d <- (s1-s2)/4
		d <- (s1-s2)/2
		if (!doEuler) {
			d <- (s1-s2)/2+ 0.1*s2
		}
	} else {
		d <- (s1+ s2)/2 -   wAB / min(s1,s2)
		if(wAB==0 & !doEuler) {
			d <- d - .1* min(s1,s2)
		}
	}

	l1 <- -d/2-s1/2; l2 <- d/2-s2/2
	r1 <- -d/2+s1/2; r2 <- d/2+s2/2

	poly.1 <- matrix(c(l1,-s1/2,l1,s1/2,r1,s1/2,r1,-s1/2),ncol=2,byrow=TRUE)
	rownames(poly.1) <- paste("s",1:4,sep="")
	poly.2 <- matrix(c(l2,-s2/2,l2,s2/2,r2,s2/2,r2,-s2/2),ncol=2,byrow=TRUE)
	rownames(poly.2) <- paste("s",2:5,sep="")
	VDP1 <- newTissueFromPolygon(points.xy=poly.1,Set=1)
	VDP2 <- newTissueFromPolygon(points.xy=poly.2,Set=2)
	TM <- addSetToDrawing (drawing1=VDP1 ,drawing2=VDP2, set2Name="Set2")
	VD <- new("VennDrawing",TM,V)
	VD <- .square.universe(VD,doWeights)
	SetLabels <- .default.SetLabelPositions(VD)
	SetLabels[1,c("x","y")] <- c(l1,s1/2)
	SetLabels[2,"hjust"] <- "right"
	SetLabels[2,c("x","y")] <- c(r2,s2/2)
	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
		if (s1 < s2) {
		dmxy <- c(l1,s1/2)
	} else if (s1>2) {
		dmxy <- c(r1,s2/2)
	} else {
		dmxy <-c (r1,s2/2 * 1.05)
	}
	FaceLabels[FaceLabels$FaceName=="DarkMatter",c("x","y")] <- dmxy
	FaceLabels[FaceLabels$FaceName=="11","x"] <- (l2+r1)/2
	FaceLabels[FaceLabels$FaceName=="10","x"] <- (l1+l2)/2
	FaceLabels[FaceLabels$FaceName=="01","x"] <- (r1+r2)/2
	VD <- VennSetFaceLabels(VD,FaceLabels)

}


################################

compute.S3 <- function(V,doWeights=TRUE) {
	stopifnot(NumberOfSets(V)==3)
	wght <- Weights(V)
	VS <- VennSignature(V)
	if ( doWeights & wght[VS=="111"]==0 ) {
		stop("Can't generate a three squares plot for a nonzero central intersection")
	}
		
	if (!doWeights) { 
		x <- 1; y <- 3
		wght <- y + 0 * wght
		wght["111"] <- wght["101"] <- wght["011"] <- wght["001"] <- x
	}

	# currently easier thant TwoSqures because we don't worry about zero weights
	IntersectionShapes <- vector("list",length(VS))
	names(IntersectionShapes) <- VS
	
	wabc <- wght[VS=="111"] 
	d1 <- 1
	wabc.width <- sqrt(wabc)/d1
	wabc.height <- sqrt(wabc)*d1
	wabc.x <- c(-wabc.width/2,wabc.width/2,wabc.width/2,-wabc.width/2)
	wabc.y <- c(-wabc.height/2,-wabc.height/2,wabc.height/2,+wabc.height/2)

	IntersectionShapes[["111"]] <- matrix(c(wabc.x,wabc.y),byrow=FALSE,ncol=2)
	
	wac <- wght[VS=="101"] 
	wac.height <- wac/wabc.width

		
	wac.x <- wabc.x
	wac.y <-  wabc.height/2+c(0,0,wac.height,wac.height)
	IntersectionShapes[["101"]] <- matrix(c(wac.x,wac.y),byrow=FALSE,ncol=2)

	wbc <- wght[VS=="011"] 
	wbc.height<- wabc.height
	wbc.width <- wbc/wbc.height
	wbc.y <- wabc.y
	wbc.x <-  wabc.width/2+c(0,wbc.width,wbc.width,0)
	IntersectionShapes[["011"]] <- matrix(c(wbc.x,wbc.y),byrow=FALSE,ncol=2)

	wab <- wght[VS=="110"]
	# (d +wabc.height)^2=wab+wabc
	#d <- 0.5 * wab/wabc.width
	# 
	d<- sqrt(wab+wabc)-wabc.height
	stopifnot(d>=0)
	wab.width <- (wab- d* wabc.width)/(wabc.height+d) + wabc.width
	wab.height <- wabc.height+d
	wab.x <-  wabc.width/2 + c(-wab.width,0,0,-wabc.width,-wabc.width,-wab.width)
	wab.y <- wabc.height/2+c(-wab.height,-wab.height,-wabc.height,-wabc.height,0,0) 
	IntersectionShapes[["110"]] <- matrix(c(wab.x,wab.y),byrow=FALSE,ncol=2)

	
	wc <- wght[VS=="001"]
	wc.height <- (wac.height+wabc.height)
	if (wc < wbc.width* wc.height - wbc ) { # not square
		dbc <- wc/(wbc.width+wc.height)
		wc.x <- wabc.width/2+  c(wbc.width,wbc.width+dbc, wbc.width+dbc,           dbc,      dbc,        0,         0,wbc.width)
		wc.y <- -wabc.height/2+c(        0,            0,wbc.height+dbc,wbc.height+dbc,wc.height,wc.height,wbc.height,wbc.height)
	} else {
		wc.width <- (wbc+wc)/wc.height # > wbc.width
		wc.x <- wabc.width/2+  c(wbc.width,wc.width,wc.width,        0,         0,wbc.width)
		wc.y <- -wabc.height/2+c(        0,                 0,          wc.height,wc.height,wbc.height,wbc.height)
	}
	IntersectionShapes[["001"]] <- matrix(c(wc.x,wc.y),byrow=FALSE,ncol=2)

	wb <- wght[VS=="010"]
	if (!doWeights) { # special is case to send 010 the other way
		wb.height <- wab.height
		wb.width <- ( wbc+wb)/wb.height
		wb.x <- wabc.width/2+c(0,wb.width,wb.width,wbc.width,wbc.width,0)
		wb.y <- wabc.height/2-wb.height+c(0,0,wb.height,wb.height,wb.height-wbc.height,wb.height-wbc.height)
	} else {
	wb.width <- wab.width + wbc.width
	wb.height <- wab.height
	if ( wb < wbc.width*(wab.height-wabc.height) ) { # need to notch 010
		# not sure this code has ever been tested, and last time I checked xa should be 1!
		xa <- 2; xb <- wab.width + wbc.width +wab.height; xc <- - wb
		db <- (-xb+sqrt(xb^2-4*xa*xc))/(2*xa)
		wb.width <- wb.width+db
		wb.height <- wb.height+db
		wb.x <- -wabc.width/2-(wab.width-wabc.width)-db+
			c(0,wab.width+2*db,wab.width+2*db,wb.width,wb.width,wab.width+db,wab.width+db,db,db,0)
		wb.y <- wabc.height/2-wb.height+
			c(0,0,wab.height-wabc.height,wab.height-wabc.height,
			wb.height-wabc.height,wb.height-wabc.height,db,db,wb.height,wb.height)
	} else {
		wb.height <- (wabc+wbc+wab+wb)/wb.width
		# that's the height giving no overlap with a
		sf <- (wb.height/wab.height) # > 1, the extent below
		db <- sqrt(sf)
		wb.height <- wb.height/db
		wb.width <- wb.width*db
		wb.x <- +wabc.width/2+wbc.width-wb.width+
			c(0,wb.width,wb.width,wb.width-wbc.width,wb.width-wbc.width,wb.width-wbc.width-wab.width,wb.width-wbc.width-wab.width,0)
		wb.y <- wabc.height/2-wb.height+
			c(0,0,wb.height-wabc.height,wb.height-wabc.height,wb.height-wab.height,wb.height-wab.height,wb.height,wb.height)
	}
	} # else (!doWeights)
	IntersectionShapes[["010"]] <- matrix(c(wb.x,wb.y),byrow=FALSE,ncol=2)

	wa <- wght[VS=="100"]

	wa.width <- wab.width 
	if ( wa < (wab.width-wabc.width)*wab.height ) {
		da <- wa/(wa.width+ wac.height)
		wa.x <- -wabc.width/2-(wab.width-wabc.width)+
			c(0,wab.width-wabc.width,wab.width-wabc.width,wab.width,wab.width,wab.width-wabc.width-da,wab.width-wabc.width-da,0)
		wa.y <- wabc.height/2+
			c(0,0,wac.height,wac.height,wac.height+da,wac.height+da,da,da)
	} else {
		wa.height <- (wac+wb)/wa.width
		wa.x <- -wabc.width/2-(wab.width-wabc.width)+
			c(0,wab.width-wabc.width,wab.width-wabc.width,wab.width,wab.width,0)
		wa.y <- wabc.height/2+
			c(0,0,wac.height,wac.height,wa.height,wa.height)
	}
	IntersectionShapes[["100"]] <- matrix(c(wa.x,wa.y),byrow=FALSE,ncol=2)

	SetShapes <- list()
	SetShapes[[1]] <- do.call(rbind,list(
		IntersectionShapes[["110"]][1:2,],IntersectionShapes[["100"]][5:nrow(IntersectionShapes[["100"]]),]
		))
	SetShapes[[2]] <- if (doWeights) {
		do.call(rbind,list(
		IntersectionShapes[["010"]][1:(nrow(IntersectionShapes[["010"]])-6),],IntersectionShapes[["011"]][3,],
			IntersectionShapes[["010"]][nrow(IntersectionShapes[["010"]]),]
		))
		} else {
		do.call(rbind,list(
		IntersectionShapes[["110"]][1,],
		IntersectionShapes[["010"]][2:3,],
		IntersectionShapes[["110"]][6,]
		))
		}
	SetShapes[[3]] <- do.call(rbind,list(
		IntersectionShapes[["111"]][1,],IntersectionShapes[["001"]][2:(nrow(IntersectionShapes[["001"]])-3),],
			IntersectionShapes[["101"]][4,]
		))
	for (ix in 1:3) {
		rownames(SetShapes[[ix]]) <- paste("c",letters[ix],1:nrow(SetShapes[[ix]]),sep="")
	}
	SetShapes[[1]] <- SetShapes[[1]][nrow(SetShapes[[1]]):1,]
	SetShapes[[2]] <- SetShapes[[2]][nrow(SetShapes[[2]]):1,]
	SetShapes[[3]] <- SetShapes[[3]][nrow(SetShapes[[3]]):1,]

	VDP1 <- newTissueFromPolygon(points.xy=SetShapes[[1]],Set=1)
	VDP2 <- newTissueFromPolygon(points.xy=SetShapes[[2]],Set=2)
	VDP3 <- newTissueFromPolygon(points.xy=SetShapes[[3]],Set=3)
	TM <- addSetToDrawing (drawing1=VDP1 ,drawing2=VDP2, set2Name="Set2")
	TM <- addSetToDrawing (drawing1=TM ,drawing2=VDP3, set2Name="Set3")
	VD <- new("VennDrawing",TM,V)
	VD <- .square.universe(VD,doWeights)
	SetLabels <- .default.SetLabelPositions(VD)
	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
	VD <- VennSetFaceLabels(VD,FaceLabels)

	VD
}



compute.S4 <- function(V,doWeights=FALSE,s=.25,likeSquares=TRUE) {
	if (doWeights) { warning("Cant do a weighted S4") }
	if (NumberOfSets(V) != 4) { stop("fournotfour")}
	
	# the nodes in an 4 set AWFE diagram see Chow Ruskey Fig 6
	top.2 <- 2+s
	px <- c(0,0,0,0, 0, 0,2,1,-1,-2,1, 1,-1,-1)
	py <- c(top.2,2,1,0,-1,-2,0,0, 0, 0,1,-1,-1, 1)
	pxy <- cbind(px,py); rownames(pxy)<- paste("p",1:14,sep="")
	

	p6p1 <- matrix(c(   0,-2-s,   -2 - 2*s,-2-s,   -2-2*s, top.2+s,0,top.2+s),ncol=2,byrow=TRUE)
	rownames(p6p1) <- paste("pa",1:nrow(p6p1),sep="")
	p1p7 <- matrix(c( 2+s,top.2, 2+s,0),ncol=2,byrow=TRUE)
	rownames(p1p7 ) <- paste("pb",1:nrow(p1p7 ),sep="")
	p10p1 <- matrix(c( -2-s,0, -2-s,top.2),ncol=2,byrow=TRUE)
	rownames(p10p1)  <- paste("pc",1:nrow(p1p7 ),sep="")

	Set <- list()
	# join all the nodes with straight lines
	Set[[1]] <- newTissueFromPolygon(Set=1,points.xy=rbind(pxy[c(1,2,3,4,5,6),],p6p1))
	Set[[2]] <- newTissueFromPolygon(Set=2,points.xy=rbind(pxy[1,,drop=FALSE],p1p7,pxy[c(7,8,4,9,10),],p10p1))
	Set[[3]] <- newTissueFromPolygon(Set=3,points.xy=pxy[c(2,11,8,12,6,13,9,14),])
	Set[[4]] <- newTissueFromPolygon(Set=4,points.xy=pxy[c(11,7,12,5,13,10,14,3),])

	diag.to.cityblock <- function(xy,first="x") {
		block.point <- if(first=="x") { c(xy[2,1],xy[1,2])} else { c(xy[1,1],xy[2,2]) }
		xy <- rbind(xy[1,,drop=FALSE], block.point,		xy[2,,drop=FALSE])
	}
	edge.diag.to.cityblock <- function(diagram,from,to,first) {
		xy <- diagram@edgeList[[1]]@xy
		fromix <- match(from,rownames(xy));
		if (rownames(xy)[fromix+1] != to) { stop("whoops")}
		nxy <- diag.to.cityblock(xy[c(fromix,fromix+1),],first=first)
		xy <- rbind(xy[seq_len(fromix),,drop=FALSE],nxy[2,,drop=FALSE],xy[(fromix+1):nrow(xy),,drop=FALSE])
		diagram@edgeList[[1]]@xy <- xy
		diagram
	}

	if (likeSquares) {
		Set[[3]] <-edge.diag.to.cityblock(Set[[3]],"p2","p11","x")
		Set[[4]]<- edge.diag.to.cityblock(Set[[4]],"p11","p7","x")
		Set[[4]]<- edge.diag.to.cityblock(Set[[4]],"p7","p12","y")
		Set[[3]]<- edge.diag.to.cityblock(Set[[3]],"p12","p6","y")
		Set[[3]]<- edge.diag.to.cityblock(Set[[3]],"p6","p13","x")
		Set[[4]]<- edge.diag.to.cityblock(Set[[4]],"p13","p10","x")
		Set[[4]]<- edge.diag.to.cityblock(Set[[4]],"p10","p14","y")
		Set[[3]]<- edge.diag.to.cityblock(Set[[3]],"p14","p2","y")

	}

	TM <- Set[[1]]
	for (ix in 2:4) {
		TM <- addSetToDrawing(TM,Set[[ix]],set2Name=paste("Set",ix,sep=""))
	}
	VD <- new("VennDrawing",TM,V)
	SetLabels <- .default.SetLabelPositions(VD)
	smidge <- 0.01 * 4
	SetLabels[1,c("x","y")] <- c(0,-2-s-smidge)
	SetLabels[2,c("x","y")] <- c(2+s,2+s+smidge)
	SetLabels[3,c("x","y")] <- c(1,-2  -smidge)
	SetLabels[4,c("x","y")] <- c(2,-1  -smidge)
	SetLabels[1,c("hjust","vjust")] <- c("right","top")
	SetLabels[2,c("hjust","vjust")] <- c("right","bottom")
	SetLabels[3,c("hjust","vjust")] <- c("right","top")
	SetLabels[4,c("hjust","vjust")] <- c("right","top")
	if (!likeSquares) {
		SetLabels[3,c("x","y")] <- SetLabels[3,c("x","y")] + c(-0.5,0.5)
		SetLabels[4,c("x","y")] <- SetLabels[4,c("x","y")] + c(-0.5,0.5)
		SetLabels[3,c("hjust","vjust")] <- c("left","center")
		SetLabels[4,c("hjust","vjust")] <- c("left","center")
	}


	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
	if (s>0) { # want the coordinates of face midpoints if s=0
		Vtemp <- compute.S4(V,doWeights=FALSE,s=0,likeSquares=likeSquares)
		FaceLabels <- VennGetFaceLabels(Vtemp)
	} else { 
		FaceLabels <- .default.FaceLabelPositions(VD)
	}
	VD <- VennSetFaceLabels(VD,FaceLabels)

	VD <- .square.universe(VD,doWeights)


}	
	
