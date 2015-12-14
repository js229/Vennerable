#warning("Entering Triangles.R")
###########################

.inscribetriangle.feasible <- function(wghts) {
	w0 <- 1- sum(wghts)
	stopifnot(all(wghts <= 1) & all(wghts>=0) & w0>=0)
	wa <- wghts[1];wb <- wghts[2]; wc <- wghts[3]
	Delta <- w0^2 - 4 * wa * wb * wc
	return (Delta>=0)
}
.inscribetriangle.compute <- function (wghts) {
	wa <- wghts[1];wb <- wghts[2]; wc <- wghts[3]
	stopifnot(.inscribetriangle.feasible(wghts)) 
	pa <- (1-wc)
	pb <- (wb+wc-wa-1)
	pc <- wa * (1-wb)
	sc <- if (wa>0) { 
		(-pb-sqrt( pb^2 - 4 * pa * pc))/(2*pa) 
		} else if (wb+wc<1) { 
		(1-wb-wc)/(1-wc) 
		} else { 
		0 
		} 
	sb <- if (sc>0 ) { 1 - wa/sc } else { wc/(1-wb) }
	sa <- wb/(1-sc)
	c(sc,sa,sb) # nb order around triangle
}
.inscribetriangle.inscribe <- function(xy,wghts) {
	scalef <- NA
	
	isfeasible <- .inscribetriangle.feasible(wghts)
	if (!isfeasible) {
		scalef <- 4 * wghts[1]*wghts[2]*wghts[3]/(1-sum(wghts))^2
		scalef <- scalef^(1/3)
		wghts <- wghts / (scalef*1.001)
		isfeasible <- .inscribetriangle.feasible(wghts)
		stopifnot(!isfeasible)
	}
	if (!isfeasible) return(list(feasible=FALSE))
	scab <- .inscribetriangle.compute (wghts)
	inner.xy <- (1-scab)*xy + scab * (xy[c(2,3,1),])
	return(list(feasible=TRUE,inner.xy=inner.xy,scalef=scalef))
}

compute.T3 <- function(V,doWeights=TRUE) {
	stopifnot(NumberOfSets(V)==3)
	wght <- Weights(V)
	VS <- VennSignature(V)
	Vorig <- V
	if (!doWeights) { 
		Wsum <- apply(Indicator(V),1,sum)
		wght[Wsum==0] <- 2
		wght[Wsum==1] <- 4
		wght[Wsum==2] <- 1
		wght[Wsum==3] <- 1
		V@IndicatorWeight[,".Weight"] <- wght
	}


	# with wa, wb, wc on the outside and w0 the universe.

	WeightUniverse <- .WeightUniverse(V)
	WeightVisible <- .WeightVisible(V)
	WeightInvisible <- WeightUniverse-WeightVisible 


	w0ratio <- WeightInvisible/WeightVisible 
	# we ignore w0 from now on and all the other weights sum to one
	wght <- wght/WeightVisible 
	
	# the inner triangle contains wab,wbc,wca and 
	# the innest contains wabc


	# the outer triangle
	wa <- wght[VS=="100"]
	wb <- wght[VS=="010"]
	wc <- wght[VS=="001"]
	outer.weights <- c(wa,wb,wc)
	outer.innerw <- 1 - sum(outer.weights)
	outer.inner.ratios <- outer.weights/outer.innerw #  ratio of each wa, wb,wc to pooled inner weights

	outer.feasible <- .inscribetriangle.feasible(outer.weights)

	# the inner triangle 
	wab <- wght[VS=="110"]
	wbc <- wght[VS=="011"]
	wca <- wght[VS=="101"]
	wabc <-  wght[VS=="111"]

	inner.weights <- c(wab,wbc,wca)
	inner.innerw <- wabc
	# we resclae the inner weights...
	sf <- (sum(inner.weights)+inner.innerw)
	Weight.Inner <- sf * WeightVisible
	if (sf>0) {	
		inner.weights <- inner.weights/sf 
		inner.feasible <- .inscribetriangle.feasible(inner.weights)
	} else {
		inner.feasible <- FALSE
	}
	
	if (inner.feasible & outer.feasible) {
		# whole triangle should have area in Weights
		side <- sqrt(4 * WeightVisible /(3*sqrt(3)))
		angles <- pi/2-c(0,2*pi/3,4*pi/3)
		outer.xy <- t(sapply(angles,function(a)c(x=side * cos(a),y= side * sin(a))))
		
		inner <- .inscribetriangle.inscribe(outer.xy,wghts=outer.weights)
		inner.xy <- inner$inner.xy
		innest <- .inscribetriangle.inscribe(inner.xy,wghts=inner.weights)
		innest.xy=innest$inner.xy
		# finally we construct the outside triangle
		# outer.xy is equilateral with centre at zero, so just scale 
		# if inner triangle has area A and rim has area A' then scaling 
		# is (A'+A)=s^2 A so s^2=1+A'/A. A'/A is the w0ratio calculated above
		outest.xy <- outer.xy * sqrt( 1+ w0ratio)
	} else { 
		if (inner.feasible) { # but not outer
			# so we make the inner equilateral of area Weight.Inner
			side <- sqrt(4 * Weight.Inner /(3*sqrt(3)))
			angles <- pi/6-c(0,2*pi/3,4*pi/3)
			inner.xy <- t(sapply(angles,function(a)c(x=side * cos(a),y= side * sin(a))))
			innest <- .inscribetriangle.inscribe(inner.xy,wghts=inner.weights)
			innest.xy=innest$inner.xy
			# 
			outer.heights <- 2 * outer.inner.ratios * Weight.Inner /side # TODO wrong?
			outer.distance <- side/(2*sqrt(3))+ outer.heights
			outer.angles <- pi/2-c(0,2*pi/3,4*pi/3)
			outer.x <- outer.distance * cos(outer.angles)
			outer.y <- outer.distance *sin(outer.angles)
			outer.xy <- matrix(c(outer.x,outer.y),ncol=2,byrow=FALSE)
			# TODO THIS IS QUITE WRONG
			outest.xy <- outer.xy * sqrt( 1+ w0ratio)
		} else { # inner and out infeasible
		stop("Can'y yet cope with inner and outer infeasible triangles")
		}
	}
	rownames(outer.xy) <- paste("to",1:3,sep="")
	rownames(inner.xy) <- paste("ti",1:3,sep="")
	rownames(innest.xy) <- paste("tt",1:3,sep="")

	outline.a.xy  <- do.call(rbind,list(outer.xy[1,,drop=FALSE],inner.xy[1,,drop=FALSE],innest.xy[1,,drop=FALSE],innest.xy[2,,drop=FALSE],inner.xy[3,,drop=FALSE]))
	outline.b.xy  <- do.call(rbind,list(outer.xy[2,,drop=FALSE],inner.xy[2,,drop=FALSE],innest.xy[2,,drop=FALSE],innest.xy[3,,drop=FALSE],inner.xy[1,,drop=FALSE]))
	outline.c.xy  <- do.call(rbind,list(outer.xy[3,,drop=FALSE],inner.xy[3,,drop=FALSE],innest.xy[3,,drop=FALSE],innest.xy[1,,drop=FALSE],inner.xy[2,,drop=FALSE]))

	VDP1 <- newTissueFromPolygon(points.xy=outline.a.xy,Set=1)
	VDP2 <- newTissueFromPolygon(points.xy=outline.b.xy,Set=2)
	VDP3 <- newTissueFromPolygon(points.xy=outline.c.xy,Set=3)
#browser()
	TM <- addSetToDrawing (drawing1=VDP1 ,drawing2=VDP2, set2Name="Set2")
	TM <- addSetToDrawing (drawing1=TM ,drawing2=VDP3, set2Name="Set3")
	VD <- new("VennDrawing",TM,Vorig)
	SetLabels <- .default.SetLabelPositions(VD)
	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
	VD <- VennSetFaceLabels(VD,FaceLabels)
	VD <- .square.universe(VD,doWeights)



}



	

