#warning("Entering Circles.R")

################################################
# some geometries 
# r is radius of circle
# d distance from origin


TwoCircles <- function(r,d,V) {
	if (length(r) !=2 ) {
		if (length(r)==1 ) {
			r <- rep(r,2)
		} else {
			stop("Need two circle radii")
		}
	}
	# Two circles, radius r1 and r2, distance d apart
	centres <- matrix(c(-d/2,0,d/2,0),ncol=2,byrow=TRUE)
	VDC1 <- newTissueFromCircle(centres[1,],radius=r[1],Set=1); 
	VDC2 <- newTissueFromCircle(centres[2,],radius=r[2],Set=2); 
	TM <- addSetToDrawing (drawing1=VDC1,drawing2=VDC2,set2Name="Set2",remove.points=TRUE)
	V2 <- new("VennDrawing",TM,V)
	SetLabels <- .circle.SetLabelPositions(V2,radii=r,centres=centres)
	V2 <- VennSetSetLabels(V2,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(V2)
	V2 <- VennSetFaceLabels(V2,FaceLabels)

	V2
}

.circle.SetLabelPositions <- function(object,radii,centres){
	yscale <- diff(VisibleRange(object)[,2]);
	smidge <- 0.01*yscale
	xy <- centres
	abovebelow <- rep(1,NumberOfSets(object))
	if (NumberOfSets(object)==3) {
		abovebelow <- c(1,-1,-1)
	}
	xy[,2] <- xy[,2]+abovebelow*(radii+smidge)
	VLabels <- data.frame(Label=rep("unset",nrow(xy)),x=NA,y=NA,hjust=I("center"),vjust=I("bottom"))
	VLabels$vjust <- ifelse(abovebelow>0,"bottom","top")
	VLabels[,2:3] <- xy
	VLabels$Label <- VennSetNames(as(object,"Venn"))
	VLabels
}


# 2-d diagrams
compute.C2 <- function(V,doWeights=TRUE,doEuler=FALSE) {
	Vcalc <- V
	if(!doWeights) {
		if (!doEuler) {
		# want each area, even if empty
			wght <- Weights(Vcalc )
			wght <- 1+ 0*wght
			Weights(Vcalc ) <- wght
		} else {
			wght <- Weights(Vcalc )
			wght [wght !=0] <- c(3,2,2,1)[wght !=0]
			Weights(Vcalc ) <- wght
		}
	}
	dList <- .Venn.2.weighted.distance (Vcalc,doEuler) # returns radii of two circles and their distance
	r1 <- dList$r1;r2 <- dList$r2; d <- dList$d; 
	C2 <- TwoCircles(r=c(r1,r2),d=d,V) # d in TwoCircles is distance of centre from origin
	C2 <- .square.universe(C2,doWeights)
	C2
	
}


#####################################
# given four weights we seek to find two circles all of whose
# interesction areas are proportional to those weights
.Venn.2.weighted.distance <- function(V,doEuler) {
	if (NumberOfSets(V) != 2) {
		stop(sprintf("Wrong number of sets (%d) in .Venn.2.weighted.distance",NumberOfSets(V)))
	}
	# cf Chow and Ruskey, 200?
	Weight <- Weights(V)
	wab <- Weight["00"] #not used 
	wAb <- as.numeric(Weight["10"])
	waB <- as.numeric(Weight["01"])
	wAB <- as.numeric(Weight["11"])

	inneroff <- 1 # set = 2 to put inner circles completely inside
	outeroff <- 1.01 # =1 to have exact adjacency
	r1 <- sqrt( (wAb+wAB)/pi)
	r2 <- sqrt( (waB+wAB)/pi) # area proportional to weights

	if (wAb==0) { #inside the other one 
		if (doEuler) {
			d <- (r2-r1)/inneroff
		} else { # bodge it outside to ensure the Venntersection..doesnt really make much sense
			d <- (r2-r1)+0.1*r1
		}
	} else if (waB==0) { # also
		if (doEuler) {
			d <- (r1-r2)/inneroff
		} else { # bodge it outside to ensure the Venntersection..doesnt really make much sense
			d <- (r1-r2)+0.1*r1
		}
	} else if (wAB==0) {
		d <- outeroff *  (r1+r2)  # no intersection
		if (!doEuler) {
			d <- 0.95 * d
		}
		
	} else  { # all nonzero
	alpha  <- function(d,r1,r2) {
		alphaD <- (d^2+r1^2-r2^2)
		alphaD <- ifelse(alphaD ==0,alphaD ,alphaD / (2 * r1 * d) )
		alphaD <- ifelse(alphaD>1,1,alphaD)
		alphaD <- ifelse(alphaD< -1,-1,alphaD)
		alpha <-  2 *acos( alphaD )
		alpha
	}
	sigma <- function(d,r1,r2) { 
		alpha <- alpha(d,r1,r2)
		betaD <- (d^2+r2^2-r1^2); 
		betaD <- ifelse(betaD==0,betaD,betaD/ (2 * r2 * d) )
		betaD <- ifelse(betaD< (-1),-1,betaD )
		betaD <- ifelse(betaD> 1,1,betaD )
		beta <-   2 *acos( betaD )
		sigma <- 0.5 * r1^2*(alpha-sin(alpha))+ 0.5 * r2^2 * (beta-sin(beta))
		sigma
	}
	sigmaminuswAB <- function(d,r1,r2,wAB) {
		sigma(d,r1,r2) - wAB
	}
	d <-	uniroot(sigmaminuswAB ,lower=r1-r2,upper=r1+r2,r1=r1,r2=r2,wAB=wAB)$root
	} # wAB != 0 
	list(r1=r1,r2=r2,d=d )

}

.twoCircleIntersectionPointsObsolete <- function(r,xy,i1=1,i2=2) {
	r <- r[c(i1,i2)] 
	xy <- xy[c(i1,i2),]
	# two circles with radius r[1], r[2], centres xy[1,], xy[2,]
	d <- sqrt(sum((xy[1,]-xy[2,])^2))
	# do they intersect? 
	if (r[1]+r[2] > d  & d > abs(r[1]-r[2])) {
		d1 <- (d^2-r[2]^2+r[1]^2)/(2*d)
		d2 <- abs(d - d1)
		xyvec <- xy[2,] - xy[1,]
		midchord <- xy[1,] + (d1/d)*xyvec
		normal <- c(-xyvec[2],xyvec[1])
		normal <- normal/sqrt(sum(normal^2))
		halfchord <-  sqrt(r[1]^2-d1^2)
		p1 <- midchord + (halfchord)*normal
		p2 <- midchord - (halfchord)*normal
		p12 <- rbind(p1,p2)
	} else {
		p12 <- matrix(NA,nrow=2,ncol=2)
	}
	rownames(p12) <- paste("p",i1,i2,c("a","b"),sep="")
	p12
}



ThreeCircles <- function(r,x,y,d,angles,V) {
	if (missing(x) | missing(y)) {
		if (missing(d)) {
			stop("Need x and y or d")
		}
		if (missing(angles)) {
			angles <- pi/2-c( 0, 2*pi/3, 4 * pi/3)
		}
		x <- d*cos(angles)
		y <- d*sin(angles)
	}
	if (length(r) !=3 ) {
			if (length(r)==1 ) {
				r <- rep(r,3)
			} else {
				stop("Need three circle radii")
			}
	}

	centres <- matrix(c(x,y),ncol=2,byrow=FALSE)
	
	nodes <- 1; TM <- NA
	while (!inherits(TM,"TissueDrawing") & nodes < 10) {
		VDC1 <- newTissueFromCircle(centres[1,],radius=r[1],Set=1,nodes=nodes); 
		VDC2 <- newTissueFromCircle(centres[2,],radius=r[2],Set=2,nodes=nodes); 
		TM <- addSetToDrawing (drawing1=VDC1,drawing2=VDC2,set2Name="Set2")
		nodes <- nodes+1
	} 
	if (nodes>=10) stop("Can't join circles")
	nodes <- 1; TM2 <- NA
	while (!inherits(TM2,"TissueDrawing")  & nodes < 10) {
		VDC3 <- newTissueFromCircle(centres[3,],radius=r[3],Set=3,nodes=nodes); 
		TM2 <- addSetToDrawing (drawing1=TM,drawing2=VDC3,set2Name="Set3")
		nodes <- nodes+1
	} 
	if (nodes>=10) stop("Still can't join circles")

	C3 <- new("VennDrawing",TM2,V)
	SetLabels <- .circle.SetLabelPositions(C3,radii=r,centres=centres)
	C3 <- VennSetSetLabels(C3,SetLabels)
}


.pairwise.overlaps <- function(V) {
	VI <- Indicator(V)
	VIsum <- apply(VI,1,sum)
	Vpairs <- unique(VI[VIsum==2,])
	Vwhich <- which(Vpairs==1,arr.ind=TRUE)
	Vindex <- Vpairs
	Vindex [Vwhich] <- Vwhich[,2]
	isdisjoint <- function(vrow) {
		vsub <- V[,vrow]
		vsum <- apply(Indicator(vsub),1,sum)
		overweight <- Weights(vsub)[vsum==2]
		overweight==0
	}
	Vpairs <- data.frame(Vpairs,check.names=FALSE)
	Vpairs$Disjoint <- apply(Vindex,1,isdisjoint)
	Vpairs
}


compute.C3 <- function(V,doWeights=TRUE) {
	doEuler <- TRUE
if (doWeights) {	
		overlaps <- .pairwise.overlaps(V)
		dList12 <- .Venn.2.weighted.distance (V[,c(1,2)],doEuler ) # returns radii of two circles and their distance
		dList23 <- .Venn.2.weighted.distance (V[,c(2,3)],doEuler ) # 
		dList31 <- .Venn.2.weighted.distance (V[,c(3,1)],doEuler ) #
		
		disjointcount <- sum(overlaps$Disjoint)
		dp <- c( cp= dList12$d, b = dList23$d, a = dList31$d)
		smidge <- 1 - 1e-4
		if (disjointcount==3) {
			# all disjoint, arrange with centres in equilateral triangle
			dmax <- max(dp)
			dp <- dp*0+dmax
		} else if (disjointcount==2) {
			#one is disjoint from both the others
			# set it at the same distance from both
			# and far enough that it will be a triangle 
			conjoint <- overlaps[!overlaps$Disjoint,-ncol(overlaps)]
			conjointix <-  ( which(conjoint==0)%%3 +1) 
			conjointdistance <-  dp[conjointix ]
			disjointdistances <- dp[-conjointix]
			disjointdistances <- max(disjointdistances,conjointdistance/2)
			dp[-conjointix] <- disjointdistances 
		} else if (disjointcount==1) {
			# only one missing intersection; we set its distance to
			# force a (near) straightline
			# nb this will still fail if the intersections are large enough
			disjoint <- overlaps[overlaps$Disjoint,-ncol(overlaps)]
			disjointix <- ( which(disjoint ==0)%%3 +1) 
			conjointdistances <-  dp[-disjointix ]
		 	dp[disjointix] <- sum(conjointdistances)*smidge
		} 
		cp= dp["cp"]; b = dp["b"]; a=dp["a"]

		# can we satisfy the triangle inequality? if not, bodge it
		if (a+b < cp) {
			cp <-  (a+b)*smidge
		}
		if (b+cp < a) {
			a <-  (b+cp)*smidge
		}
		if (cp+a < b ) {
			b <- (cp+a)*smidge
		}
		# pick a centre for the first one	
		c1 <- c(0, dList12$r1)
		# then place the others
		if (a==0 || b==0 || cp==0) {
			o21 <- cp *c(1,0) ; o31 <- a * c(0,1)
		} else {
			# use the SSS rule to compute angles in the triangle
			CP <- acos( (a^2+b^2-cp^2)/(2*a*b))
			B <- acos( (a^2+cp^2-b^2)/(2*a*cp))
			A <- acos( (cp^2+b^2-a^2)/(2*b*cp)) 
			# arbitrarily bisect one and calculate xy offsets from first centre
			theta <- B/2
			o21 <- cp * c( sin(theta), - cos(theta))
			o31 <- a * c( -sin(B-theta),	-cos(B-theta))
		}
		c2 <- c1 + o21
		c3 <- c1 + o31
		C3 <- ThreeCircles(r=c(dList12$r1,dList12$r2,dList31$r1),
			x=c(c1[1],c2[1],c3[1]),y=c(c1[2],c2[2],c3[2]),V=V)
	} else {
		C3 <- ThreeCircles(r=0.6,d=0.4,V=V)
	}
	C3 <- .square.universe(C3,doWeights=doWeights)
	FaceLabels <- .default.FaceLabelPositions(C3)
	C3 <- VennSetFaceLabels(C3,FaceLabels)

	C3
}

