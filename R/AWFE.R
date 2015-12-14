
compute.AWFE <- function(V,doWeights=FALSE,type="battle") {
	if(doWeights) { warning("Weights ignored by AWFE algorithm") }
	n <- NumberOfSets(V)
	# load from cached data, into an environment so we can persuade RCMD check
	# we have it
	env <- new.env()
	loaded <- data(VennDiagrams,package="Vennerable",envir=env)
	stopifnot("VennDiagrams" %in% loaded)
	VennDiagrams <- get("VennDiagrams",envir=env)
	if (! type %in% names(VennDiagrams)) {
		stop(sprintf("No AWFE diagram of type %s known\n",type))
	}
	if (length(VennDiagrams[[type]])< n ) {
		stop(sprintf("No AWFE diagram of type %s and size %d known (max %d)\n",type,n,length(VennDiagrams[[type]])))
	}
	TM <- VennDiagrams[[type]][[n]]
	VD <- new("VennDrawing",TM,V)
	SetLabels <- .default.SetLabelPositions(VD)
	VD <- VennSetSetLabels(VD,SetLabels)
	FaceLabels <- .default.FaceLabelPositions(VD)
	VD <- VennSetFaceLabels(VD,FaceLabels)
	VD <- .square.universe(VD,doWeights=FALSE)
	VD

}

################
# the Edwards construction for n sets;

makeAWFESets <- function(nmax,type="AWFE",hmax) {
	Setlist<- list()
	for (ix in 1:nmax) {
		Setlist[[ix]] <- Setfun(n=ix,type=type,nmax=nmax,hmax=hmax) 
	}
	Setlist
}



###############################
makeAWFE  <- function(n,AWFEnminus1,Setlist) { 	
	if (missing(AWFEnminus1)) {
		thisSet <- Setlist[[1]]
		TM <- thisSet 
		if (n>=2) {
			for (ni in 2:n) {
				TM <- addSetToDrawing(TM,Setlist[[ni]],remove.points=TRUE  )
			}
		}
	} else {
		TM <- addSetToDrawing(AWFEnminus1,Setlist[[n]] ,remove.points=TRUE )
	}
	TM
}


Setfun <- function(n,nmax,hmax,type="AWFE",npoints=200) {
	if (is.na(nmax)) { nmax <- n }
	if(type=="battle") {
		rsize <-  1.5 *(1+   max(0, floor(nmax-3))/(2+2^(nmax-6)) )*1.1
	} else {
		rsize <- 4
	}
	rect1 <- rsize *matrix(c (-1.05,-1,-1.05,1,0,1,0,-1),ncol=2,byrow=TRUE); 
	rect2 <- rsize *matrix(c (-1,0,-1,1.05,1,1.05,1,0),ncol=2,byrow=TRUE); 
	if (n==1) res <- newTissueFromPolygon(points.xy=rect1,Set=1) 
	else if (n==2) res <- newTissueFromPolygon(points.xy=rect2,Set=2)
	else if (n==3 & type %in% c("AWFE","AWFEscale","cog")) {
		 res <- newTissueFromCircle(centre.xy=c(0,0),radius=2,Set=3)
	}
	else {
		if (missing(hmax)) {
			if (type=="AWFE") { hmax <- 0.4 }
		}
		Snf.xy <- set.function(n=n,nmax=nmax,s=(0:npoints)/npoints,hmax=hmax,type=type)
		res <- newTissueFromPolygon(points.xy=Snf.xy,Set=n)
	}
	res
}	


###################################################
# in package use these functions are not normally used because the
# diagrams created by them are cached in a .rda file.....



 
set.function <-  function(n,offset,hmax,nmax,s,type="AWFE") {
	if (type=="AWFE" | type=="AWFEscale") {
		res <- Smithn.function (n-1,offset=offset,hmax=hmax,nmax,s=s,type=type) 
	} else if (type =="cog") {
		res <- cog.function(n,hmax=hmax,nmax=nmax)
	} else if (type=="battle") {
		res <- battle.function(n,nmax=nmax)
	}
	res
}
# if type=cog, s is ignored and the required poolygon is returned

cog.function <- function(n,nmax,hmax=1) {
	if (missing(nmax)) { 
		sf <-  hmax * 1/2^(n-4) # can only scale geometrically 
	} else {
	# want h smoothly spaced from outerlimit hmax down to zero
		# as n increases from 4 to nmax
		deltah <- (hmax) /(nmax-2)
		sf <- (hmax - (n-3)*deltah ) 
# except that n=3 has no zero on line y=0, so need to push it out more
		if (n==3) { sf <- sf +  deltah }

	}
	zeros <- zerotheta(n); zeros <- rev(sort(zeros))
	zeros <- c(zeros,zeros[1])
		phi <- asin(sf) 
		thetaphilist <- list()
		for (ix in seq(1,length(zeros)-1,by=2)) {
			 thetaphilist [[ix]] <- matrix(
			c(zeros[ix],phi ,
			  zeros[ix],0,
			  zeros[ix],-phi ,
			  zeros[ix+1],-phi ,
			  zeros[ix+1],0,
			  zeros[ix+1],phi ),ncol=2,byrow=TRUE)
		}
		thetaphi <- do.call(rbind,thetaphilist )
		xy <- projection.thetaphi(thetaphi,projection="PS")
	xy
}

battle.function <- function(n,nmax) {
	if (missing(nmax)) { 
		stop("nmax must be specified in battle.function") 
	} 
	if (n>nmax) {
		stop("%d is too large for %d\n",n,nmax)
	}
	maxbumps <- 2^(nmax-6)
	r3 <- 1 + 4*maxbumps+1
	r4 <- r3 + (nmax-3)
	dx5 <- 2^(nmax-5)
	if (n==3) { 
		res <- matrix(c( -r3,-r3,-r3,r3,r3,r3,r3,-r3),ncol=2,byrow=TRUE)
	} else if (n==4) {
		res <- matrix(c(-r4,-r4,-r4,r4,-dx5,dx5,dx5,dx5,r4,r4,r4,-r4,dx5,-dx5,-dx5,-dx5),ncol=2,byrow=TRUE)
	} else {
		dxn <- 2^(nmax-n)
		dyn <- nmax-n+1
		if (n==5) {
			xy <- matrix(c(dxn,r3+dyn,dxn,dxn,r3+dyn,dxn),ncol=2,byrow=TRUE)
		} else {
			xy <- matrix(c(dxn,r3+dyn),ncol=2,byrow=TRUE)
			nbumps <- 2^(n-6) 
			xoff <- dxn
			for (ib in seq_len(nbumps)) {
				bxy <- matrix(c(xoff,r3-dyn,xoff+2*dxn,r3-dyn,xoff+2*dxn,r3+dyn,xoff+4*dxn,r3+dyn),ncol=2,byrow=TRUE)
				xy <- rbind(xy,bxy)
				xoff <- xoff + 4 * dxn
			}
			xy <- xy[-nrow(xy),]
			cxy <- matrix(c(r3+dyn,r3+dyn),ncol=2,byrow=TRUE)
			xy <- rbind(xy,cxy)
			mirrorxy <- xy[,c(2,1)]
			mirrorxy <- mirrorxy[nrow(mirrorxy):1,]
			mirrorxy <- mirrorxy[-1,] 
			xy <- rbind(xy,mirrorxy)
		}
		q2 <- xy
		q2[,2] <- -q2[,2]
		q2  <- q2[nrow(q2):1,]; 
		q3 <- xy
		q3[,1] <- -q3[,1]
		q3[,2] <- -q3[,2]
		q4 <- xy
		q4[,1] <- -q4[,1]
		q4 <- q4[nrow(q4):1,]

		res <- do.call(rbind,list(xy,q2,q3,q4))
	}
	# finally scale so set 3 has half width 2
	res <- 2* res/r3
	res
}



# this creates a function of s=(0,1) on a sphere, then projects it as reequested
Smithn.function <- function(n,offset=0,hmax=4,nmax,s,type="AWFE") {
	if (missing(offset)) offset <- 0
# first gets called for 4th set with n=3, will go up to nmax set with n=nmax-1
	if (type=="AWFE") {
		sf <-  hmax * 1/2^(n-3)
	} else if (type=="AWFEscale") {
		if (missing(nmax)) stop("Need to specify nmax for AWFEscale")
		# want h smoothly spaced from outerlimit hmax down to zero
		# as n increases from 3 to nmax-1
		deltah <- (hmax) /(nmax-2)
		sf <- (hmax - (n-2)*deltah ) 
#cat(sf,"\n")
	}
	theta <-  -(s * (2 * pi) ) #clockwise polygons please
	h <-  sf * cos(2^(n-2) * (theta-offset))
	phi <- asin(h)
	thetaphi <- cbind(theta,phi)
	xy <- projection.thetaphi(thetaphi,projection="PS")
	xy
}
 


# the equatorial zeroes of the nth Smith function on (0:1)
zeropos<- function(n) {
	nzeroes <- 2^((n-1))
	zerospacing <- 1/nzeroes
	zerotheta <- -(zerospacing/2)+ zerospacing*( 1:nzeroes)
	zerotheta
	}
zerotheta <- function(n) {  2 * pi* zeropos(n) }
                                    




projection.thetaphi <- function(thetaphi,projection="PS") {
	theta <- thetaphi[,1]; phi <- (thetaphi[,2])
	if (projection=="PS") {
		rho <- cos(phi)/(1-sin(phi))* 2
		x <- rho * cos(theta)
		y <- rho * sin(theta)
		xy <- cbind(x,y)
	} else if (projection=="EC") {
		y <- sin(phi)
		x <- theta
		xy <- cbind(x,y)
	}
	xy
}


thetah.to.xy <- function(thetah,projection) {
	xy <- if (projection=="EC") { 
		thetah
	} 	else { # PS
		phi <- asin(thetah[,2])
		thetaphi <- cbind(thetah[,1],phi)
		projection.thetaphi(thetaphi,projection)
	}
	xy
}



