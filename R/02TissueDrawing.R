#warning("Entering 02TissueDrawing")
# # $Id: 02Drawing.R,v 1.1 2007/10/16 10:59:52 js229 Exp $

# $Log: 02Drawing.R,v $

setGeneric("PlotSetBoundaries",function(drawing,gp){standardGeneric("PlotSetBoundaries")})
setGeneric("PlotNodes",function(drawing,gp){standardGeneric("PlotNodes")})
setGeneric("PlotFaces",function(drawing,faceNames,gp,arrow,colourAlgorithm){standardGeneric("PlotFaces")})



#####
# todo: fix adding edges to faces with repeated points
# wholly internal new sets (ie non simple faces)

##########################
# we start with points
# they dont currently have a class of their own, but they have unique
# names and unique xy coordinates and are stored as a 1x2 matrix with the row name as the name of the point
# and will be kept in a named nodeList by any drawing


###########################
# points are joined by edges of course,
# they always go to and from named points stored in the class,
# and so the set of them (will be in a named edgeList) define the topological graph structure
# but each edge also defines the xy-coordinates of the line in the drawing
# not true that the only named points on an edge are at the beginning and end
# but any internal named points are guaranteed not to be at intersections of edges

setClass("VDedgeDrawn",representation(from="character",to="character",visible="logical",bb="matrix"),
	prototype(from=character(0),to=character(0),visible=TRUE))

## these generics have methods for all edge types inheriting from VDedgeDrawn
# called only within this file and documentation vignette
setGeneric(".identical",function(edge1,edge2){standardGeneric(".identical")})
setGeneric(".reverseEdge",function(edge){standardGeneric(".reverseEdge")})
setGeneric(".midpoint",function(edge)standardGeneric(".midpoint"))
setGeneric(".checkPointOnEdge",function(edge,point.xy)standardGeneric(".checkPointOnEdge"))
setGeneric(".splitEdgeAtPoint",function(edge,point.xy)standardGeneric(".splitEdgeAtPoint"))
setGeneric(".findIntersectionByType",function(edge1,edge2){standardGeneric(".findIntersectionByType")})
setGeneric(".edge.to.xy",function(edge,dx){standardGeneric(".edge.to.xy")})


setClass("VDedgeSector",representation("VDedgeDrawn",
	radius="numeric",fromTheta="numeric",toTheta="numeric",centre="numeric",hand="numeric"))
setClass("VDedgeLines",representation("VDedgeDrawn",xy="matrix"))
# not implemented:
#setClass("VDedgeFunction",representation(Set="numeric",s="numeric",f="function"))


setMethod("show","VDedgeDrawn",function(object) {
	res <- c(from=object@from,to=object@to,visible=object@visible)
	show(res)
	invisible(	object)
})


setMethod("show","VDedgeSector",function(object) {
	show(as(object,"VDedgeDrawn"))
	res <- c(centre=paste(object@centre,collapse=","),
		radius=object@radius,fromTheta=object@fromTheta,toTheta=object@toTheta,hand=object@hand)
	show(res)
	invisible(	object)
})


setMethod(".identical",c("VDedgeSector","VDedgeSector"), function(edge1,edge2) {
	fequal(edge1@radius,edge2@radius) & 
	fequal(edge1@fromTheta,edge2@fromTheta) &
	fequal(edge1@toTheta,edge2@toTheta) &
	all(fequal(edge1@centre,edge2@centre)) &
	fequal(edge1@hand,edge2@hand) 
})

setMethod(".identical",c("VDedgeLines","VDedgeLines"), function(edge1,edge2) {
	if (nrow(edge1@xy) != nrow(edge2@xy)) return(FALSE)
	all(fequal(edge1@xy,edge2@xy))
})
setMethod(".identical",c("VDedgeSector","VDedgeLines"), function(edge1,edge2) {
	return(FALSE) 
})
setMethod(".identical",c("VDedgeLines","VDedgeSector"), function(edge1,edge2) {
	return(FALSE) 
})

.findIntersection<- function(edge1,edge2) {
	bb1 <- edge1@bb; bb2 <- edge2@bb;	smudge <-  (.Machine$double.eps ^ 0.5)

	x1 <- bb1[,1]; x2 <- bb2[,1]
	xok <- x2[2] >= x1[1] - smudge & x2[1] <= x2[2] +smudge 
	y1 <- bb1[,2]; y2 <- bb2[,2]
	yok <- y2 [2] >= y1 [1] - smudge & y1 [1] <= y2 [2] +smudge 
	if (xok & yok) {
		found <- .findIntersectionByType(edge1,edge2)
	} else {
		found <- matrix(ncol=2,nrow=0)
	}
	found
}



newEdgeSector <- function(centre,hand=1,from,to,fromTheta,toTheta,radius,visible=TRUE) {
	if (missing (from)) { from <- "p1" }
	if (missing (to)) if ( fequal( (fromTheta-toTheta) %% (2 * pi),0)) { to <- from } else { to <- "p2"}
	bb <- rbind( centre-radius,centre+radius)
	lh <- new("VDedgeSector",centre=centre,hand=hand,from=from,to=to,fromTheta=fromTheta,toTheta=toTheta,radius=radius,visible=visible,bb=bb)
	lh
}

newEdgeLines <- function(from,to,xy,visible=TRUE) {
	bb <- rbind(apply(xy,2,min),apply(xy,2,max))
	lh  <- new("VDedgeLines",from=from,to=to,xy=xy,visible=visible,bb=bb)
	lh
}



setMethod("show","VDedgeLines",function(object) {
	show(as(object,"VDedgeDrawn"))
	res <- c(npoints=nrow(object@xy))
	show(res)
	invisible(	object)
})

	
.find.sector.sector.intersection <- function(edge1,edge2) {
	r1 <- edge1@radius
	r2 <- edge2@radius
	centre.diff <- edge2@centre - edge1@centre
	d <- sqrt(sum(centre.diff^2))
	if (d > (r1+ r2)) { # disjoint
		intersectionPoints  <- (matrix(numeric(0),ncol=2))
	} else if (d < abs(r1-r2)) { # one inside the other
		intersectionPoints  <- (matrix(numeric(0),ncol=2))
	} else if (d == abs(r1-r2)) { # one inside, tangent touching
		p1 <- matrix((edge1@centre+ r1*(r1-r2)*centre.diff/d^2),ncol=2) # from simplification of below
		rownames(p1)<-"i1"
		intersectionPoints  <- (p1)
	} else if (d == (r1+r2)) { # adjacent, tangent touching
		p1 <- matrix((edge1@centre+ (r1/(r1+r2))*centre.diff),ncol=2) # 
		rownames(p1)<-"i1"
		intersectionPoints  <- (p1)
	} else { # two circles intersect in two points 
		if (d==0 & r1==r2){stop("Identical circles")}
		d1 <- (d^2 - r2^2+ r1^2) /( 2* d)
	 	d2 <- d - d1
		# half-height from the chord to the intersection
		yh <- (1/(2*d))* sqrt(4*d^2*r1^2-(d^2-r2^2+r1^2)^2)

		dvec <- centre.diff/d # normalised vector from 1 to 2
		nvec <- c(-dvec[2],dvec[1]) # normal to that
		p1 <- matrix(edge1@centre + (d1)*dvec + yh * nvec,ncol=2);rownames(p1)<-"i1"
		p2 <- matrix(edge1@centre + (d1)*dvec - yh * nvec,ncol=2);rownames(p2)<-"i2"
		if (yh==0) {
			intersectionPoints  <- p1
		} else {
			intersectionPoints  <- (rbind(p1,p2))
		}
	}
	if (nrow(intersectionPoints)>0) {
		onedge1 <- sapply(seq_len(nrow(intersectionPoints  )),function(ix){
		.checkPointOnEdge(edge1,intersectionPoints  [ix,,drop=FALSE])})
		onedge2 <- sapply(seq_len(nrow(intersectionPoints  )),function(ix){
		.checkPointOnEdge(edge2,intersectionPoints  [ix,,drop=FALSE])})

		intersectionPoints <- intersectionPoints [onedge1 & onedge2 ,,drop=FALSE]
	}
	dimnames(intersectionPoints) <- NULL

	intersectionPoints
}
setMethod(".findIntersectionByType",c("VDedgeSector","VDedgeSector"), .find.sector.sector.intersection)

.det2 <- function(m) { m[1,1] * m[2,2] - m[1,2]* m[2,1] }

.find.linear.intersection <- function(xy1,xy2) {

	# See Mathworld line-line intersection

	denommat <- .det2(rbind( xy1[1,]-xy1[2,],  xy2[1,]-xy2[2,]))
	if (fequal(0,denommat)) { # parallel 
		# it may still be colinear, but in that case we choose to describe this as no intersection,
		# as the endpoints of the edges will also be in other edges which won't be colinear with one of them 
		#p1 <- xy1[1,];
		#p2 <- xy1[2,];
		#p3 <- xy2[1,]
		#detm <- p1[1]*(p2[2]-p3 [2])+p2[1]*(p3 [2]-p1[2])+p3 [1]*(p1[2]-p2[2])
		# no, so no intersection
		#if (!fequal(0,denommat)) {
			return(c(NA,NA))
		#}
	}
	xm <- .det2( matrix( c(.det2(xy1) , xy1[1,1]-xy1[2,1],.det2(xy2), xy2[1,1]-xy2[2,1]),nrow=2,byrow=TRUE))/denommat
	ym <- .det2( matrix( c(.det2(xy1) , xy1[1,2]-xy1[2,2],.det2(xy2), xy2[1,2]-xy2[2,2]),nrow=2,byrow=TRUE))/denommat
	
	# parameter s1 goes from 0 to 1 from first to second point

	if (!fequal(xy1[2,1]-xy1[1,1],0)) {
		s1 <- (xm - xy1[1,1])/(xy1[2,1]-xy1[1,1]); 
	} else {
		s1 <- (ym - xy1[1,2])/(xy1[2,2]-xy1[1,2]);
	}
	if (!fequal(xy2[2,1]-xy2[1,1],0)) {
		s2 <- (xm - xy2[1,1])/(xy2[2,1]-xy2[1,1])
	} else {
		s2 <- (ym - xy2[1,2])/(xy2[2,2]-xy2[1,2])
	}
	smudge <-  (.Machine$double.eps ^ 0.5)

	if ( s1< 0-smudge | s1> 1 +smudge | s2 < 0-smudge | s2 > 1+smudge) {
		return(c(NA,NA))
	} else {
		return(c(xm,ym))
	}


}

.find.linear.circle.intersection <- function(xy,sector) {
	intersectionPoints <- matrix(nrow=0,ncol=2)
	xr <- sector@centre[1]; yr <- sector@centre[2]
 	r <- sector@radius

	# line is y=mx+c
	m <- ( xy[2,2] - xy[1,2] )/(xy[2,1] - xy[1,1]) 
	if (is.infinite(m)) { # must be easier way than special casing eg x<->y
		v <- as.numeric(xy[1,1]) # vertical line at x=v
		if ( v <  xr - r | v > xr + r ) {
			return(intersectionPoints)
		} 
		if ( v== xr - r) {
			intersectionPoints <- rbind(intersectionPoints,c(xr-r,yr))
			return(intersectionPoints)
		}
		if (v == xr+r ) {
			intersectionPoints <- rbind(intersectionPoints,c(xr+r,yr))
			return(intersectionPoints)
		}
		hoff <- sqrt(r*r - (xr-v)^2)
		intersectionPoints <- rbind(intersectionPoints,c(v,yr+hoff))
		intersectionPoints <- rbind(intersectionPoints,c(v,yr-hoff))
		return(intersectionPoints)
	}


		
	cc <-  xy[1,2] - m * xy[1,1] 
	# circle is (x-xr)^2+(y-yr)^2 = r^2
	# (x-xr)^2 + (mx+c-yr)^2 = r^2
	# x^2(1+m^2) + 2 x (m*c-m*yr-xr) + xr^2 + (c-yr)^2 - r^2 = 0
	
	qa <- (1+m^2)
	qb <- 2 * (m *cc - m * yr - xr)
	qc <- xr^2 + (cc-yr)^2 - r^2


	det <- (qb*qb - 4 * qa* qc)
	if (det <0) { }
	if (det ==0) {
		x0 <- (-qb/(2*qa))
		y0 <- m * x0 + cc
		intersectionPoints <- rbind(intersectionPoints,c(x0,y0))
	}
	if (det >0) {
		x1 <- (-qb + sqrt(det))/(2*qa)
		x2 <- (-qb - sqrt(det))/(2*qa)
		y1 <- m * x1 + cc
		y2 <- m * x2 + cc
		intersectionPoints <- rbind(intersectionPoints,c(x1,y1))
		intersectionPoints <- rbind(intersectionPoints,c(x2,y2))
	}
	intersectionPoints 
}


setMethod(".findIntersectionByType",c("VDedgeLines","VDedgeSector"), function(edge1,edge2) {
	# first find all the intersection points assuming infinite line and complete circle
	intersectionPoints <- matrix(nrow=0,ncol=2)
	for(i1 in 1:(nrow(edge1@xy)-1)) {
		ict <-.find.linear.circle.intersection(xy=edge1@xy[i1:(i1+1),],sector=edge2)
		if (nrow(ict)>0) {
			intersectionPoints <- rbind(intersectionPoints,ict)
		}
	}
	intersectionPoints <- .removeDuplicates(intersectionPoints )
	if (nrow(intersectionPoints)>0) {
		onedge1 <- sapply(seq_len(nrow(intersectionPoints)),function(ix){
			.checkPointOnEdge(edge1,intersectionPoints[ix,,drop=FALSE])})
		onedge2 <- sapply(seq_len(nrow(intersectionPoints)),function(ix){
			.checkPointOnEdge(edge2,intersectionPoints[ix,,drop=FALSE])})
		intersectionPoints <- intersectionPoints [onedge1 & onedge2 ,,drop=FALSE]
	}
	dimnames(intersectionPoints) <- NULL
	intersectionPoints
}
)
setMethod(".findIntersectionByType",c("VDedgeSector","VDedgeLines"), function(edge1,edge2).findIntersectionByType(edge2,edge1))

.find.linear.linear.intersection <-  function(edge1,edge2) {
	intersectionPoints <- matrix(nrow=0,ncol=2)
	for(i1 in 1:(nrow(edge1@xy)-1)) {
		for (i2 in 1:(nrow(edge2@xy)-1)) {
			ict <-.find.linear.intersection(xy1=edge1@xy[i1:(i1+1),],xy2=edge2@xy[i2:(i2+1),])
			if (!any(is.na(ict))) {
				intersectionPoints <- rbind(intersectionPoints,ict)
			}
		}
	}
	intersectionPoints <- .removeDuplicates(intersectionPoints )
	intersectionPoints 
}
setMethod(".findIntersectionByType",c("VDedgeLines","VDedgeLines"),.find.linear.linear.intersection)

.removeDuplicates <- function(xypoints) {
	if (nrow(xypoints)<2) {return(xypoints)}
	for (ix in 1:(nrow(xypoints)-1)) {
		for (jx in (ix+1):nrow(xypoints)) {
#cat(ix,jx,"\n")
#show(xypoints)
			if (!any(is.na(xypoints[c(ix,jx),1]))) {
				dist2 <- sum((xypoints[ix,]-xypoints[jx,])^2)
				if (dist2 < .Machine$double.eps ) {
					xypoints[jx,] <- NA
				}
			}
		}
	}
	xypoints <- xypoints[!is.na(xypoints[,1]),,drop=FALSE]
}


setGeneric("joinEdges",function(object1,object2)standardGeneric("joinEdges"))
setMethod("joinEdges",c("VDedgeLines","VDedgeLines"),function(object1,object2).join.lines(object1,object2))
.join.lines <- function(object1,object2) {
	# assumes they do join!
	visibility <- c(object1@visible,object2@visible);stopifnot(visibility[1]==visibility[2])
	xy <- rbind(object1@xy,object2@xy[-1,,drop=FALSE])
	newEdge <- newEdgeLines(from=object1@from,to=object2@to,visible=visibility[1],xy=xy)
	newEdge
}
setMethod("joinEdges",c("VDedgeSector","VDedgeSector"),function(object1,object2).join.arcs(object1,object2))
.join.arcs <- function(object1,object2) {
	# assumes they do join and have same radius centre and hand (and as a side effect will have the same bb)
	visibility <- c(object1@visible,object2@visible);stopifnot(visibility[1]==visibility[2])
	if (!fequal((object1@toTheta - object2@fromTheta) %% (2 * pi),0)) {
		stop(sprintf("Sectors joined at different thetas %g %g\n",object1@toTheta, object2@fromTheta))
	}
	if (object1@toTheta <= 0 & object2@fromTheta > 0 ) {
		object2@fromTheta <- object2@fromTheta - 2*pi # not used but you get the idea
		object2@toTheta <- object2@toTheta - 2 *pi
	}
	newEdge <- object1; 
	newEdge@toTheta <- object2@toTheta
	newEdge@to <- object2@to
	newEdge
}
##################
# for when we actually want to plot an edge:


###
# sector related xy calculations
# normal angles and atan2 goes anticlockwise from 3pm
# but we want +ve handed arcs to go clockwise

.point.xy.to.theta <- function(point,centre) {
	pointCentre <- as.numeric(point)-centre
	theta <- atan2(pointCentre[2],pointCentre[1])
#	theta <- -theta
	theta <- theta %% ( 2* pi)
}
.theta.to.point.xy <- function(theta,r,centre) {
#	theta <- -theta
	x <- r* cos(theta)+centre[1];y <- r*sin(theta)+centre[2]
	cbind(x,y)
}

sector.to.xy <-  function(edge,dx=.05) {
	r <- edge@radius;
	hand <- edge@hand
		thetafrom <- edge@fromTheta;thetato <- edge@toTheta
		# if hand > 0 we always go anti-clockwise, decreasing thetafrom 
		# so thetato must be less  than thetafrom
		if (hand>0) {
			if (thetato>thetafrom)  { thetato <- thetato - 2* pi } 
			thetadist <- thetafrom - thetato
		}
		else {
			if (thetato<thetafrom)  { thetato <- thetato + 2* pi } 
			thetadist <- thetato - thetafrom
		}
	arclength <- (thetadist) * r
	nintervals <- max(3,arclength/dx)
	theta <- seq(from=thetafrom,to=thetato,length=nintervals)
	xy <- .theta.to.point.xy(theta,r,edge@centre)
}

.normalise.sector <- function(edge) {
	# we guarantee 2pi>from>0 and from>to>-2pi
	stopifnot(edge@fromTheta>=0 & edge@fromTheta<=2*pi)
	if (edge@toTheta>edge@fromTheta) {
		edge@toTheta <- edge@toTheta- (2*pi)
	}
	stopifnot(edge@fromTheta>=edge@toTheta & edge@toTheta>=-2*pi)
	edge
}

setMethod(".edge.to.xy",c("VDedgeSector","numeric"),
	function(edge,dx){sector.to.xy(edge)})
setMethod(".edge.to.xy",c("VDedgeSector","missing"),
	function(edge,dx){sector.to.xy(edge)})
setMethod(".edge.to.xy",c("VDedgeLines","numeric"),
	function(edge,dx) {	edge@xy})
setMethod(".edge.to.xy",c("VDedgeLines","missing"),
	function(edge,dx) {edge@xy})

############################
# reversing edges:

setMethod(".reverseEdge","VDedgeLines",function(edge){
	edge@xy <- edge@xy[ rev(seq_len(nrow(edge@xy))),]
	temp <- edge@from
	edge@from <- edge@to
	edge@to <- temp
	edge
})
setMethod(".reverseEdge","VDedgeSector",function(edge){
	edge@hand <- - edge@hand
	temp <- edge@from
	edge@from <- edge@to
	edge@to <- temp
	temp <- edge@fromTheta
	edge@fromTheta<- edge@toTheta
	edge@toTheta <- temp
	edge
})




#######################################
# check if a new point is on an existing edge

setMethod(".checkPointOnEdge",c("VDedgeSector"),function(edge,point.xy) {
	r <- edge@radius
	rp <- sqrt( sum( (point.xy - edge@centre)^2))
	theta <- .point.xy.to.theta(point.xy,edge@centre)
	if (edge@from==edge@to) { 
		intheta <- TRUE
	} else {
		if (edge@toTheta>=0) {
			intheta <- (edge@fromTheta >= theta & theta >= edge@toTheta)
		} else {
#			intheta <- (0<=theta & theta <= edge@fromTheta ) | ( 2*pi+edge@toTheta <= theta & theta <= 2*pi)
			intheta <- ( edge@fromTheta >= theta) | ( 2*pi+edge@toTheta <= theta )
		}
		if (edge@hand<0) intheta <- !intheta
	}
	intheta & isTRUE( all.equal(r,rp))

	}
)

fequal <- function(x,y) {
	abs(x-y) < (.Machine$double.eps ^ 0.5)
}

.find.point.on.EdgeLines <- function(edge,point.xy) {
	ison <- FALSE
	for (nix in 1:(nrow(edge@xy)-1)) {
		line <- edge@xy[c(nix,nix+1),]
		p1 <- line[1,];
		p2 <- line[2,];
		# first check if points are colinear
		detm <- p1[1]*(p2[2]-point.xy[2])+p2[1]*(point.xy[2]-p1[2])+point.xy[1]*(p1[2]-p2[2])
		if(fequal(0,detm)) {
			# then that point is on segment
			# is it vertical? if so choose y
			if ( fequal(p1[1],p2[1]) && fequal(point.xy[1],p1[1])) {
				s <- ( point.xy[2]-p1[2])/(p2[2]-p1[2])
			} else {
				s <- ( point.xy[1]-p1[1])/(p2[1]-p1[1])
			}
			smudge <-  (.Machine$double.eps ^ 0.5)			
			if (s>=0-smudge && s<= 1+smudge) {
				ison <- TRUE;
				break;
			}
		}
	}
	if (ison) { return(nix) } else { return(NA) }
}

setMethod(".checkPointOnEdge",c("VDedgeLines"),function(edge,point.xy) {
	return(!is.na(.find.point.on.EdgeLines(edge,point.xy)))
}
)


.face.midplace <- function(drawing,faceName) {
	# try forming the centroid of the points at the midplace of each edge
#	browser()
#	edgeClasses <- .faceEdgeClasses(drawing,faceName)
#	stopifnot( (length(edgeClasses)==2 & all(edgeClasses=="VDedgeSector")))
	edges <- .face.to.faceEdges(drawing,faceName)
	midpoints <- t(sapply(edges,.midpoint))
	midpoints.centroid <- matrix(apply(midpoints,2,mean),ncol=2,byrow=TRUE)
	midpoints.centroid
}

########################################
# if a point is on an edge, we will split want to the edge into two
#
setMethod(".splitEdgeAtPoint",c("VDedgeSector"),function(edge,point.xy) {
	new1 <- edge
	new2 <- edge
	pointName <- rownames(point.xy)
	new1@to <- pointName
	new2@from <- pointName
	pointTheta <- .point.xy.to.theta(point.xy,edge@centre)
#print(edge)
#print(pointTheta)
	if (edge@toTheta > 0) {
		stopifnot(edge@fromTheta>= pointTheta & pointTheta <= edge@fromTheta)
		new1@toTheta <- pointTheta
		new2@fromTheta <- pointTheta
	} else {
		if (pointTheta> edge@fromTheta) {
			new1@toTheta <- pointTheta - 2*pi
			new2@fromTheta <- pointTheta
			new2@toTheta <- edge@toTheta+ 2*pi
		} else {
			new1@toTheta <- pointTheta
			new2@fromTheta <- pointTheta
		}
	}
		new2@fromTheta <- pointTheta
#print(new1)
#print(new2)
	new1 <- .normalise.sector(new1)
	new2 <- .normalise.sector(new2)
	list(new1,new2)
})
	

setMethod(".splitEdgeAtPoint",c("VDedgeLines"),function(edge,point.xy) {
	nix <- .find.point.on.EdgeLines(edge,point.xy)
	if (is.na(nix)){ stop(sprintf("%s is not in edge",rownames(point.xy))) }
	new1 <- edge
	new2 <- edge
	pointName <- rownames(point.xy)
	new1@to <- pointName
	new2@from <- pointName

	new1@xy <- rbind(new1@xy[1:nix,],as.numeric(point.xy))
	n2rest <- new2@xy[(nix+1):nrow(new2@xy),,drop=FALSE]
#show(n2rest)
	if (all(fequal(as.numeric(point.xy),n2rest[1,]))) {
		new2@xy <- n2rest
	} else {
		new2@xy <- rbind(as.numeric(point.xy),new2@xy[(nix+1):nrow(new2@xy),])
	}
#show(edge)
#show(pointName)
#show(new1@xy)
#show(new2@xy)
	list(new1,new2)
})


spliceinstead <- function(vec,old,new) {
		medge <- match(old,vec)
		if (all(is.na(medge))) {
			return(vec)
		}
		upto <- if (medge[1]>1) { vec[1:(medge[1]-1)]} else { character(0) }
		endm <- medge[length(medge)]
		after <- if (endm  < length(vec)) { vec[(endm +1):length(vec)] } else { character(0) }
		 c(upto,new,after)
	}



############################
# when we add a new edge, we want to know which face it is in,
# which we do by finding a 'midpoint' on the edge and then locating which face that is in

# now it gets tricky
setMethod(".midpoint",c("VDedgeLines"),function(edge){
	edgexy <- .edge.to.xy(edge)
	if (nrow(edgexy)%%2 == 1) {
		midn <- (nrow(edgexy)+1)/2 
		midn <- c(midn,midn) 
	}	else {
		midn <- nrow(edgexy)/2
		midn <- c(midn,midn+1)
	}
	midx <- edgexy[midn,]
	midmean <- matrix(apply(midx,2,mean),ncol=2)
	midmean
})
setMethod(".midpoint",c("VDedgeSector"),function(edge){
	theta <- (edge@fromTheta+edge@toTheta)/2
	point.xy <- .theta.to.point.xy(theta,r=edge@radius,centre=edge@centre)
	point.xy
})

##############################################################################
##
##
##  faceList
##
##
#########################################################################

setClass("TDEdgeList",representation(edgeList="list"))
setClass("TDFaceList",representation(faceList="list",faceSignature="list")) # of the same length; the name of the face is its name in the list

.faceNames <- function(drawing,onlyVisible=FALSE) {
	faceNames <- names(drawing@faceList )
	if (onlyVisible) {
		faceNames <- setdiff(faceNames,"DarkMatter")
	}
	faceNames
}
.faceSignatures <- function(drawing,onlyVisible=FALSE) {
	faceSignatures <- drawing@faceSignature 
	if (onlyVisible) {
		faceSignatures[["DarkMatter"]] <- NULL
	}
	faceSignatures 
}
	

updateSignature <- function(drawing,faceNames,suffix) {
	for (faceName in faceNames) {
		sig <- drawing@faceSignature[[faceName]]
		if (sig=="DarkMatter" & suffix =="1") {
			sig <-  paste(c(rep("0",length(drawing@setList)-1),"1"),collapse="")
		} else {
			sig <- paste(sig,suffix,sep="")
		}
		drawing@faceSignature[[faceName]] <- sig 
	}
	drawing
}
setSignature <- function(drawing,faceName,signature) {
	drawing@faceSignature[[faceName]] <- signature
	drawing
}



.faceEdgeNames <- function(drawing,faceName,unsigned=FALSE,type="face") { 
	if (type=="set") {
		edges <- drawing@setList[[faceName]]
	} else {
		edges <- drawing@faceList[[faceName]] 
	}
	if (unsigned) { edges <- sub("^-","",edges);  }
	
	edges
}

.faceEdgeClasses <- function(drawing,faceName,type="face") { 
	if (type=="set") {
		edges <- drawing@setList[[faceName]]
	} else {
		edges <- drawing@faceList[[faceName]] 
	}
	unsigned.edges <- sub("^-","",edges); 
	res <- sapply(drawing@edgeList[unsigned.edges],function(x)as.character(class(x)))
	res
}


renameFaces <- function(drawing,oldName,newName) {
	stopifnot(length(oldName) == length(newName)) 
	oldFaceNames <- names(drawing@faceList) 
	names(oldFaceNames) <- oldFaceNames
	newFaceNames <- oldFaceNames
	newFaceNames[oldName] <- newName
	newix <- which(duplicated(newFaceNames))
	for (ix in seq_along(newix)) {
		six <- 1; 
		while(TRUE) {
			newName <- paste(newFaceNames[newix[ix]],six,sep="-")
			if (! newName %in% newFaceNames) {
				break
			}
			six <- six + 1
		}
		newFaceNames[newix[ix]] <- newName
	}
	
	for (i in seq_along(oldName)) {
		names(drawing@faceList) <- newFaceNames
		names(drawing@faceSignature) <- newFaceNames[oldFaceNames]
	}
	drawing
}


setMethod("show","TDFaceList",function(object){
	facedf <- do.call(rbind,lapply(object@faceList,function(face){
	data.frame(faces=paste(face,collapse=";"))}))
	print(facedf)
	sdf <- data.frame(sig=do.call(c, object@faceSignature))
	print(sdf)
})




spliceEdgeIntoFace <- function(drawing,faceName,edgeName,edgeNames,doReverse=TRUE) {
	revEdgeName <- sub("^--","",paste("-",edgeName,sep=""))
	revEdgeNames <- rev(sub("^--","",paste("-",edgeNames,sep="")))
	thisFaceList <- drawing@faceList[[faceName]]
	thisFaceList  <- spliceinstead(thisFaceList  ,edgeName,edgeNames)
	if (doReverse) {
		thisFaceList  <- spliceinstead(thisFaceList  ,revEdgeName,revEdgeNames)
	}
	drawing@faceList[[faceName]] <- thisFaceList
	drawing@recentChanges <- edgeNames
	drawing
}


.startFaceAtPoint <- function(drawing,faceName,from) {
	FacePoints <-  .points.of.face (drawing,faceName)
	fromix=match(from,FacePoints)
	if (is.na(fromix)) { stop(sprintf("%s is not a point in face %s",from,faceName)) }
	oldFace <- drawing@faceList[[faceName]]
	face <- c(oldFace[(fromix):length(oldFace)],oldFace[seq(length=fromix-1)])
	drawing@faceList[[faceName]] <- face
	drawing
}
	


addFace <- function(drawing,faceName,faceSignature,face,edit=FALSE) {
#cat(sprintf("Adding face %s\n",faceName))
	newfaceName <- faceName
	if (!edit & !faceName=="DarkMatter") { # a new face.. but we won't overwrite an existing one with the same name unless it is dark matter
		if (faceName %in% .faceNames(drawing)) {
			ix <- 1; 
			while(TRUE) {
				newfaceName <- paste(faceName,ix,sep="-")
				if (!newfaceName %in% .faceNames(drawing)) { break }
				ix <- ix+1
			}
#cat(sprintf("...changed to  %s\n",newfaceName))
		}
	} 
#if (newfaceName %in% .faceNames(drawing)) { cat(sprintf("Adding %s which already exists with edit %s\n",newfaceName,edit)) }
	drawing@faceList[[newfaceName ]] <- face
	drawing@faceSignature[[newfaceName ]] <- faceSignature
	list(drawing=drawing,faceName=newfaceName)
}

deleteFace <- function(drawing,faceName) {
	drawing@faceList[[faceName]] <- NULL
	drawing@faceSignature[[faceName]] <- NULL
	drawing
}

getFace <- function(drawing,faceName,reverse=FALSE) {
	res <-drawing@faceList[[faceName]] 
	if (reverse) {
		res <- rev(	sub("^--","",paste("-",res,sep="")))
	}
	res
		
}


#####################################################################################
setClass("TissueDrawing",representation(
	"TDEdgeList", # named list of VDedgeDrawns
	"TDFaceList", # named list of (edge,hand) pairs
	setList="list", # named list of (edge, hand) pairs
	nodeList="list", # named list of xy coordinates
	recentChanges="character" # used in gross point injection code
	)
)


# nb dont currently require that all edges in edge list are actually used in faceList

.validateFaces <- function(drawing) {
cat(sprintf("Validating a drawing on %d sets...",length(drawing@setList)))

	for (faceName in .faceNames(drawing)) {
		# check we have all the edges
		for (edge.name in .faceEdgeNames(drawing,faceName)) {
			if (! .edge.in.Drawing(drawing,edge.name)) {
				stop(sprintf("Face %s has unknown edge %s",faceName,edge.name))
			}
		}
		# then that points are visited in sequence
		thisFaceEdges<- .face.to.faceEdges(drawing,faceName)
		fromtomat <- sapply(thisFaceEdges,function(y){c(y@from,y@to)})
		for (ix in seq_len(ncol(fromtomat))) {
			if (ix==ncol(fromtomat)){jx<-1}else{jx<-ix+1}
			if (!(fromtomat[2,ix]==fromtomat[1,jx])	) {
				print(fromtomat)
				stop(sprintf("Face %s is wrong",faceName))
			}
		}

	}
cat(sprintf("...done\n"))
}

.validateDrawing <- function(drawing) {
	.validateFaces(drawing)
	for (set.name in names(drawing@setList )) {
		for (edge.name in drawing@setList [[set.name]]) {
			if (! .edge.in.Drawing(drawing,edge.name)) {
				stop(sprintf("Set %s has unknown edge %s",set.name,edge.name))
			}
		}
	}
	for (edge.name in names(drawing@edgeList)) {
		edge = drawing@edgeList[[edge.name]]
		if ( !(edge@from %in% names(drawing@nodeList))) {
			stop(sprintf("Edge %s has unknown from node %s",edge.name,edge@from))
		}
		if ( !(edge@to %in% names(drawing@nodeList))) {
			stop(sprintf("Edge %s has unknown to node %s",edge.name,edge@to))
		}
	}
	faceSignatures <- unlist(.faceSignatures(drawing))
	if (faceSignatures[names(faceSignatures)=="DarkMatter"] != "DarkMatter") {
		warning(sprintf("Dark matter face has sig %s\n"),faceSignatures[names(faceSignatures)=="DarkMatter"])
	}
	nsig <- faceSignatures[!faceSignatures=="DarkMatter"]
	issig <- sapply(nsig,function(s){regexpr("^[01]*$",s)>0})
	if (!all(issig)) {
		warning(sprintf("Faces %s have sigs %s\n",paste(names(nsig)[!issig],collapse=";"),paste((nsig)[!issig],collapse=";")))
	}
	siglen <- sapply(nsig,nchar)
	if (length(unique(siglen))>1) {
		warning(sprintf("sigs %s differ in length\n",paste(nsig,collapse=";")))
	}
	siglen <- unique(siglen)[1]
	if (any(duplicated(nsig))) {
		dupsig <- unique(nsig[duplicated(nsig)])
		sapply(dupsig,function(sig){
			cat(sprintf("sig %s duplicated in faces %s\n",sig,paste(names(nsig)[nsig==sig],collapse=";")))
		})
	}
	Indicator <- ifelse((data.matrix(do.call(expand.grid,lapply(seq(1,length=unique(siglen)),function(x){c(0,1)}))))==1,1,0)
	inn <- apply(Indicator ,1,paste,collapse="")
	inn[inn==paste(rep("0",siglen),collapse="")] <- "DarkMatter"
	signotindiag <- faceSignatures[!faceSignatures %in% inn ]
	if (length(signotindiag )>0) {
		warning(sprintf("Signatures %s not present in diagram",paste(signotindiag ,collapse=";")))
	}
}

##############
# show method for debugging

setMethod("show","TissueDrawing",function(object){
	edgedf <- do.call(rbind,lapply(object@edgeList,
		function(edge) {
			df <- data.frame(from=edge@from,to=edge@to,type=as.character(class(edge)))

			if (inherits(edge,"VDedgeSector")) {
				df <- cbind(df,data.frame(npoints=NA,centre=paste(edge@centre,collapse=","),hand=edge@hand))
			} else {
				df <- cbind(df,data.frame(npoints=nrow(edge@xy),centre=NA,hand=NA))
			}
		}
	))
	print(edgedf)
	nodedf <- do.call(rbind,lapply(object@nodeList,function(node){
		dimnames(node)<- NULL; df <- data.frame(node); }))
	print(nodedf)
	show(as(object,"TDFaceList"))
	setdf <- do.call(rbind,lapply(object@setList,function(face){
		data.frame(paste(face,collapse=";"))}))
	print(setdf)
	
})


####################################################################
# plotting functions



.VDPlotArcs <- function(drawing,arcnames,arrowfunc=NULL,gp=gpar()) {
	drawing<-  as(drawing,"TissueDrawing")
	edgeList <- drawing@edgeList
	if (missing(arcnames)) {
		arcnames <- names(edgeList) 
	}
	for (arcname in arcnames) {
		if (arcname %in% names(edgeList))
			arcxy <- .edge.to.xy(edgeList[[arcname]])
		else {
			revEdge <- sub("^--","",paste("-",arcname,sep="")) 
			if (revEdge %in% names(edgeList) ) {
				arcxy <- .edge.to.xy(.reverseEdge(edgeList[[revEdge ]]))
			} else { 
				stop(sprintf("Can't find edge %s in drawing\n",arcname))
			}
		}
		grid.lines(arcxy[,1],arcxy[,2],default.units="native",arrow=arrowfunc,gp=gp)
	}
}
	


.PlotSetBoundaries.TissueDrawing <- function(drawing,gp){
#browser()
	drawing <- as(drawing,"TissueDrawing")
	setList <- drawing@setList	
	setNames <- names(drawing@setList)
	if (missing(gp)) {
		gp <- SetColours(drawing)
	}
	for (setName in setNames) {
		.VDPlotArcs(drawing,setList[[setName]],gp=gp[[setName]])
	}
}

setMethod("PlotSetBoundaries","TissueDrawing",.PlotSetBoundaries.TissueDrawing)


setMethod("PlotNodes","TissueDrawing",function(drawing,gp){
	dv <- as(drawing,"TissueDrawing")
	pxy <- dv@nodeList
	xy <- do.call(rbind,pxy)
	grid.text(x=xy[,1],y=xy[,2],label=names(pxy),default.units="native",
		just=c("left","bottom"))
})

.face.to.faceEdges <-  function(drawing,faceName,type="face") {
	if (type=="set") {
		faceEdgeNames <- drawing@setList[[faceName]]	
	} else {
		faceEdgeNames <- .faceEdgeNames(drawing,faceName)
	}
	faceSign <- substring(faceEdgeNames,1,1)=="-"
	unsignedEdgeNames <- sub("^-","",faceEdgeNames)
	notinedgeList <- unsignedEdgeNames [! unsignedEdgeNames %in% names(drawing@edgeList)] 
	if (length(notinedgeList)>0) {
		stop(sprintf("Face %s has unknown edges %s",faceName,paste(notinedgeList, collapse=" ")))
	}
	faceEdges <- drawing@edgeList[unsignedEdgeNames ]
	names(faceEdges) <- faceEdgeNames
	faceEdges[faceSign] <- lapply(faceEdges[faceSign],.reverseEdge)
	faceEdges
}

.face.toxy <- function(drawing,faceName,dx=0.05,type="face") {
	faceEdges <- .face.to.faceEdges(drawing,faceName,type=type)
	face.Sxy <- lapply(faceEdges,function(x).edge.to.xy(x,dx=dx))
	face.Sxy <- lapply(face.Sxy,function(x)x[-nrow(x),,drop=FALSE])
	
	all.xy <- do.call(rbind,face.Sxy)
	all.xy
}


.face.area <- function(drawing,faceName) {
	all.xy <- .face.toxy(drawing,faceName);
	.polygon.area(all.xy)
}

.polygon.area <- function(xy) {
	# this area is negative for clockwise polygons, sigh
	xy1 <- xy; xy2 <- xy[ c(2:nrow(xy),1),]
	x1 <- xy1[,1];y1 <-xy1[,2];x2<-xy2[,1];y2<- xy2[,2]
	det <- x1*y2 - x2*y1
	(sum(det)/2)

}

faceAreas <- function(drawing) {
	sapply(.faceNames(drawing),function(faceName)abs(.face.area(drawing,faceName)))
}

.polygon.centroid <- function(all.xy) {
	xy1 <- all.xy; xy2 <- all.xy[ c(2:nrow(all.xy),1),]
	x1 <- xy1[,1];y1 <-xy1[,2];x2<-xy2[,1];y2<- xy2[,2]
	area <- .polygon.area(all.xy)
	if (area>0) {
	cx <- sum((x1+x2) * ( x1 * y2 - x2 * y1))/(6* area)
	cy <- sum((y1+y2) * ( x1 * y2 - x2 * y1))/(6* area)
	} else {
		cx <- mean(x1);cy <- mean(y1)
	}
	centroid.xy <- matrix(c(cx,cy),ncol=2)
	centroid.xy
}

.face.centroid <- function(drawing,faceName) {
	all.xy <- .face.toxy(drawing,faceName)
	res <- .polygon.centroid(all.xy)
	names(res) <- "centroid"
	res
}

.PlotFace.TissueDrawing <- function(drawing,faceName,dx=0.05,gp=gpar(),doDarkMatter=FALSE) {
#cat(faceName,"\n")
	if (!doDarkMatter & faceName=="DarkMatter") {
		return()
	} 
	all.xy <- .face.toxy(drawing,faceName,dx=dx)
	#  first we do the fill using the gps
	if (!is.null(all.xy)) {
		grid.polygon(all.xy[,1],all.xy[,2],default.units="native",gp=gp)
	}
	if (!is.null(gp$lty)) { # have a line spec for the outside which we no longer index by set number
		faceEdges <- .face.to.faceEdges(drawing,faceName)
		face.Sxy <- lapply(faceEdges,.edge.to.xy,dx=dx)

		lapply(1:length(face.Sxy),function(edgeix) {
			faceEdge <-  faceEdges[[edgeix]]
			face.xy <- .edge.to.xy(faceEdge)
			if (nrow(face.xy)>0) {
				grid.lines(face.xy[,1],face.xy[,2],
					default.units="native",gp=gp)
			}
		})
	}
}

.PlotFaceNames.TissueDrawing <- function(drawing,faceNames,signature=TRUE){
	if(missing(faceNames)) {	
		faceNames <- .faceNames(drawing,onlyVisible=TRUE)
	}
	if (signature) {
		faceSignatures <- .faceSignatures(drawing,onlyVisible=TRUE)
		label <- (faceSignatures[faceNames])
	} else {
		label <- faceNames
	}
	for (faceName in faceNames) {
		text.xy <- .find.point.within.face (drawing,faceName)
#		grid.points(x=text.xy[,1],y=text.xy[,2],default.units="native")
		grid.text(x=text.xy[,1],y=text.xy[,2],label=label[[faceName]],default.units="native")
	}
}


.PlotFaces.TissueDrawing<- function(drawing,faceNames,gp,arrow,colourAlgorithm){
	if(missing(faceNames)) {	
		faceNames <- .faceNames(drawing)
	}
	if (missing(gp)){ 
		gp <- FaceColours(drawing=drawing,faceNames=faceNames,colourAlgorithm=colourAlgorithm)
	}
	for (face.ix in seq_along(faceNames)) {
		faceName <- faceNames[face.ix]
		.PlotFace.TissueDrawing(drawing,faceName,gp=gp[[faceName]])
	}
}

setMethod("PlotFaces","TissueDrawing",.PlotFaces.TissueDrawing)

.find.point.on.face <- function(drawing,faceName) {
	edgeName <- .faceEdgeNames(drawing,faceName)[1]
	edgeName <- sub("^-","",edgeName)
	edge <- drawing@edgeList[[edgeName]]
	point <- .midpoint(edge)
	point
}

.face.maxradius <- function(drawing,faceName) {
	edgebb <- lapply(.face.to.faceEdges(drawing,faceName),function(x)x@bb)
	absbb <- max(sapply(edgebb,function(x)max(abs(x))))
	maxradius <- sqrt(2)*absbb
	maxradius 
}


.find.point.within.face <- function(drawing,faceName,treat.dark.matter.as.face=FALSE) {
	if (faceName=="DarkMatter" & treat.dark.matter.as.face) {
		faceName <- ".fpwf"
		drawing <- renameFaces(drawing,"DarkMatter",faceName)
		# has the effect of treating as an ordinary face
	}
	edgeClasses <- .faceEdgeClasses(drawing,faceName)
	if (all(edgeClasses=="VDedgeSector") & length(edgeClasses)==2) {
		aPoint <- 	.face.midplace(drawing,faceName)
		if (.is.point.within.face(drawing,faceName,aPoint )) {
			return(aPoint )
		} else {
			aPoint <- .face.centroid(drawing,faceName=faceName)
			if (.is.point.within.face(drawing,faceName,aPoint )) {
				return(aPoint )
			}
		} 
	} else { #otherwise test in other order
		aPoint <- 	.face.centroid(drawing,faceName=faceName)
		if (.is.point.within.face(drawing,faceName,aPoint )) {
			return(aPoint )
		} else {
			aPoint <- .face.midplace(drawing,faceName=faceName);
			if (.is.point.within.face(drawing,faceName,aPoint )) {
				return(aPoint )
			}
		} 
	}

	# ok, just try and find some corner guaranteed to be inside

	ear.triangle <- .find.triangle.within.face(drawing,faceName)
	earCentroid <- .polygon.centroid(ear.triangle)
	if (!	.is.point.within.face(drawing,faceName,earCentroid )) {
		stop(sprintf("Ear method failed in face %s\n",faceName))
	}
	return(earCentroid)
}

.find.triangle.within.face <- function(drawing,faceName) {
	# poor mans triangulation... subtracting ear method cf wikipedia polygon triangulation
	xy <- .face.toxy(drawing,faceName)
	A <- .face.area(drawing,faceName)
	if (A==0) { # collinear
		return(xy[1:3,])
	}
	xy <- rbind(xy,xy[1:2,])
	fix <- NA
	for (ix in 2:(nrow(xy)-1)) {
		from <- xy[ix-1,,drop=FALSE]
		to <- xy[ix+1,,drop=FALSE]
		pt <- xy[ix,,drop=FALSE]
		thetafrom <- atan2( from[,2]-pt[,2],from[,1]-pt[,1])
		thetato <-   atan2(   to[,2]-pt[,2],to[,1]-pt[,1])
		thetato <- thetato - thetafrom
		thetato <- thetato %% (2 * pi)
		if (thetato > pi) { # not a convex point
			next
		}
		npoints <- .probe.chord.intersections(drawing,faceName,from,to)
		fromdist <- ((npoints[,1]-from[1])^2+(npoints[,2]-from[2])^2 ) * sign(npoints[,1]-from[1])
		npoints <- npoints[order(fromdist),]; fromdist <- sort(fromdist)
		fromix <- min(which(fequal(fromdist,0)))
		if (fromix !=1) {
			fromdist <- -fromdist
			npoints <- npoints[order(fromdist),];fromdist <- sort(fromdist)
			fromix <- min(which(fequal(fromdist,0)))
			stopifnot(fromix==1)
		}
		# the next point along has to be the to point otherwise intersection
		nextix <- min(which(!fequal(fromdist,0)))
		nextpt <- npoints[nextix,,drop=FALSE]
		nointersect <- all(fequal(nextpt,to))
		if (nointersect) {
			fix <- ix
			break
		}
	}
	if (is.na(fix)) {
		stop(sprintf("Can't find ears for face %s\n",faceName))
	}
	return(xy[ (fix-1):(fix+1),])
}


internalPointsofFaces <- function(drawing) {
	fNames <-setdiff(.faceNames(drawing),"DarkMatter") 
	res <- lapply(fNames ,function(x).find.point.within.face(drawing=drawing,x))
	res <- do.call(rbind,res)
	res <- rbind(res,c(NA,NA))
	rownames(res) <- c(fNames,"DarkMatter")
	
	res
}

#############
.edge.in.Drawing <- function(drawing,edge.name) {
	unsignedEdgeName <- sub("^-","",edge.name)
	unsignedEdgeName %in% names(drawing@edgeList)
}

######################################################
# we will build up drawings by adding edges and points:



injectPoint <- function(drawing,edgeName,newPoint) {
#print(edgeName);

	pointName <- rownames(newPoint)
	if (pointName %in% names(drawing@nodeList)) {
		existingPoint <- drawing@nodeList[[pointName]]
		if (!all(fequal(newPoint-existingPoint,0))) {
			stop(sprintf("Point %s already exists at a different place",pointName))
		}
	} else {
		drawing@nodeList[[pointName]] <- newPoint
	}

	if (is.null(edgeName)) {
		# just an isolated point
		return()
	}

	thisEdge <- drawing@edgeList[[edgeName]]
	if (is.null(thisEdge)) {stop(sprintf("Edge %s not found",edgeName))	}
	if (!.checkPointOnEdge(edge=thisEdge,point.xy=newPoint)) { 
		show(drawing)
		stop(sprintf("Point %s not on edge %s",pointName,edgeName))
}

	new12 <- .splitEdgeAtPoint(thisEdge,newPoint)
	splitName <- strsplit(edgeName,split="|",fixed=TRUE)[[1]]
	if (length(splitName>= 3)) { setName <- splitName[3] } else {setName <- "?"}

	edgeNames <- c(paste(thisEdge@from,pointName,setName ,sep="|"),paste(pointName,thisEdge@to,setName ,sep="|"))
	# remove the old one
#cat(sprintf("Replacing edge %s by %s and %s\n",edgeName,edgeNames[1],edgeNames[2]))
	drawing@edgeList[[edgeName]] <- NULL
	drawing@edgeList[[edgeNames[1]]] <- new12[[1]]
	drawing@edgeList[[edgeNames[2]]] <- new12[[2]]

	# now insert into the setlist

	for (setName in names(drawing@setList)) {
		if (edgeName %in% drawing@setList[[setName]]) {
			drawing@setList[[setName]] <- spliceinstead(drawing@setList[[setName]],edgeName,edgeNames)
		}
	}

	# and now into the facelist
	for (faceName in .faceNames(drawing)) {
		drawing <- spliceEdgeIntoFace(drawing,faceName,edgeName,edgeNames)
	}

	drawing	
}

##########################################
# code for detecting if a point is within a face
#################
# W Randolph Franklin's PNPOLY
# http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
# int pnpoly(int npol, float *xp, float *yp, float x, float y)
#    {
#      int i, j, c = 0;
#      for (i = 0, j = npol-1; i < npol; j = i++) {
#        if ((((yp[i]<=y) && (y<yp[j])) ||
#             ((yp[j]<=y) && (y<yp[i]))) &&
#            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
#
#          c = !c;
#      }
#      return c;		
	

pnpoly <- function(xp,yp,x,y) {
	npol <- length(xp); stopifnot(npol==length(yp))
	c= 0
	for (i in seq_len(npol)) {
		if (i==1) { j<- npol } else {j<-i-1}
# cat(i," ",j,"\n")
       if ((((yp[i]<=y) && (y<yp[j])) ||
             ((yp[j]<=y) && (y<yp[i]))) &&
	           (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i])) {
		  c <- c+1
		}
	}
	return ( c %%2 == 1)
}



pnpolytest <- function() {
	pnpoly( xp=c(-1,1,1,-1),yp=c(-1,-1,1,1),x=0,y=0)
	pnpoly( xp=c(-1,1,1,1,-1),yp=c(-1,-1,0,1,1),x=0,y=0)
	pnpoly( xp=c(-1,1,1,1,1,-1),yp=c(-1,-1,0,0,1,1),x=0,y=0)
	pnpoly( xp=c(1,1,1,-1,-1,1),yp=c(0,0,1,1,-1,-1),x=0,y=0)
	pnpoly( xp=c(1,1,-1,-1,1,1),yp=c(0,1,1,-1,-1,0),x=0,y=0)
	pnpoly( xp=c(-3,-2,-3,-3),yp=c(2,1,0,2),x=-4.5,y=0)
	pnpoly( xp=c(-3,-2,-3),yp=c(2,1,0),x=-4.5,y=0)
}


.is.point.within.face <- function(drawing,faceName,point.xy,type="face") {
	.face.xy <-  .face.toxy(drawing,faceName,type=type)
	inFace <- pnpoly(.face.xy[,1],.face.xy[,2],point.xy[,1],point.xy[,2])
	if (faceName=="DarkMatter") {
		inFace <- !inFace
	}
	return(inFace)
}

.is.face.within.set <- function(drawing,faceName,setName) {
		aPoint <- .find.point.within.face(drawing,faceName)
		.is.point.within.face(drawing,faceName=setName,point.xy=aPoint,type="set")

}


.points.of.face <- function(drawing,faceName,from) {
	# these are returned in sequence, starting from (the first) FROM if specified
	thisFaceEdges<- .face.to.faceEdges(drawing,faceName)
	fromtomat <- sapply(thisFaceEdges,function(y){c(y@from,y@to)})
	thisFacePoints <- as.vector(c(fromtomat[1,1],fromtomat[2,-ncol(fromtomat)]))
	if (!missing(from)) {
		thisFacePoints2 = c(thisFacePoints,thisFacePoints)
		mix = match(from,thisFacePoints2)
		if (is.na(mix)){stop(sprintf("%s is not in face %s",from,faceName))}
		thisFacePoints <- thisFacePoints2[ (mix): (mix+ length(thisFacePoints)-1)]
	}
	thisFacePoints
}

.find.face.containing.edge <- function(drawing,edgeList) {
# normally faces _dont_ contain edges but when we are injecting a new one
# we want to find which face it is crossing

# a new edge must be inside one of the faces which have from and to as members 
# or outside all of them (equivalently in the dark matter)
	from= edgeList[[1]]@from
	to = edgeList[[length(edgeList)]]@to
	pointsPerFace <- lapply(.faceNames(drawing),function(faceName){
		.points.of.face(drawing,faceName)
	})
	names(pointsPerFace) <- .faceNames(drawing)
	fromtoInFace <- sapply(pointsPerFace ,function(x)from %in% x & to %in% x)
	fromtoFaceNames <- names(fromtoInFace)[fromtoInFace]
	
	inFaceName <- "DarkMatter"
	edge.midpoint <- .midpoint(edgeList[[1]])
	for (faceName in setdiff(fromtoFaceNames,"DarkMatter")) {
		if (.is.point.within.face(drawing,faceName,edge.midpoint)) {
			inFaceName <- faceName
			break
		}
	}
	inFaceName
}

injectEdge <- function(drawing,newEdgeList,inFaceName,set2Name,addToList=TRUE) {
	if(addToList) {	drawing@edgeList <- c(drawing@edgeList,newEdgeList) }

	inFaceName <- .find.face.containing.edge(drawing=drawing,edgeList=newEdgeList)
	edgeNames <- names(newEdgeList); 
	edge.from <- newEdgeList[[1]]@from
	edge.to <- newEdgeList[[length(newEdgeList)]]@to
#cat(inFaceName,"\n")
	drawing <- .startFaceAtPoint(drawing,faceName=inFaceName,from=edge.from)
	oldFacePoints <- .points.of.face (drawing,faceName=inFaceName)

	# 
	fromix.possible = which(edge.from==oldFacePoints)
	toix.possible = which(edge.to==oldFacePoints)
	if (length(toix.possible)==0) { browser();stop(sprintf("%s is not in face %s",edge.to,inFaceName)) }


	if (length(fromix.possible)==1) {
		fromix <- fromix.possible
	} else  {
		warning(sprintf("Face %s has points %s; can't cope well with joining repeated FROM point %s yet",
		inFaceName,paste(oldFacePoints,collapse=","),edge.from))
		fromix <- fromix.possible[1]
	# we know we want to split the face between FROM and TO, but
	# the trouble is that eg from may be repeated in the edge of the
	# face and we need to pick the right one
	# first find the right FROM ; take the last possible TO
	}

	if (length(toix.possible)==1) {
		toix<- toix.possible
	} else  {
		warning(sprintf("Face %s has points %s; can't cope well with joining repeated TO point %s yet",
		inFaceName,paste(oldFacePoints,collapse=","),edge.to))
		toix<- toix.possible[length(toix.possible)]
	}

	inFaceSignature <- .faceSignatures(drawing)[[inFaceName]]
	nSets <- length(drawing@setList)
	if (inFaceName=="DarkMatter") {
		newFaceSigP <- paste(c(rep("0",nSets-1),"1"),collapse="")
		newFaceSigM <- "DarkMatter"
	} else {
		if (nchar(inFaceSignature)!=nSets) {
			newFaceSigP <- paste(inFaceSignature ,"1",sep="")
			newFaceSigM <- paste(inFaceSignature ,"0",sep="")
		} else {
			# that means this face has already been subdivided,
			# so further subdivision takes us back out of the Set
			newFaceSigP <- inFaceSignature 
			newFaceSigM <- paste(substr(inFaceSignature ,1,nSets-1),"0",sep="")
	#		if (newFaceSigM==paste(rep("0",nSets),collapse="")) {
	#			newFaceSigM <- "DarkMatter"
	#		}
		}
	
	}
	newFaceNamep <- newFaceSigP; newFaceNamem <- newFaceSigM

	# we duplicate so we can wrap around 
	oldFace <- .faceEdgeNames(drawing,inFaceName)
	oldFace2 <- c(oldFace,oldFace)
	reverseEdgeNames <- rev(sub("--","",paste("-",edgeNames,sep="")))


	# using the forward new edge will cut out FROM to TO and replace it by our new edgeNames
	# and then from the old face TO, back round the long way, to FROM
	# the only time that is wrong is if toix==fromix and  the new set completely encircles the old one
	if (toix==fromix) {
		new1 <- edgeNames
	}	else {
		f1start <- toix
		f1end <- fromix+length(oldFacePoints)
		new1 <- c(edgeNames,oldFace2[ f1start:(f1end-1)])
	}
	# that will be the new +ve entry
	res <- addFace(drawing,newFaceNamep,newFaceSigP,new1)
	drawing <- res$drawing; newFaceNamep <- res$faceName
	# the exception is when the new edge completely encircles the old one (and so fromix==toix as well)
	isEncircled <- FALSE
	if (fromix==toix) {
		oldPoint <- .find.point.on.face(drawing,inFaceName)
		oldinnew <- .is.point.within.face(drawing,faceName=newFaceNamep,point.xy=oldPoint)
		isEncircled <-  oldinnew 
	}
	if (isEncircled) {
		newFace <- c(new1,oldFace)
		drawing <- addFace(drawing,newFaceNamep,newFaceSigP,newFace,edit=TRUE)$drawing
	}
			
	# using the reverse edge keeps FROM to TO, then returns with the new edge	
	if(!isEncircled) {
		if (fromix==toix) {
			f2start <- fromix
			f2end <- fromix+length(oldFace)-1
		} else {
			f2start <- fromix
			f2end <- toix - 1
		}
		stopifnot(f2start*f2end!=0)
		new2 = c(oldFace2[ f2start:f2end],reverseEdgeNames)
	} else {
		new2 <- reverseEdgeNames
	}
	res <- addFace(drawing,newFaceNamem,newFaceSigM ,new2)
	drawing <- res$drawing; newFaceNamem <- res$faceName

	if (inFaceName != newFaceNamem & inFaceName != newFaceNamep) {
		drawing <- deleteFace(drawing,inFaceName)
	}


#cat(sprintf("Split face %s into %s and %s\n",inFaceName,newFaceNamep ,newFaceNamem))
	drawing
}


##############################################
# Make a drawing from a circle
newTissueFromCircle <- function(centre.xy,radius,Set=1,nodes=1) {
	nodeangles <-   2 * pi -  ( 0:(nodes-1) * 2 * pi / (nodes))
	p1.xy <- radius * cbind(cos(nodeangles),sin(nodeangles))
	p1.xy <- rep(centre.xy,each=nrow(p1.xy))+p1.xy

	rownames(p1.xy) <- paste("c",Set,1:nrow(p1.xy),sep="")

	nodeName <- rownames(p1.xy)[1] 
	lhcircle <- newEdgeSector(from=nodeName,to=nodeName,radius=radius,fromTheta=2*pi,toTheta=0,centre=centre.xy,hand=1)
	edgeName <- paste(nodeName,nodeName,Set,sep="|")
	edgeList <- list(lhcircle); names(edgeList) <- edgeName
	faceName <- "1"
	faceList <- list(edgeName,paste("-",edgeName,sep=""));names(faceList)<-c(faceName,"DarkMatter")
	faceSignature <- lapply(names(faceList),function(x){x}); names(faceSignature) <- names(faceList)
	nodeList <- list(p1.xy[1,,drop=FALSE])
	names(nodeList) <- rownames(p1.xy)[1]
	setList <- faceList[1]; names(setList) <- paste("Set",Set,sep="")
	VD <- new("TissueDrawing",edgeList=edgeList,faceList=faceList,faceSignature=faceSignature ,
		nodeList=nodeList,setList=setList)
	if  (nrow(p1.xy)>1) {
		VD <- injectPoints(VD,edgeName,p1.xy[-1,,drop=FALSE])
	}
	VD
}

newTissueFromEllipse <- function(f1,phi,e,a,Set,dx=0.05) {
	arclength <- (4*pi) 
	nintervals <- arclength/dx

		twoc <- a* e* 2
		f2 <- f1+  twoc*c(-cos(phi),sin(phi))
		theta <- seq(0,2*pi,length=nintervals)
		r <- (a * (1-e^2))/(1+e*cos(theta+phi))
		x <- f1[1]+r*cos(theta)
		y <- f1[2]+r*sin(theta)
		
	points.xy <- cbind(x,y)	
	rownames(points.xy) <- paste("e",letters[Set],1:nrow(points.xy),sep="")
	newTissueFromPolygon(points.xy=points.xy,Set=Set)
}

##############################################
# Make a drawing from a rectangle
newTissueFromPolygon <- function(points.xy,Set=1) {
	points.xy <- 	.removeDuplicates (points.xy)
	if (nrow(points.xy)<3) {stop("Not enough distince points for a polygon")}
	isclockwise <- (.polygon.area(points.xy) < 0) 
	if (!isclockwise) {
		points.xy <- points.xy[nrow(points.xy):1,]
	}

	frompoint <- points.xy[1,,drop=FALSE]
	if (!is.null(rownames(frompoint))) {
		from = rownames(frompoint) ; 
	} else {
		from <- paste("p",letters[Set],"1",sep="")
	}
	to = from;
	points.xy <- rbind(points.xy,frompoint )
	ledge<- newEdgeLines(from=from,to=from,xy=points.xy)
	edgeName <- paste(from,to,Set,sep="|")
	edgeList <- list(ledge); names(edgeList) <- edgeName
	faceName <- "1"; 
	faceList <- list(edgeName,paste("-",edgeName,sep=""));names(faceList)<-c(faceName,"DarkMatter")
	nodeList <- list(frompoint); names(nodeList) <- from
	faceSignature <- lapply(names(faceList),function(x){x}); names(faceSignature) <- names(faceList)
	setList <- faceList[1]; names(setList) <- paste("Set",Set,sep="")
	VD <- new("TissueDrawing",edgeList=edgeList,faceList=faceList,faceSignature=faceSignature ,
		nodeList=nodeList,setList=setList)
	VD
}


injectPoints <- function(drawing,edgeName,newPoints) {
	if (nrow(newPoints)<1) { return(drawing) }
	newPoint <- newPoints[1,,drop=FALSE]
	drawing <- injectPoint(drawing=drawing ,edgeName=edgeName,newPoint=newPoint) 
	newEdges <- drawing@recentChanges
	otherPoints <- newPoints[-1,,drop=FALSE]
	if(nrow(otherPoints)>0) {

		arepointsOn2 <- sapply(seq_len(nrow(otherPoints)) ,function(pointix){
			.checkPointOnEdge(edge=drawing@edgeList[[newEdges[2]]],otherPoints[pointix,,drop=FALSE])})
		pointsOn2 <- otherPoints[arepointsOn2,,drop=FALSE]
		if (nrow(pointsOn2)>0) {
			drawing <- injectPoints(drawing,newEdges[2],pointsOn2 )
		}

		arepointsOn1 <- sapply(seq_len(nrow(otherPoints)) ,function(pointix){
			.checkPointOnEdge(edge=drawing@edgeList[[newEdges[1]]],otherPoints[pointix,,drop=FALSE])})
		pointsOn1 <- otherPoints[arepointsOn1,,drop=FALSE]
		if (nrow(pointsOn1)>0) {
			drawing <- injectPoints(drawing,newEdges[1],pointsOn1 )
		}
	}

	drawing
		
}

##########################################
# combine two drawings...first by adding a face of the second
	


addSetToDrawing <- function(drawing1,drawing2,set2Name,remove.points=FALSE) {
		

	# ensure no clash of face names
	face2Name <- .faceNames(drawing2,onlyVisible=TRUE)
	stopifnot(length(face2Name )==1 )

	tempface2Name <- ".aSTD.temp";
	drawing2 <- renameFaces(drawing2,face2Name,tempface2Name)

	# check for name clashes, and also for identical points with different names
	new2 <- .check.clashes(drawing1,drawing2)

	# calculate a list of all the intersections and add them in as points
	res <- .add.intersection.points(drawing1=drawing1,drawing2=new2)
	new1=res$new1;new2=res$new2;intersectionPoints=res$intersectionPoints

	# we may have two edges with different names but identical geometry, if so use the name from diagram1
	res <- .check.duplicate.edges(drawing1=new1,drawing2=new2)
	new1<- res$new1; new2 <- res$new2

	# edges only in diagram2 need to be added to diagram1
	posnames <- sub("^-","",names(new2@edgeList)); negnames <- paste("-",posnames,sep="")
	newEdges <- new2@edgeList[! (posnames %in% names(new1@edgeList) | negnames %in% names(new1@edgeList)) ]
	new1@edgeList <- c(new1@edgeList,newEdges)
	# for lookups below, need to do the same in diagram2 in case we relabelled an edge that a face relies on
	new2@edgeList <- c(new2@edgeList,new1@edgeList[!names(new1@edgeList) %in% names(new2@edgeList)])


	# that will have updated the single set in new2
	new1@setList <- c(new1@setList,new2@setList)
	
	newNodes <- setdiff(names(new2@nodeList),names(new1@nodeList))
	new1@nodeList <- c(new1@nodeList,new2@nodeList[newNodes])


	# now split boundary of face2 into its edges
	face2Points <- .points.of.face(new2,tempface2Name )
	# are any of them intersection points?
	face2IntersectionPoints <- intersect(face2Points,intersectionPoints)

	if (length(face2IntersectionPoints)==0) {
		if (length(newEdges)==0) {
			new1 <- .addSetWithExistingEdges(new1,drawing2,tempface2Name)
		} else {
			new1 <- .addNonintersectingFace(new1,drawing2,tempface2Name) 
		}
		if (!inherits(new1,"TissueDrawing")) { return(NA) } # when we can't build an invisible edge joining the two faces
	} else { 
		new1 <- .addIntersectingFace(new1,new2,tempface2Name,face2IntersectionPoints)
	} 
	new1 <- .merge.faces.invisibly.split(new1)

	if (remove.points) { 
		new1 <- remove.nonintersectionpoints(new1)
	}
	new1


}


.addSetWithExistingEdges<- function(drawing1,drawing2,tempface2Name) {
	# the new Set is solely comprised of existing edges 
	# so there are no new faces
	# all we need to do is update the signatures
	# probably an easier way to do this, by tracking which edges came from which sets...
	# but we just cycle through the faces and update through brute force


	for (faceName in setdiff(.faceNames(drawing1),"DarkMatter")) {
		aPoint <- .find.point.within.face(drawing1,faceName)
		inFace <- .is.point.within.face(drawing=drawing1,faceName=tempface2Name,point.xy=aPoint,type="set")
		inFaceSig <- if (inFace) { "1"} else {"0"}
		drawing1<- updateSignature(drawing1,faceName,inFaceSig )
	}
	drawing1
}

.addNonintersectingFace <- function(drawing1,drawing2,tempface2Name) {
	# drawing2 contains a single face 
	#must be inside one of the faces or outside them all
	#either way can add the new face unchanged to the faceList
	
	res <- addFace(drawing=drawing1,faceName=tempface2Name,faceSignature="dummy",face=getFace(drawing2,tempface2Name))
	new1 <- res$drawing; tempface2Name<- res$faceName

	# then find which face it is within, first so we can set the signature correctly
	aPoint <- .find.point.within.face(drawing2,tempface2Name)
	outerFaceName <- ""
	for (faceName in .faceNames(new1)) {
		if(.is.point.within.face(drawing=new1,faceName=faceName,point.xy=aPoint)) {
			outerFaceName <- faceName
			break
		}
	}
	new1 <- setSignature(new1,tempface2Name,.faceSignatures(new1)[[outerFaceName]])

	# but then need to add invisible edges to join in to rest of drawing
	res <- .create.edge.joining.faces(drawing=new1,outerFaceName=outerFaceName ,innerFaceName=tempface2Name )
		if (!res$ok) { 
			return(NA)
		}
	iedgeName <- res$edgeName; new1 <- res$drawing;
	redgeName <- paste("-",iedgeName,sep="")
		
	# point to attach to (called outer because I imagined the face being inside, 
	# but also works when inserting a new face into dark matter and this point 
	outerPoint <-  new1@edgeList[[iedgeName ]]@from
	innerPoint <- new1@edgeList[[iedgeName ]]@to
	new1 <- .startFaceAtPoint(new1,tempface2Name,innerPoint)
	newEdges <- c(iedgeName,getFace(new1,tempface2Name,reverse=TRUE),redgeName)
	# find the edge going in to it
	outerEdges <- .face.to.faceEdges(new1,outerFaceName)
	edgetoPoint <- names(outerEdges)[min(which(sapply(outerEdges,function(edge)edge@to==outerPoint)))]
	# normally when replacing edges we want to do it for both faces containing the edge,
	#  but not in this case hence doReverse=FALSE
	new1 <- spliceEdgeIntoFace (drawing=new1,faceName=outerFaceName,edgeName=edgetoPoint,edgeNames=c(edgetoPoint,newEdges),doReverse=TRUE) 
	# now calculate the names
	oldFaceNames <- .faceNames(new1); faceNames <- oldFaceNames
	notInvolved <-  !faceNames %in% c(tempface2Name,"DarkMatter")
	faceNames[ notInvolved ] <- paste(faceNames[ notInvolved ],"0",sep="")
	if (outerFaceName=="DarkMatter") {
		face2Name <- paste(c(rep("0",length(new1@setList)-1),"1"),collapse="")
	} else {
		face2Name <- paste(outerFaceName,"1",sep="")
	}
	faceNames[ faceNames ==tempface2Name] <- face2Name
	new1 <- renameFaces(new1,oldFaceNames,faceNames)
	new1 <- updateSignature(new1,faceNames[notInvolved],"0")
	new1 <- updateSignature(new1,face2Name,"1")
	new1
}

.find.point.in.diagram <- function(drawing,aPoint) {
	xy <- do.call(rbind,drawing@nodeList)
	dist <- (xy[,1]-aPoint[1])^2 + (xy[,2]-aPoint[2])^2
	isEqual <- fequal(dist,0)
	if (!any(isEqual)) { return(NA)}
	if (length(which(isEqual))>1) stop("A third nonuique point")
	pointName <- rownames(xy)[isEqual]
	return(pointName)
}

.create.edge.joining.faces <- function(drawing,outerFaceName,innerFaceName) {
	# if outerFaceName is DarkMatter, then we really want to connect any
	# one of the other faces to innerFaceName, and the idea is to draw a line joining the centres of the two faces
	# that must have at least one segment which joins (something connected to the first face) to (something connected to) the second face
	# if outerFaceName is a regular face, (and then the innerFace is actually nested inside it, hence the names
	# then a point on its boundary will do as well
	# in practice the second face is always a single set though

	if (outerFaceName=="DarkMatter") {
		outerFaceForPoint <- setdiff(.faceNames(drawing),"DarkMatter")[1]
		outerPoint <- .find.point.within.face(drawing,outerFaceForPoint )
	} else {
		outerFaceForPoint <- outerFaceName
		outerPointName <- .points.of.face(drawing,outerFaceName)
		outerPoint <- drawing@nodeList[[outerPointName]]
	}
	innerPoint <- .find.point.within.face(drawing,innerFaceName)
	rownames(outerPoint) <- ".cejf"
	
	# find all the places it hits the 'inner' ie single set face 
	innerpoints <- .probe.chord.intersections(drawing,innerFaceName,outerPoint ,innerPoint )
	innerpoints <- innerpoints [rownames(innerpoints ) != ".cejf",,drop=FALSE]
	# and all the points it hits things connected to the outer face, so have to look through all the other faces too
	outerpoints <- matrix(nrow=0,ncol=2)
	for (faceName in setdiff(.faceNames(drawing),c("DarkMatter",innerFaceName))) {
		opoints <- 	.probe.chord.intersections(drawing,faceName ,outerPoint ,innerPoint )
		outerpoints <- rbind(outerpoints,opoints )
		outerpoints <- unique(outerpoints [rownames(outerpoints ) != ".cejf",,drop=FALSE])
	}
	# code the points by 1 or 2 depending on whether they are on the inner or outer faces
	linepoints <- rbind(cbind(outerpoints,1),cbind(innerpoints,2))
	dist <- (linepoints[,1]-outerPoint [1])^2+(linepoints[,2]-outerPoint [2])^2
	# sort by distance from the outer set
	linepoints <- linepoints[order(dist),,drop=FALSE]
	lastOuter <- max(which(linepoints[,3]==1))
	stopifnot(lastOuter < nrow(linepoints)) # because there should be at least one inner intersection
	# now we have a point that will work..it may already be in the diagram but if not must inject it
	op <- linepoints[lastOuter,1:2,drop=FALSE]
	opName <- .find.point.in.diagram(drawing,op)
	if (is.na(opName)) {
		opEdge <- strsplit(rownames(op),";")[[1]][1]
		nix <- .node.number.unused(drawing)
		opName <- paste("e",nix,sep="")
		rownames(op) <- opName 
		drawing <- injectPoint(drawing,opEdge,op)
	} else {
		rownames(op) <- opName
	}
	ip <- linepoints[lastOuter+1,1:2,drop=FALSE]
	ipName <- .find.point.in.diagram(drawing,ip)
	if (is.na(ipName)) {
		ipEdge <- strsplit(rownames(ip),";")[[1]][1]
		nix <- .node.number.unused(drawing)
		ipName <- paste("e",nix,sep="")
		rownames(ip) <- ipName 
		drawing <- injectPoint(drawing,ipEdge,ip)
	} else {
		rownames(ip) <- ipName
	}

	xy <- do.call(rbind,drawing@nodeList[c(opName,ipName)])
	testEdge <- newEdgeLines(from=opName,to=ipName,xy=xy,visible=FALSE)	
	stopifnot(!.internal.edge.drawing.intersection(drawing,testEdge)) 

	edgeName <- paste(opName,ipName,"invisible",sep="|")
	tel <- list(testEdge); names(tel) <- edgeName
	drawing@edgeList <- c(drawing@edgeList,tel)
	return(list(edgeName=edgeName,drawing=drawing,ok=TRUE))
	
}

.probe.chord.intersections <- function(drawing,faceName,chord.from.xy,chord.to.xy)  {
	# given two points, chord.from.xy outside the face, and a second point chord.to.xy,
	# draw a line between the two, and see where the line intersects the face. 
	# then arrange all of these intersection points in the order they appear in along the line,
	# including the chord.from.xy point

	

	# names pc and pmid not used
	chord <- newEdgeLines(from="pc",to="pmid",xy=rbind(chord.from.xy,chord.to.xy))

	foundList <- list()
	
	for ( edgeName in unique(.faceEdgeNames(drawing,faceName,unsigned=TRUE)) ) {
		faceEdge <- drawing@edgeList[[edgeName]]
		found <- .findIntersection(chord,faceEdge)
		foundList[[edgeName]] <- found
	}	
	# now we have a collection of points at which the chord crosses the face
	ipoints <- do.call(rbind,lapply(names(foundList),
		function(x){
			y<-foundList[[x]];
			if(nrow(y)>0){rownames(y)<-paste(x,seq_len(nrow(y)),sep=";")};y})
		)
	# we want to order them along the line of the chord
	npoints <- rbind(ipoints,chord.from.xy)

	bottom <- npoints[npoints[,2]==min(npoints[,2]),,drop=FALSE]
	bottomleft <- bottom[bottom[,1]==min(bottom[,1]),,drop=FALSE]

	dist <- (npoints[,1]-bottomleft[1])^2+(npoints[,2]-bottomleft[2])^2
	npoints <- npoints[order(dist),,drop=FALSE]
	npoints

}



.addIntersectingFace <- function(new1,new2,tempface2Name,face2IntersectionPoints) {
	# for each intersection point there is a (set of) edges of face2 to the next intersection 
	# point that we need to add as a (multiple) edge
	# 
	faceEdgeList <- .SplitFaceAtintersections(new2,tempface2Name,face2IntersectionPoints)
	set2Name  <- names(new2@setList)[1]
	
	# now each element of faceEdgeList is a set of edges from one intersection point to another
	# if we had edges duplicated between drawing1 and 2, dont need to add them in ( I think)
	
	currentFaceEdges <- do.call(c,lapply( .faceNames(new1), .faceEdgeNames,drawing=new1))
	seenEdges <- lapply(faceEdgeList,function(flist){flist %in% currentFaceEdges})
	lapply(seenEdges,function(seen){if (length(seen)>1 & !all(seen==seen[1])){stop("Some edges seen others not")}})
	seenEdges <- sapply(seenEdges,all)	
	faceEdgeList <- faceEdgeList[!seenEdges]

	for (edgeSet in faceEdgeList) { 
#cat("Adding edge", edgeSet,"\n")
		new1 <- injectEdge(drawing=new1,newEdgeList=new2@edgeList[edgeSet],set2Name=set2Name,addToList=FALSE)
#show(new1)
	}
	# that will have split and correctly renamed all the faces of drawing1 that the edges passed through
	# we need to catch any others
	faceNames <- .faceNames(new1,onlyVisible=TRUE) 
	faceSignatures <- .faceSignatures(new1,onlyVisible=TRUE) 
	nSets <- length(new1@setList)
	unchangedFaces <- faceNames[!sapply(faceSignatures ,function(x)nchar(x)==length(new1@setList))]
	for (faceName in unchangedFaces) {
		if (.is.face.within.set(drawing=new1,faceName=faceName,setName=set2Name)) {
			suffix <- "1"
		} else {
			suffix <- "0"
		}
		newFaceName<- paste(faceName,suffix,sep="")
		new1 <- renameFaces(new1,faceName,newFaceName)
		new1 <- updateSignature(new1,newFaceName,suffix)
	}
	new1
}

.SplitFaceAtintersections <- function(new2,tempface2Name,face2IntersectionPoints) {
	face2Points <- .points.of.face(drawing=new2,tempface2Name)
	ixpoints <- sapply(face2Points,function(x){x%in%face2IntersectionPoints})

	# breaks the face into a list of edge-lists, each starting and finishing at an intersection point as defined by logical vector ixpoints
	ix <- min(which(ixpoints))
	# rearrange edges so the first one starts at an intersection point
	new2 <- .startFaceAtPoint(drawing=new2,faceName=tempface2Name,from=.points.of.face(new2,tempface2Name)[ix])
	# then recalculate
	face2Points <- .points.of.face(drawing=new2,tempface2Name)
	# now first member of point list is an intersection point
	ixstart <- which(sapply(face2Points,function(x){x%in%face2IntersectionPoints}))
	ixend <- (ixstart-1)
	ixend <- c(ixend[-1],length(face2Points))
	faceEdgeNames <- .faceEdgeNames(drawing=new2,tempface2Name)
	faceEdgeList <- list()
	for (ix in 1:length(ixstart)) {
		faceEdgeList[[ix]] <- faceEdgeNames[ seq(ixstart[ix],ixend[ix]) ] 
	}
	faceEdgeList
}




.internal.edge.drawing.intersection <- function(drawing,edge) {
	for (edgeName in names(drawing@edgeList)) {
		found <- .findIntersection(edge1=drawing@edgeList[[edgeName]],edge2=edge)
		# all the intersections, want to exclude from and to points
		fxy <- drawing@nodeList[[edge@from]]; txy <- drawing@nodeList[[edge@to]]
		if(nrow(found)>0) {
			isf <- sapply(seq_len(nrow(found)),function(ix){isTRUE(all.equal(as.numeric(found[ix,,drop=FALSE]),as.numeric(fxy)))})
			ist <- sapply(seq_len(nrow(found)),function(ix){isTRUE(all.equal(as.numeric(found[ix,,drop=FALSE]),as.numeric(txy)))})		
			found <- found[!isf & !ist,,drop=FALSE]
			}
		if (nrow(found)>0) {
			return(TRUE)
		}
	}
	return(FALSE)
}

.find.point.in.nodelist <- function(drawing,point.xy) {
	existing.xy <- do.call(rbind,drawing@nodeList)
	dist <- existing.xy;
	dist[,1] <- dist[,1]-point.xy[1]; dist[,2] <- dist[,2]-point.xy[2]
	d2 <- apply(	dist^2,1,sum)
	iszero <- which(d2 <  (.Machine$double.eps ^ 0.5))
	if (length(iszero)==0) { return(NA) } 
	return(names(drawing@nodeList)[iszero])
}
	
.find.duplicate.point <- function(drawing1,drawing2,pointName) {
	.find.point.in.nodelist(drawing1,drawing2@nodeList[[pointName]])
}

.node.number.unused <- function(drawing) {
	max(as.numeric(gsub("[^0-9]","",c(names(drawing@nodeList)))))+1
}

.add.intersection.points <- function(drawing1,drawing2) {
	new1 <- drawing1; new2 <- drawing2;
	intersectionPoints <- character(0);
	# max number used in nodenames to avoid clashes
	nix <- max(.node.number.unused(drawing1),.node.number.unused(drawing2))
	fres <- NULL
	for (edgeName1 in names(drawing1@edgeList)) {
		for (edgeName2 in names(new2@edgeList)) {
#{cat(sprintf("looking for intersections between %s and %s\n",edgeName1,edgeName2));} 
			found <- .findIntersection(edge1=drawing1@edgeList[[edgeName1]],edge2=new2@edgeList[[edgeName2 ]])
#if(nrow(found)>0){cat(sprintf("Intersections between %s and %s\n",edgeName1,edgeName2));show(found)} 
			# check if any of those were intersections we already had names for
			if(nrow(found)>0) {
				n1nodes <- apply(found,1,.find.point.in.nodelist,drawing=new1)
				n2nodes <- apply(found,1,.find.point.in.nodelist,drawing=new2)
				rownames(found) <-  paste("i",nix+seq_len(nrow(found)),sep="")
				nix <- nix + nrow(found)
				colnames(found)<- NULL
				rownames(found)[!is.na(n2nodes)] <- n2nodes[!is.na(n2nodes)]
				rownames(found)[!is.na(n1nodes)] <- n1nodes[!is.na(n1nodes)]
				new2nodes <- unique(rownames(found)[is.na(n2nodes)])
				new1nodes <- unique(rownames(found)[is.na(n1nodes)])
				found <- found[!duplicated(rownames(found)),,drop=FALSE]
#show(found)
#show(n1nodes)
#show(n2nodes)
				found.df <- data.frame(found); colnames(found.df) <- c("x","y")
				found.df$nodeName <- rownames(found)
				found.df$edgeName1 <- edgeName1
				found.df$edgeName2 <- edgeName2
				if (is.null(fres)) { fres <- found.df } else {
					fres <- rbind(fres,found.df)
				}
			}
		}
	} #edgeName1
	# we have to store them up and then add them all at once so the edge names dont change as we are finding the
	if (is.null(fres)) { 
		res <- list(new1=new1,new2=new2,intersectionPoints=character(0));
	} else {
#browser()
		rownames(fres) <- 1:nrow(fres)
		fres <- fres[!duplicated(fres$nodeName),]
		for (edgeName1 in unique(fres$edgeName1)) {
			found1 <- fres[fres$edgeName1==edgeName1,]
			found1 <- subset(found1, !found1$nodeName %in% names(new1@nodeList))
			newPoints <- data.matrix(found1 [,c("x","y")]); rownames(newPoints ) <- found1$nodeName
			new1 <- injectPoints(drawing=new1 ,edgeName=edgeName1,newPoints=unique(newPoints)) 
		}
		for (edgeName2 in unique(fres$edgeName2)) {
			found2 <- fres[fres$edgeName2==edgeName2,]
			found2 <- subset(found2, !found2$nodeName %in% names(new2@nodeList))
			newPoints <- data.matrix(found2 [,c("x","y")]); rownames(newPoints ) <- found2$nodeName
			new2 <- injectPoints(drawing=new2 ,edgeName=edgeName2,newPoints=unique(newPoints) ) 
		}
		intersectionPoints <- unique(fres$nodeName)
		res <- list(new1=new1,new2=new2,intersectionPoints=intersectionPoints);
	}
	res
}
#debug(.add.intersection.points)

.node.distance <- function(xy1,xy2) {
	distmat <- matrix(NA,nrow=nrow(xy1),ncol=nrow(xy2))
	for (i in seq_len(nrow(xy1))) {
		for (j in seq_len(nrow(xy2))) {
			distmat[i,j] <- sqrt(sum((xy1[i,] - xy2[j,])^2))
		}
	}
	distmat
}
.nodes.identical <- function(xy1,xy2) {
	stopifnot(nrow(xy2)==1)
	distmat <- .node.distance(xy1,xy2)
	idmat <- abs(distmat) < (.Machine$double.eps ^ 0.5)
	wix <- which(idmat)
	if(length(wix)==0) { return(NA) }
	return(min(wix))
}


.check.clashes <- function(drawing1,drawing2) {
	if (any(names(drawing2@setList) %in% names(drawing1@setList))) {
		stop("Clashing set names")
	}
	# now we add all the nodes from new2 in
	new2 <- drawing2
	# do we have nodes with the same name but at different places?
	newNodes <- names(new2@nodeList)
	oldNodes <- newNodes[newNodes %in% names(drawing1@nodeList)]
	newNodes <- newNodes[!newNodes %in% names(drawing1@nodeList)]
	for (node in oldNodes) {
		p.1 <- drawing1@nodeList[[node]];
		p.2 <- new2@nodeList[[node]]
		nix <- .nodes.identical(p.1,p.2);
		if (is.na(nix)) {
			stop(sprintf("Node %s is different in two drawings",node))
		}
	}
	# or nodes which have different names but in the same place?
	for (node in newNodes) {
		p.1 <- do.call(rbind,drawing1@nodeList);
		p.2 <- new2@nodeList[[node]]	
		nix <- .nodes.identical(p.1,p.2);
		if (!(is.na(nix))) {
#			cat(sprintf("Node %s is the same as %s;\n",names(drawing1@nodeList)[nix],node))
			new2 <- rename.node(drawing=new2,oldName=node,newName=names(drawing1@nodeList)[nix])
		}
	}
	new2
}

rename.node <- function(drawing,oldName,newName) {
	names(drawing@nodeList)[names(drawing@nodeList)==oldName] <- newName
	for (ix in seq_len(length(drawing@edgeList))) {
		if (drawing@edgeList[[ix]]@from == oldName) {
			drawing@edgeList[[ix]]@from <- newName
		}
		if (drawing@edgeList[[ix]]@to == oldName) {
			drawing@edgeList[[ix]]@to <- newName
		}
	}
	drawing
}

.check.duplicate.edges <- function(drawing1,drawing2) {
	# first we check the forward edges....

	fromtomat1 <- sapply(drawing1@edgeList,function(y){paste(y@from,y@to)})
	fromtomat2 <- sapply(drawing2@edgeList,function(y){paste(y@from,y@to)})
	commonfromto <- fromtomat2 %in% fromtomat1
	common2 <- fromtomat2[commonfromto]
	common1 <- lapply(common2,function(x)fromtomat1[fromtomat1==x])
	for (edgeName2 in names(common1)) {
		edgeNames1 <- names(common1[[edgeName2]])
		edges1 <- drawing1@edgeList[edgeNames1]
		edge2 <- drawing2@edgeList[[edgeName2]]
		for (edgeName1 in names(edges1)) {
			edge1 <- edges1[[edgeName1]]
			if (.identical(edge1,edge2)) {
#cat(sprintf("Replacing edge %s by %s\n",edgeName2,edgeName1))
				drawing2 <- .rename.edge(drawing2,edgeName2,edgeName1)
				if (edgeName2 %in% names(drawing1@edgeList)) {
					drawing1@edgeList[[edgeName2]] <- NULL
				}
#show(names(drawing1@edgeList))
#cat("==----\n")
#show(names(drawing2@edgeList))
#cat("======\n")
			}
		}
	}
	# then the reverse 
	fromtomat1 <- sapply(drawing1@edgeList,function(y){paste(y@from,y@to)})
	fromtomat2 <- sapply(drawing2@edgeList,function(y){paste(y@to,y@from)}); # nb reverse fromto
	commonfromto <- fromtomat2 %in% fromtomat1
	common2 <- fromtomat2[commonfromto]
	common1 <- lapply(common2,function(x)fromtomat1[fromtomat1==x])
	for (edgeName2 in names(common1)) {
		edgeNames1 <- names(common1[[edgeName2]])
		edges1 <- drawing1@edgeList[edgeNames1]
		edge2 <- .reverseEdge(drawing2@edgeList[[edgeName2]]) # nb reverse
		for (edgeName1 in names(edges1)) {
			edge1 <- edges1[[edgeName1]]
			if (.identical(edge1,edge2)) {
#cat(sprintf("Replacing edge %s by %s\n",edgeName2,edgeName1))
				revName1<- sub("^--","",paste("-",edgeName1,sep=""))
				drawing2 <- .rename.edge(drawing2,edgeName2,revName1)
				if (edgeName2 %in% names(drawing1@edgeList)) {
					drawing1@edgeList[[edgeName2]] <- NULL
				}
#show(names(drawing1@edgeList))
#cat("==----\n")
#show(names(drawing2@edgeList))
#cat("======\n")
			}
		}
	}


	list(new1=drawing1,new2=drawing2)
}

.rename.edge <- function(drawing,oldEdgeName,newEdgeName) {
	oldEdgeName <- sub("^-","",oldEdgeName)
	revEdgeName <- paste("-",oldEdgeName,sep="")
	revNewEdgeName <- sub("^--","",paste("-",newEdgeName,sep=""))
	names(drawing@edgeList)[names(drawing@edgeList)==oldEdgeName] <- newEdgeName
	names(drawing@edgeList)[names(drawing@edgeList)==revEdgeName ] <- revNewEdgeName 

	faceNames <- .faceNames(drawing)
	for (faceName in faceNames) {
		drawing <- spliceEdgeIntoFace (drawing,faceName,oldEdgeName,newEdgeName,doReverse=TRUE) 
	}

	drawing@setList <- lapply(drawing@setList ,function(flist) {
		flist[flist==oldEdgeName ] <- newEdgeName
		flist[flist==revEdgeName ] <- revNewEdgeName 
		flist
	})
	drawing
}

remove.nonintersectionpoints <- function(drawing) {
	# a non intersection point is a named point (usually from when the face was a single edge)
	# at which there are no intersections with any other sets/invisible edges
	# ie it is at the end of exactly one edge
	#browser()
	# rely on having only the edges in the drawing in the edgelist, and only have one of each orientation
	toNodes <- lapply(drawing@edgeList,function(x)x@to)
	fromNodes <- lapply(drawing@edgeList,function(x)x@from)
	toNode.count <- table(as.character(toNodes))
	fromNode.count <- table(as.character(fromNodes))
	noni.nodes <- intersect(names(toNode.count)[toNode.count==1],names(fromNode.count)[fromNode.count==1])
	for (node in noni.nodes) {
		inedgeName <- names(toNodes)[toNodes==node]
		outedgeName<- names(fromNodes)[fromNodes==node]
		drawing <- joinEdgesInDrawing(drawing,inedgeName ,outedgeName)
		}
	drawing
}

getEdge <- function(drawing,edgeName) {
	edgeUnsigned <- sub("^-","",edgeName)
	edge <- drawing@edgeList[[edgeUnsigned]]
	if (edgeUnsigned != edgeName) {
		edge <- .reverseEdge(edge)
	}
	edge
}

joinEdgesInDrawing <- function(drawing,inedgeName ,outedgeName) {
	inrev <- substr(inedgeName,1,1)=="-"
	outrev <- substr(outedgeName,1,1)=="-"
	if (inrev !=outrev) {stop("Cant cope with joining edges of opposite polarity")}
	if (!inrev){
		inpos <- inedgeName; inneg <- paste("-",inedgeName,sep="")
		outpos <- outedgeName; outneg <- paste("-",outedgeName,sep="")
	}	else { 
		inpos <- sub("^-","",outedgeName);inneg<- outedgeName
		outpos <- sub("^-","",inedgeName); outneg <- inedgeName
	}
	inedge <- getEdge (drawing,inpos) 
	outedge <-getEdge (drawing,outpos) 
	if (inedge@to != outedge@from) { stop(sprintf("Edges do not joint at single point: %s,%s\n",inedge@to,outedge@from))}
	if (class(inedge) != class(outedge)) { stop(sprintf("Cant join edges of different classes")) }
	#if (class(inedge) != "VDedgeLines") {stop(sprintf("Cant join edges of non edgeLine classes")) }


	newEdge <- joinEdges(inedge,outedge)

	inEdgeSplit <- strsplit(inedgeName ,split="|",fixed=TRUE)[[1]];
	if (length(inEdgeSplit)>=3) {
		Set <- inEdgeSplit[3]
	} else {
		Set <- "X"
	}
#cat(sprintf("Joining %s and %s in Set %s\n",inedgeName,outedgeName,Set))
	newEdgeName <- paste(inedge@from,outedge@to,gsub("[A-Za-z]","",Set),sep="|")
	
	drawing@edgeList[[ inpos]] <- NULL
	drawing@edgeList[[ outpos]] <- NULL
	drawing@edgeList[[ newEdgeName]]  <- newEdge
	for (setName in names(drawing@setList)) {
		fedges <- .faceEdgeNames(drawing,setName,type="set")
		if (fedges[1]==outpos) { fedges <- fedges[c(2:length(fedges),1)] }
		fedges <- spliceinstead(fedges,c(inpos,outpos),newEdgeName)
		drawing@setList[[setName]] <- fedges
	}
	for (fname in .faceNames(drawing)) {
		fedges <- .faceEdgeNames(drawing,fname,type="face")
		# first look for +ve pairs
		if (fedges[1]==outpos) { fedges <- fedges[c(2:length(fedges),1)] }
		fedges <- spliceinstead(fedges,c(inpos,outpos),newEdgeName)
		# then negative ones
		if (fedges[1]==inneg){ fedges <- fedges[c(2:length(fedges),1)] }
		revName <- paste("-",newEdgeName,sep="")
		fedges <- spliceinstead(fedges,c(outneg,inneg),revName)
		drawing@faceList[[fname]] <- fedges
	}
	drawing@nodeList[[inedge@to]] <- NULL
	drawing	
}

.merge.faces.invisibly.split <- function(diagram) {
	doneamerge<- TRUE
	while (doneamerge) {
		res <- .try.merge.faces.invisibly.split(diagram)
		doneamerge <- res$merged
		diagram <- res$diagram
	}
	diagram
}

.try.merge.faces.invisibly.split <- function(diagram) {
	# first we identify multople faces with the same signature
	fsigs <- data.frame(cbind(Name=unlist(.faceNames(diagram)),Signature=unlist(.faceSignatures(diagram))),stringsAsFactors=FALSE);
	rownames(fsigs)<- 1:nrow(fsigs)
	nSets <- unique(nchar(setdiff(fsigs$Signature,"DarkMatter")))
	fsigs$Signature[fsigs$Signature==paste(rep("0",nSets),collapse="")] <- "DarkMatter"
	FacesPerSignature <- lapply(split(fsigs$Name,fsigs$Signature),length)
	FacesPerSignature <- FacesPerSignature [FacesPerSignature >1]
	doingamerge <- FALSE
	if (length(FacesPerSignature )>0) {
		for (wn in names(FacesPerSignature )) {
			wnames <- fsigs$Name[fsigs$Signature==wn]
			if (length(wnames)!=2) { warning("Can't merge multiple invisibly split faces, just trying first two")}
			faceEdges1 <- .face.to.faceEdges(diagram,wnames[1])
			faceVisible1 <- sapply(faceEdges1,function(x)x@visible)	
			if (all(faceVisible1)) { break; }		
			faceEdges2 <- .face.to.faceEdges(diagram,wnames[2])
			faceVisible2 <- sapply(faceEdges2,function(x)x@visible)	
			if (all(faceVisible2)) { break; }		
			commonInvisibleEdges <- intersect(sub("^-","",names(faceEdges1)),sub("^-","",names(faceEdges2)))
			if (length(commonInvisibleEdges)==0) break;
			doingamerge <- TRUE
			firstVisibleIx <- min(which(faceVisible1))
			firstVisibleFrom <- faceEdges1[[firstVisibleIx ]]@from
			diagram<- .startFaceAtPoint(diagram,wnames[1],firstVisibleFrom )
			faceEdges1 <- .face.to.faceEdges(diagram,wnames[1])
			faceVisible1 <- sapply(faceEdges1,function(x)x@visible)	
			faceEdges2 <- .face.to.faceEdges(diagram,wnames[2])
			faceVisible2 <- sapply(faceEdges2,function(x)x@visible)	
			Invisibles1 <-names( faceEdges1)[!faceVisible1]
			Invisibles2 <- rev(sub("^--","",paste("-",Invisibles1,sep="")))
			faceEdgeNames1 <- .faceEdgeNames(diagram,wnames[1])
			faceEdgeNames2 <- .faceEdgeNames(diagram,wnames[2])
			firstInvisibleIx1 <-  min(which(!faceVisible1))
			lastInvisibleIx1 <-  max(which(!faceVisible1))
			beforeNames1 <- faceEdgeNames1[1:(firstInvisibleIx1-1)]
			afterNames1 <- if(lastInvisibleIx1 <length(faceEdgeNames1)) { faceEdgeNames1[(lastInvisibleIx1+1) :length(faceEdgeNames1)] } else {character(0)}
			firstInvisibleIx2 <-  min(which(!faceVisible2))
			lastInvisibleIx2 <-  max(which(!faceVisible2))
			beforeNames2 <- if(firstInvisibleIx2 >1)faceEdgeNames2[1:(firstInvisibleIx2 -1)]else { character(0)}
			afterNames2 <- if(lastInvisibleIx2 <length(faceEdgeNames2)) { faceEdgeNames2[(lastInvisibleIx2+1) :length(faceEdgeNames2)]	} else {character(0)}

			newEdgeNames <- c(beforeNames1,afterNames2,beforeNames2,afterNames1)
			diagram@faceList[[wnames[1]]] <- newEdgeNames
			diagram@faceList[[wnames[2]]] <- NULL
			diagram@faceSignature[[wnames[2]]] <- NULL
			lostedges <- unique(sub("^-","",c(Invisibles1 ,Invisibles2)))
			diagram@edgeList[lostedges] <- NULL
		}
	}
	return(list(diagram=diagram,merged=doingamerge))
}
	
