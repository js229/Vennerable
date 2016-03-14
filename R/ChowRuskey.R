#warning("Entering ChowRuskey")

TDtograph <- function(TD) {
	nodes <- names(TD@nodeList)
	# need to delete nonintersection points recognised by two successive 
	# edges in a face corresponding to the same set
	gr <- new("graphNEL",nodes=nodes,edgemode="directed")
	nodeDataDefaults(gr,"x") <- NA; nodeDataDefaults(gr,"y") <- NA; 

	for (node in nodes) {
		nodeData(gr,node,"x") <- TD@nodeList[[node]][,1]
		nodeData(gr,node,"y") <- TD@nodeList[[node]][,2]
	}
	edgeDataDefaults(gr,"Set") <- NA
	edgeDataDefaults(gr,"Name") <- NA
	edgeDataDefaults(gr,"Face+") <- NA
	edgeDataDefaults(gr,"Face-") <- NA
	for (sname in names(TD@setList)) {
		sedges <- .face.to.faceEdges(drawing=TD,faceName=sname,type="set")
		for (edgeName in names(sedges)) {
			edge <- sedges[[edgeName]]	
			gr <- addEdge(from=edge@from,to=edge@to,graph=gr)
			edgeData(gr,from=edge@from,to=edge@to,"Set") <- which(sname==names(TD@setList))
			edgeData(gr,from=edge@from,to=edge@to,"Name") <- edgeName
		}
	}
	for (fname in .faceNames(TD)) {
		fedgeNames <-  .faceEdgeNames(drawing=TD,faceName=fname,unsigned=FALSE)
		fedges <- .face.to.faceEdges(drawing=TD,faceName=fname)
		for (edgeName in fedgeNames ) {
			edge <- fedges[[edgeName ]]
			faceReverse <-  substr(edgeName,1,1)=="-"
			if (!faceReverse) {
				edgeData(gr,from=edge@from,to=edge@to,"Face+") <- fname
			} else {
				edgeData(gr,from=edge@to,to=edge@from,"Face-") <- fname
			}
		}
	}

	gregions <- list()
	for (fname in .faceNames(TD)) {
		fedges <- .face.to.faceEdges(drawing=TD,faceName=fname)
		fnodes <- sapply(fedges,function(x)x@to)
		gregion <- subGraph(gr,snodes=fnodes)
		gregions[[fname]] <- gregion
	}


	list(graph=gr,regions=gregions)
}

scythegr <- function(TDgr) {
	faces <- unlist(edgeData(TDgr,attr="Face+"))
	vfaces <- setdiff(faces,"DarkMatter")
	nSets <- unique(nchar(vfaces)); stopifnot(length(nSets)==1)
	dualPath <- sapply(0:(nSets),function(i){paste(paste(rep("1",nSets-i),collapse=""),paste(rep("0",i),collapse=""),sep="")})
	dualPath[length(dualPath)] <- "DarkMatter"
	if (nSets==4) {dualPath[2] <- "1101"} # to agree with Chow Ruskey example
	ed <- edgeData(TDgr)
	faces1 <- sapply(ed,function(x)x$"Face+")
	faces2 <- sapply(ed,function(x)x$"Face-")
	for (ix in 1:(length(dualPath)-1)) {
		face1 <- dualPath[ix]; face2 <- dualPath[ix+1]
		edge1 <- faces1[faces1 == face1]
		edge2 <- faces2[faces2 == face2]
		cutedge <- intersect(names(edge1),names(edge2))	; stopifnot(length(cutedge)==1)
		fromto <- strsplit(cutedge,split="|",fixed=TRUE)[[1]]
		TDgr <- removeEdge(graph=TDgr,from=fromto[1],to=fromto[2])
	}
	TDgr
}


node.to.ray <- function(pnodes,atsort) {
	by.tsort <- sapply(pnodes,function(pnode){
		which(sapply(atsort,function(nodes)pnode%in%nodes))}
	)
	names(by.tsort) <- pnodes
	2* (by.tsort)-1
}

.remove.loose.nodes <- function(G) {
	loose.nodes <- nodes(G)[( sapply(edges(G),length)==0 ) & (sapply(inEdges(G),length)==0)]	
	if (length(loose.nodes)>0) {
		G <- removeNode(loose.nodes,G)
	}
	G
}

.cycle.to.path <- function(face) {
	node.path <- nodes(face)[1]
	repeat {
		next.edge <- edges(face,node.path[length(node.path)])[[1]]
		node.path <- c(node.path,next.edge)
		if (next.edge==nodes(face)[1]) { break }
	}
	node.path
}

make.setlist <- function(region,atsort) {
	# we know for each region the fixed and free  edges around it
	# we convert the node names to ray numbers
	# and use the Set attribute to assign ray points to each set
	# on the free boundary
	rays <- node.to.ray(nodes(region),atsort) #the mapping from the topo sort 

	twok <- 2*length(atsort)

	free.region <- region
	fixed.region <- region
	fixed.edges <- unlist(edgeData(region,attr="fixed"))
	for (edgebar in names(fixed.edges)) {
		fromto <- strsplit(edgebar,split="|",fixed=TRUE)[[1]]
		if (fixed.edges[edgebar]) {
			free.region <- removeEdge(fromto[1],fromto[2],free.region)
		} else {
			fixed.region <- removeEdge(fromto[1],fromto[2],fixed.region)
		}
	}
	free.region <- .remove.loose.nodes(free.region) 
	fixed.region <-.remove.loose.nodes(fixed.region) 

	isStart <- (length( fixed.edges[fixed.edges])==0)
	if (isStart) { 
		free.node.path <- .cycle.to.path(free.region)
	} else  {
		free.node.path <- tsort(free.region)
		fixed.node.path <- tsort(fixed.region) 
	}

	EdgeLabels <- lapply(1: (length(free.node.path)-1), function(n) {
		from <- free.node.path[n]
		to <- free.node.path[n+1]
		rayfrom <- rays[from]
		rayto <- rays[to]
		if (rayto<rayfrom) {
			rayrange <- c( (rayfrom:(twok)),1:rayto) 
		} else {
			rayrange <- rayfrom:rayto
		}
		Set <-edgeData(free.region,from=from,to=to,attr="Set")[[1]]
		rayptnames <- rep(NA,length(rayrange))
		list(SetNumber=Set,Raypoints=rayrange,RayPointNames=rays[c(to,from)])	
	})

	if (!isStart) {
		EdgeLabels[[1]]$Raypoints <- EdgeLabels[[1]]$Raypoints[-1] # the first point is on the already drawn set
		EL <- EdgeLabels[[length(EdgeLabels)]]
		EL$Raypoints <- EL$Raypoints[-length(EL$Raypoints) ] # and the last .
		EdgeLabels[[length(EdgeLabels)]] <- EL
	}

	EdgeLabels
}

make.setlist.from.AWFE <- function(G,regions,regionOrder,atsort) {
	setlist <- list()
	edgeDataDefaults(G,"fixed") <- FALSE

	for (vs in regionOrder[-length(regionOrder)]){ # last entry is "0000..."
#cat(vs,"\n");
		region <- regions[[vs]]
		edgeDataDefaults(region,"fixed") <- FALSE
		for (from in names(edges(region))) {
			for (to in edges(region)[[from]]) {
				edgeData(region,from,to,"fixed") <- edgeData(G,from,to,"fixed")
			}
		}
		setlist[[vs]] <- make.setlist(region,atsort)
		
		vs.edges <- edges(region)
		for (from in names(vs.edges)) {
			for (to in edges(region)[[from]]) {
				edgeData(G,from,to,"fixed") <- TRUE
			}
		}
	}
	setlist
}



.unify.rays <- function(asetlines,twok) {
	# given the start and end points of each Set in the rayPoints
	# our task is to join them up and supply the outer two rays
	# (unless we already have the whole circle)
	
	raystarts <- sapply(asetlines,function(x)x$Raypoints[1])
	rayends <- sapply(asetlines,function(x)x$Raypoints[length(x$Raypoints)])
	thru <- rayends < raystarts
	raystartsabs <- ifelse(thru,raystarts-twok,raystarts)
	ors <- raystarts[order(raystartsabs)]
	ore <- rayends[order(raystartsabs)]
	# should all be overlapping 
	stopifnot( length(ors==1) | (ors[2:length(ors)] == ore[1:(length(ore)-1)]))
	rs <- ors[1]
	re <- ore[length(ore)]
	mod2k1 <- function(n) { 1+ (n-1)%% twok }
	if (rs==re) {
		allrays <- 1:twok
	} else {
		rs <- mod2k1(rs-1); 
		re <- mod2k1(re+1);
	
		if (rs>=re) {
			allrays <- c(rs:twok,1:re)
		} else {
			allrays <- (rs:re)
		}
	} 
	allrays
}


makesrp <- function(vs,asetlines,SetRayPoints,Weight,outerRay,angleray) {
	# helper function for compute.CR
	# given an intersection set (vs), with a topology specified in asetlines,
	# first of all compute the (outer) free edges of the intersection set
# { cat(sprintf("%s \n",vs))}
	
	# the rays on which the set is drawn
	raypts <- .unify.rays(asetlines,ncol(SetRayPoints))

	# the offset from the fixed (inner) edges to the free ones
	delta <- deltagivenouter (outerRay[raypts],Weight[vs],angleray)
	# then assign the correct labels to the outer edges
	deltaix <- 1
	for (EL in asetlines) {
		rp <- EL$Raypoints
		# special case for the central point
		if (length(raypts)==ncol(SetRayPoints)) {
			addto <- unique(delta) ; 
			stopifnot(length(addto)==1)
		} else {
			addto <- delta[ seq(deltaix,deltaix+length(rp)-1)]
		}
		SetRayPoints[EL$SetNumber,rp] <- 	outerRay[rp]+addto
		deltaix <- deltaix+length(rp)-1
	}
	SetRayPoints
}

makeirs <- function(vs,asetlines,SetRayPoints,IntersectionRaySets,outerRay) {
	# once the SetRayPoints for the current intersection set have been computed,
	# use them, plus the topology, to compute the perimeter of the intersection set by
	# taking the right bits from the right set

	# the intersection polygon goes out along the inside
	raypts <- .unify.rays(asetlines,ncol(SetRayPoints))
	isStart = (length(raypts)==ncol(SetRayPoints))
	if (!isStart) {
		IntersectionRaySets[[vs]] <- cbind(raypts,outerRay[raypts])
		# then back (hence the rev) along the outside
		for (setline in rev(asetlines)) {
			IntersectionRaySets[[vs]] <- rbind(IntersectionRaySets[[vs]],
				cbind(rev(setline$Raypoints),SetRayPoints[setline$SetNumber,rev(setline$Raypoints)])
			)
		}
	} else {
		IntersectionRaySets[[vs]] <- cbind(numeric(0),numeric(0)) # inner polygon has no inner
		for (setline in (asetlines)) {
			IntersectionRaySets[[vs]] <- rbind(IntersectionRaySets[[vs]],
				cbind((setline$Raypoints),SetRayPoints[setline$SetNumber,(setline$Raypoints)])
			)
		}
	}
	IntersectionRaySets
}

make.maxiray <- function(irs,twok) {
	irs.start <- irs[1,1]
	mod2k1 <- function(n) { 1+ (n-1)%% twok }
	irs[,1] <- mod2k1(irs[,1]-irs[1,1])
	irsrange <- lapply(split(irs[,2],irs[,1]),function(x)(range(x)))
	irsdiff <- sapply(irsrange,function(x)diff((x)))
	if (length(irsdiff)==twok) { # the inner ring
		midray <- 1
		midpoint <- 0
	} else  {
		maxirays <- names(irsdiff)[irsdiff==max(irsdiff)]
		midray <- as.vector(quantile(as.numeric(maxirays),0.5,type=1))
		midpoint <- mean(irsrange[[as.character(midray)]])
		midray <- mod2k1(midray+irs.start)
	}
	res <- matrix(c(midray,midpoint),ncol=2)
	return(res)
}	

.dual.region.order <- function(regionNames ) {
	nnsplit <- strsplit(regionNames ,split="")
	nncount <- sapply(nnsplit,function(x)length(x[x=="1"]))
	regionNames [order(-nncount)]
}

compute.CR <- function(V,doWeights=TRUE) {
	nSets <- NumberOfSets(V)
	Vorig <- V
	if (!doWeights) {	Weights(V) <- 1 + 0 * Weights(V)}
	
	AWFE.diagram <- compute.AWFE(V,type="battle")
	AWFE.graphregions <- TDtograph(AWFE.diagram)
	AWFE.graph <- AWFE.graphregions$graph
	AWFE.regions <- AWFE.graphregions$regions

	sTDgr <- scythegr(AWFE.graph)
	ray.to.point.data <- my.tsort(sTDgr)
	k <- 	length(ray.to.point.data)
	twok <- 2*k

	regionNames <-	.faceNames(AWFE.diagram)
	regionNames[regionNames=="DarkMatter"] <- dark.matter.signature(V)
	regionOrder <- .dual.region.order(regionNames )
	
	####################
	# then prepare the graph we want to construct
	#################
	Weight <- Weights(V)
		
	
	angleray <- 2*pi / twok 
	outerRay <- rep(0,twok )
	SetRayPoints <- matrix(NA,nrow=nSets,ncol=twok )

	# we work in ray points with coordinates measured out along each ray
	# at each stage, the farthest known point out on each ray is in outerRay
	# SetRayPoints, the point at which each Set is on the ray will be built up as we work through 
	# the regions in topological order
	# IntersectionRaySets, is the (Set-labelled) list of ray points specifying the edges of each face
	
	setlines <- make.setlist.from.AWFE (G=AWFE.graph,regions=AWFE.regions,regionOrder,atsort=ray.to.point.data)

	
	# this is the loop in topological order of G*
	for (vs in names(setlines)) { 
#cat(outerRay,"\n")
		SetRayPoints <-  makesrp(vs,setlines[[vs]],SetRayPoints,Weight,outerRay,angleray)		
		outerRay <- apply(SetRayPoints,2,max,na.rm=TRUE)
	}

	srp.df <- reshape::melt(SetRayPoints); colnames(srp.df) <- c("SetNumber","Ray","r")
	srp.df$Name <- NA
	# that's all the work done, just encode this as a tissue diagram now
	for (vs in names(setlines)) {
		for (ix in 1:length(setlines[[vs]])) {
			asetlines <- setlines[[vs]][[ix]]
			Set <- asetlines$SetNumber
			rpn <- asetlines$RayPointNames
			for (rp in rpn) {
				srp.df[srp.df$SetNumber==Set & srp.df$Ray ==rp,"Name"] <- names(rpn)[rpn==rp]
			}
		}
	}
	xdf <- .raypoint.to.xy(srp.df[,c("Ray","r")],angleray); colnames(xdf) <- c("x","y")
	srp.df <- cbind(srp.df,xdf)
	points.df <- unique(srp.df[!is.na(srp.df$Name),c("Ray","r","Name","x","y")])	


	# we already have the topology and edge patterns correctly set up in AWFE.diagram
	# all we have to do is adjust the position of its nodes, and redraw its edges
	TD <- AWFE.diagram
	TD@nodeList <- lapply(split(points.df[,c("x","y","Name")],points.df$Name),function(xydf){ w<- data.matrix(xydf[,c("x","y")]);rownames(w)<-xydf$Name;w})

	for (six in 1:length(names(TD@setList))) {
		setName <- names(TD@setList)[[six]]
		sedges <- .face.to.faceEdges(drawing=TD,setName,type="set")
		setray <- srp.df[srp.df$SetNumber==six,]
		setray <- rbind(setray,setray)
		for (sex in 1:length(sedges)) {
			ename <- names(sedges)[sex]
			edge <- sedges[[sex]]
			fromix <- which(setray$Name==edge@from)[1]
			tempray <- setray; tempray$Name[1:fromix]<- NA
			toix <- which(tempray$Name==edge@to)[1]
			xymat <- data.matrix(setray[fromix:toix,c("x","y")])
			edge@xy <- xymat 
			edge@bb <- rbind(apply(xymat,2,min),apply(xymat,2,max))
			# this assumes the previous AWFE edge was a line, which it wouldnt be if we used classic rather than battle AWFE
			TD@edgeList[[ename]] <- edge
		}
	}
	# also delete the empty faces
	emptyFaces <- names(Weight)[Weight==0]
	for (fname in emptyFaces) {
		TD@faceList[[fname]] <- NULL
		TD@faceSignature[[fname]] <- NULL
	}

	TD <- as(TD,"TissueDrawing")
	VD <- new("VennDrawing",TD,Vorig)
	SetLabels <- .default.SetLabelPositions(VD)
	VD <- VennSetSetLabels(VD,SetLabels)
	VD <- .square.universe(VD ,doWeights=doWeights)
	FaceLabels <- .default.FaceLabelPositions(VD)
	VD <- VennSetFaceLabels(VD,FaceLabels)
	VD
}

.raypoint.to.xy <- function(raypoints,angleray) {
	angles <- -angleray * (raypoints[,1])
	xy <- matrix( c(cos(angles), sin(angles)),ncol=2,byrow=FALSE)
	xy <- raypoints[,2] * xy
	xy
}	

.raypoint.area <- function(raypoints,angleray) {
	xy <- .raypoint.to.xy(raypoints,angleray)
	.polygon.area(xy)
}

compute.delta <- function( outervals, wght,angleray) {
	nrays <- length(outervals)
	#ca <- 2
	# try change
	ca <- (nrays-3) 
	singix <- c(1,2,nrays -1,nrays )
	cb <- sum(outervals[singix]) + 2* sum(outervals[-singix])
	cc <- -wght/(sin(angleray)/2)
	delta <- (-cb+sqrt(cb^2-4*ca*cc))/(2*ca)
	delta <- rep(delta,nrays -2)
	delta
}

deltagivenouter <- function( outervals, wght,angleray,smidge=1e-3) {
	# given n outervals as ray points on n successive rays separated by angleray
	# return n-2 deltas corresponding to the middle n-2 rays with area wght
	nray <- length(outervals); stopifnot( nray>2)
	if (wght==0) {
		delta <- rep(0,nray-2)
		return(delta)
	} 
	if (length(unique(outervals[2:(nray-1)]))==1) {
		# a nice uniform inside
		# special case: first time through all the outervals are zero and we have a polygon
		if (all(outervals==0)) {
			delta <- 	sqrt(wght/(length(outervals)*sin(angleray)/2))
			delta <- rep(delta,nray-2)
		} else {
			delta <- compute.delta ( outervals, wght,angleray) 
		}
		return(delta)
	} 
	# not uniform. Can we make it so?
	# where.max <- which(outervals==max(outervals))
	rmax <- max(outervals[-c(1,nray)])
	smidge <- min(smidge,wght/(rmax^2*nray*angleray) )
	rmax <- rmax*(1+smidge)

	feasible.shape.raypoints <- rbind(
		cbind(1:nray,outervals),
		cbind((nray-1):2,rmax))
	feasible.shape.area <- .raypoint.area(feasible.shape.raypoints,angleray )
	if (feasible.shape.area <= wght  ) {
		# that means we can fill in the feasible shape; create a uniform outside and call again to make the outer
		# shape (which will be again uniform ) 
		# nice and smooth; fill in the inner bit
			newouter <- c(outervals[1],rep(rmax,nray-2),outervals[nray])
			delta <- deltagivenouter(newouter, wght- feasible.shape.area,angleray)
			delta <- delta+rmax
			delta <- delta-outervals[2:(nray-1)]
	} else {
		if ( wght < feasible.shape.area  & feasible.shape.area < wght + 2*pi*rmax^2*smidge ) {
			# it's possible that the smidge is too large relative to the weight
			# to allow any feasible solution 
			# if this happens should choose a smaller smidge but we bodge it
			# the effect 
			warning("Try a different smidge; bodging")				
		}
		if (nray>5) {
			# if the maxes occur at the edges, take those out
			if (outervals[2]==max(outervals[-c(1,nray)]))	{
				delta.2 <- smidge
				smidge.triangle <- cbind(c(1,2,2),c(outervals[1],outervals[2]+smidge,outervals[2]))
				smidge.area <- .raypoint.area(smidge.triangle,angleray )
				sub.outer <- outervals[-1]
				sub.outer[1] <- sub.outer[1]+smidge
				sub.delta <- deltagivenouter(sub.outer, wght-smidge.area,angleray)
				delta <- c(delta.2,sub.delta)
				return(delta)
			} else if (outervals[nray-1]==max(outervals[-c(1,nray)])) {
				delta.2 <- smidge
				smidge.triangle <- cbind(c(nray,nray-1,nray-1),c(outervals[nray],outervals[nray-1],outervals[nray-1]+smidge))
				smidge.area <- .raypoint.area(smidge.triangle,angleray )
				sub.outer <- outervals[-nray]
				sub.outer[length(sub.outer)] <- sub.outer[length(sub.outer)]+smidge
				sub.delta <- deltagivenouter(sub.outer, wght-smidge.area,angleray)
				delta <- c(sub.delta,delta.2)
				return(delta)
			}
			# if the max is in the middle, want to divide the weights
			# but lets not bother for now
		} # nray>5

		# mostly though, will need to notch.
		# same delta for each ray (could do much better than this...)
		delta <- compute.delta ( outervals, wght,angleray) 
	}
	return(delta)
}


deltasmooth <- function( outervals, wght,angleray) {
	delta <- compute.delta (outervals, wght,angleray)
	if (length(outervals)!=5) {return(delta)}
	s <- (outervals[1]+outervals[2])/outervals[3]
	if (s>1) {
		delta <- delta * c(1/s,s,1/s)
	}
}

