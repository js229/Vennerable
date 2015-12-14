

my.tsort <- function(graph,type="wrong",res=list()) {
	if (type=="right") {
		incount <- sapply(inEdges(graph),length)
		top.nodes <- names(incount)[incount==0]
		if (length(top.nodes)==0) {
			return(list(nodes(graph)))
		}
		frompts <- acc(graph,top.nodes)
		acc.from <- do.call(rbind,lapply(names(frompts),function(x)data.frame(from=x,to=names(frompts[[x]]),acc=frompts[[x]])))
		acc.from <- subset(acc.from, ! acc.from$to %in% top.nodes)
		acc.dist <- sapply(split(acc.from$acc,acc.from$to),min)
		res <- c(list(top.nodes),split(names(acc.dist),acc.dist))
		names(res) <- NULL
		res
	} else {
		incount <- sapply(inEdges(graph),length)
		topts <- names(incount)[incount==0]
		res <- c(res,list(topts))
		graph <- removeNode(topts,graph)
		if (numNodes(graph)>0) {
			res <- my.tsort(graph,type,res)
		}
		res
	}
}


###########################
# general dual graph on 4 Sets
# ie the ncube Qn
makeQn <- function(nSets) {
	E4.Vs <- new("graphNEL",
		nodes=VennSignature(Venn(numberOfSets=nSets)),
		edgemode="directed"
		)
	edgeDataDefaults(E4.Vs,"Set") <- NA

	for (node in nodes(E4.Vs)) {
		for (Set in 1:nSets) {
			if (substr(node,Set,Set)=="1") { 
				nextnode <- node
				substr(nextnode,Set,Set) <- "0"
				E4.Vs <- addEdge(node,nextnode,E4.Vs )
				edgeData(E4.Vs,node,nextnode,"Set") <- Set
			}
		}
	}
	E4.Vs
}


matched.parentheses <- function(Vs) {
	# 0 is a (
	# 1 is a )
	# matched parentheses are returned as (), unmatched as []
	x <- strsplit(Vs,split="")[[1]]
	rhp <- match("1",x )
	while( any(x=="0" | x=="1")) {
		rhp <- match("1",x )
		if (is.na(rhp)) { #no more, so everything left is unmatched
			x[x=="0"] <- "["
		} else if (rhp==1) {# unmatched rp 
			 x [rhp] <- "]" 
		} else { # where was its lhp? 
			lphs <- which(x[1:(rhp-1)]=="0")
			if (length(lphs)>0) {
				lhp <- max(lphs)
				x [c(lhp,rhp)] <- c("(",")") 
			} else { # no lhp, so rhp is unmatched
				x[rhp] <- "]"
			}
		}
	}
	x
}

makeSCD <- function(Qn) {	
	# Greene-Kleitman cited in Ruskey Savage & Wagon
	nn <- nodes(Qn)
	SCD <- new("graphNEL",nodes=nn,edgemode="directed")
	edgeDataDefaults(SCD,"Set") <- NA
	mp <- sapply(nn,matched.parentheses)
	first.unmatched.1 <- apply(mp,2,function(x)match("]",x))
	first.unmatched.0 <- apply(mp,2,function(x)match("[",x))

	# nodes with no unmatched 1 and an unmatched 0
	no1un0 <- is.na(first.unmatched.1) & !is.na(first.unmatched.0)
	from <- nn[no1un0 ]
	fm0 <- first.unmatched.0[no1un0 ]

	to <- from
	substr(to,fm0,fm0) <- "1"
	SCD <- addEdge(from,to,SCD)
	edgeData(SCD,from,to,"Set") <- fm0
	
	# now repeat the 0-1 change
	while(length(to)>0) {
		mp2 <- sapply(to,matched.parentheses)
		first.unmatched.0.2 <- apply(mp2,2,function(x)match("[",x))
		# but just nodes with an unmatched 0
		un0 <-  !is.na(first.unmatched.0.2)
		from.2 <- to[un0]
		fm0.2 <- 	first.unmatched.0.2 [un0]
		to.2 <- from.2
		substr(to.2,fm0.2,fm0.2) <- "1"
		if (length(from.2)>0) {
			SCD <- addEdge(from.2,to.2,SCD)
			edgeData(SCD,from.2,to.2,"Set") <- fm0.2
		}
		to <- to.2
	}
	SCD		
}

addSedge <- function(G,from,to) {
	Set.changes <- rep(NA,length(from))
	for (i in seq_along(from)) {
		Set.changes[i]  <- which(strsplit(from[i],split="")[[1]]!=strsplit(to[i],split="")[[1]])
		if (length(Set.changes[i])!=1) {
			stop(sprintf("Number of Set changes is not one from %s to %s\n",from[i],to[i]))
		}
	}
	G <- addEdge(from,to,G)
	edgeData(G,from,to,"Set") <- Set.changes
	G
}


addcovers <- function(QSCD,layout="radial") {
	# QSCD is  disconnected collection of symmetric chain decompositions
	# of the hypercube Q[n]. We want to join them together to make
	# a maximal spanning subgraph of Q[n]
	# Moreover we compute coordinates for each point in an embedding
	# that is guaranteed to be planar
	# heads of chains
	incount <- sapply(inEdges(QSCD),length)
	heads <- nodes(QSCD)[ (incount==0)]
	headcovers <- sapply(heads,function(x){
		last.1 <- which(strsplit(x,split="")[[1]]==1)
		if (length(last.1)>0) {
			last.1 <- max(last.1)
			substr(x,last.1,last.1)	<- "0"
		}
		x
	})
	from <- headcovers[heads!=headcovers]
	to <- heads[heads!=headcovers]

	res <- QSCD
	res <- addSedge(res,from,to)


	head.to.tail <- acc(QSCD,heads)
	head.to.tail <- sapply(names(head.to.tail),function(n){
		s <-0; names(s)<-n
		c(s,head.to.tail[[n]])
	})
	chain.length <- sapply(head.to.tail,length)+1
	tails <- sapply(head.to.tail,function(x)names(x)[which.max(x)])
	source.node <- names(chain.length)[which.max(chain.length)]
	sink.node <- tails[source.node]
	inner.tails <- tails[tails!=sink.node]
	inner.heads <- names(inner.tails)
	inner.covers <- headcovers[match(inner.heads,heads)]
	cover.tails <- tails[inner.covers]
	res <- addSedge(res,inner.tails,cover.tails)
	# finally, we want to ensure the symmetric chains are ordered 
	# in such a way that the graph is planar
	# build a little tree of the heads of the chains
	headgraph  <- subGraph(heads,res)
	# visit it in depthfirst order
	dfs.order <- dfs(headgraph)$discovered
	heads.order <- match(heads,dfs.order)
	names(heads.order) <- heads
	# how deep did we have to go?
	heads.height <- acc(headgraph,source.node)[[1]]
	hs <- 0;names(hs)<- source.node
	heads.height <- c(hs,heads.height)
	# then apply that order to every member of the chain
	nodeDataDefaults(res,"x") <- NA
	nodeDataDefaults(res,"y") <- NA

	nChains <- length(heads)
	longestChain <- length(head.to.tail[[source.node]])
	for (h in names(head.to.tail)) {
		SCDorder <- as.numeric(heads.order[h])
		chain <- names(head.to.tail[[h]])
		r <- longestChain-as.numeric(seq_along(chain)+heads.height[h])
		if (layout=="radial") {
			angles <- seq(0,nChains-1)*(2*pi)/nChains
			x <- r * cos(angles[SCDorder])
			y <- r * sin(angles[SCDorder])
		} else {
			x <- SCDorder
			y <- longestChain-r
		}
		nodeData(res,chain,"x") <- x
		nodeData(res,chain,"y") <- y	
	}


	res
}

makePMSGn <- function(nSets) {
	Qn <- makeQn(nSets) # general hypercube Qn
	QSCD <- makeSCD(Qn) # symmetric chain decomposition gives a planar spanning subgraph
	QPSM <- addcovers(QSCD) # add covers makes it monotone
	QPSM
}



plotxygraph <- function(G) {
	nn <- nodes(G)
	x <- unlist(nodeData(G,attr="x"));names(x)<- nn
	y <- unlist(nodeData(G,attr="y"));names(y)<- nn

	grid.newpage()
	pushViewport(dataViewport(x,y,name="plotxygraph") )
	ee <- edges(G)
	eft <- do.call(rbind,lapply(names(ee),function(en){ees<- ee[[en]];cbind(rep(en,length(ees)),ees)}))
	row.names(eft) <- 1:nrow(eft)
	eft <- data.frame(eft,stringsAsFactors=FALSE)
#	eft$Edge <- paste(eft[,1],eft[,2],sep="|")
	eft.x <- cbind((x[eft[,1]]),(x[eft[,2]]))
	eft.y <- cbind((y[eft[,1]]),(y[eft[,2]]))
	eft$Set <- unlist(edgeData(G,from=eft[,1],to=eft[,2],attr="Set"))
	eft$Set <- as.numeric(as.factor(eft$Set))
	eft$Colour <- c("red","green","black","blue","brown")[eft$Set]
	grid.segments(x0=eft.x[,1],x1=eft.x[,2],y0=eft.y[,1],y1=eft.y[,2],gp=gpar(col=eft$Colour),default.units="native")
	grid.circle(x=x,y=y,r=0.05,gp=gpar(fill="white"),default.units="native")
	grid.text(x=x+0.1,y=y,gp=gpar(fill="white",cex=0.5),default.units="native",label=nn,just="left")
}


#PMSG4 <- makePMSGn(4)
#plotAWFE(Qn )
#plotAWFE(QSCD)
#plotAWFE(QPSM )
#plotAWFE(PMSG4)
#plotxygraph (QPSM)

	