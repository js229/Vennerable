#warning("Entering 00Venn")
#**************************************************************************
#
# # $Id: 00Venn.R,v 1.19 2007/10/16 10:59:52 js229 Exp $

# $Log: 00Venn.R,v $

setClass("Venn",
	representation( IndicatorWeight="matrix",IntersectionSets="list")
)

Venn <- function(Sets,Weight,SetNames,numberOfSets) {
	if (!missing(Sets)) {
		if (is.null(names(Sets)) & !missing(SetNames)) {
			names(Sets) <- SetNames
		}
		return(VennFromSets(Sets))
	}
	if(missing(numberOfSets)) {
		numberOfSets <- if (missing(SetNames)) 0  else length(SetNames)
	}
	if (missing(SetNames)) {
		SetNames <- seq_len(numberOfSets)
	}


	if(numberOfSets==0) {
		Indicator <- matrix(nrow=0,ncol=0)
	} else {
		Indicator <- (data.matrix(do.call(expand.grid,lapply(seq(1,length=numberOfSets),function(x){c(0,1)}))))==1
	}
	rownames(Indicator) <- apply((Indicator),1,function(x){paste(as.numeric(x),collapse="")})
	colnames(Indicator) <- SetNames

	if(missing(Weight)) {
		Weight <- rep(1,nrow(Indicator))
	} else if( length(Weight) >0 & length(Weight) < nrow(Indicator) & is.null(names(Weight)) ) {
		stop("Weight length does not match number of intersections")
	}
	if (is.null(names(Weight))) {
		names(Weight) <- rownames(Indicator)
	}
	
	IndicatorWeight <- cbind(Indicator,.Weight=Weight[rownames(Indicator)])

	new("Venn",IndicatorWeight =IndicatorWeight )
}

VennFromSets <- function(setList) {
	SetNames <- names(setList) 
	if (is.null(SetNames)) { SetNames <- seq_along(setList) }
	VN <- Venn(SetNames=SetNames)
	VIsig <- VennSignature(VN)
	IntersectionSets <- sapply(1:length(VIsig),function(x)NULL)
	names(IntersectionSets) <- VIsig
	universe <- NULL
	for ( iset in setList) { universe <- union (universe,iset) }
	for (element in universe) {
		sig <- as.numeric(!is.na(sapply(setList,match,x=element)))
		sig <- paste(sig,collapse="")
		IntersectionSets[[sig]] <- union(IntersectionSets[[sig]],element)
		}
	Weights <- sapply(IntersectionSets,length)
	
	Vres <- Venn(SetNames=SetNames,Weight=Weights)
	Vres@IntersectionSets <- IntersectionSets
	Vres
}

setMethod("show","Venn",function(object){
	cat(sprintf("A Venn object on %d sets named\n",NumberOfSets(object)))
	cat(paste(VennSetNames(object),collapse=","),"\n")
	show(Weights(object))
})

#setGeneric("NumberOfSets",function(object){standardGeneric("NumberOfSets")})
#setMethod("NumberOfSets","Venn",function(object){ncol(object@IndicatorWeight)-1})
NumberOfSets <- function(object){ncol(object@IndicatorWeight)-1}
Indicator <- function(object){
	object <- as(object,"Venn")
	object@IndicatorWeight[,-ncol(object@IndicatorWeight),drop=FALSE]}

#setGeneric("SetNames",function(object){standardGeneric("SetNames")})

VennSignature<- function(object){
	ind <- Indicator(object)
	inn <- apply(ind,1,paste,collapse="")
	inn
}
Weights <- function(object) {
	V <- as(object,"Venn")
	wght <- V@IndicatorWeight[,".Weight"]
	names(wght) <- VennSignature(V)
	wght
}

"Weights<-" <- function(object,value) {
	V <- as(object,"Venn")
	VS <- VennSignature(V)
	if (is.null(names(value))) {
		names(value) <- VS
	}
	value <- value[match(VS,names(value))]
	object@IndicatorWeight[,".Weight"] <- as.numeric(value)
	object
}

dark.matter.signature <- function(object) {
	V <- as(object,"Venn")
	VS <- VennSignature(V)
	VS[regexpr("1",VS)<0]
}	


#setMethod("SetNames","Venn", function(object) {cn <- colnames(object@IndicatorWeight);cn[cn!=".Weight"] })
VennSetNames <-  function(object) {cn <- colnames(object@IndicatorWeight);cn[cn!=".Weight"] }
setMethod("[","Venn", function(x,i,j,...,drop) {
	if (!missing(i)) {
		stop("Can't subset on rows")
	}
	if (!missing(j)) {
		Indicator <- x@IndicatorWeight
		Signature <- apply(Indicator [,setdiff(colnames(Indicator ),".Weight")],1,paste,collapse="")
		Indicator.df <- data.frame(Indicator)
		Indicator.df$Signature <- Signature

		newIndicator <- aggregate(Indicator[,".Weight",drop=FALSE],
			by=data.frame(Indicator[,j,drop=FALSE]),FUN=sum)
		for (col in setdiff(colnames(newIndicator ),".Weight")) {
			newIndicator [,col] <- as.numeric(as.character(newIndicator [,col])) # sigh
		}
		x@IndicatorWeight <- data.matrix(newIndicator)	
	
		newIndicator.df <- newIndicator[,setdiff(colnames(newIndicator),".Weight")]
		newSignature <- apply(newIndicator.df,1,paste,collapse="")
		newIndicator.df <- data.frame(newIndicator.df)
		newIndicator.df$newSignature <- newSignature
		newIndicator.df

		IntersectionSets <- x@IntersectionSets
		if (length(IntersectionSets )>0) {
			dfm <- merge(newIndicator.df,Indicator.df)
			newsetNames <- split(dfm$Signature,dfm$newSignature)
			newIntersectionSets <- sapply(newsetNames,function(sigs)as.vector(do.call(c,IntersectionSets[sigs])))
			x@IntersectionSets <- newIntersectionSets 
		}
	}
	x
	})

.WeightVisible <- function(V) {# excluding  elements not in any sets
	wght <- V@IndicatorWeight
	wghtsum <- apply(wght[,-ncol(wght)],1,sum)
	sum(wght[wghtsum>0,ncol(wght)]) 
}
.WeightUniverse <- function(V) {
	wght <- Weights(V)
	sum(wght) 
}

################################

####################

compute.Venn <- function(V,doWeights=TRUE,doEuler=FALSE,type) {
	nSets <- NumberOfSets(V)
	if (nSets < 2) {
		stop("Not enough sets")
	} 
	if (missing(type)) {
		type <- if (nSets==2) {
			"circles"
		} else if (nSets==3) {
			"circles"
		} else if (nSets==4) {
			if (doWeights) "ChowRuskey" else "squares"
		} else {
			if (doWeights) "ChowRuskey" else "AWFE"	
		}
	}        
	C3 <-switch(type,
		AWFE=,AWFEscale=,battle=,cog=compute.AWFE(V,type=type),
		ChowRuskey=,compute.CR(V,doWeights),
		circles=
			if (nSets==2) { compute.C2(V,doWeights,doEuler) 
			}else if (nSets==3) {compute.C3(V,doWeights)
			} else { stop(sprintf("Type %s not implemented for %d sets",type,nSets))
			} ,
		squares=
		  	if (nSets==2) { compute.S2(V,doWeights,doEuler) 
			}else if (nSets==3) {compute.S3(V,doWeights)
			}else if (nSets==4) {compute.S4(V,doWeights)
		} else { stop(sprintf("Type %s not implemented for %d sets",type,nSets))
			} ,
		triangles= if (nSets==3) { compute.T3(V,doWeights) 
			} else { stop(sprintf("Type %s not implemented for %d sets",type,nSets))
			} ,
		ellipses= if (nSets==4) { compute.E4(V,doWeights) 
			} else { stop(sprintf("Type %s not implemented for %d sets",type,nSets))
			} 
	)	
	C3
}



plotVenn <- function(V,doWeights=TRUE,doEuler=FALSE,type,add=FALSE,
			show=list(FaceText="weight",Faces=TRUE),
		gpList){
	C3 <- compute.Venn(V,doWeights=doWeights,doEuler=doEuler,type=type)
	if (!add) {
		grid.newpage()
	}
	PlotVennGeometry(C3,gpList=gpList,show=show)
}

setGeneric("plot")
setMethod("plot", signature(x="Venn", y="missing"),function(x,y,...)plotVenn(V=x,...))

#plot.Venn <- plotVenn
#setMethod("plot",signature(x="Venn",y="missing"),function(x,y,...)plotVenn(V=x,...))


