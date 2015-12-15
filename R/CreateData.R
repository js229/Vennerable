# This function is not called as normal part of package use (& it is not exported)
# It creates (or updates) the VennDiagrams.rda data file
# ...takes ~1 hr to run

# The package developer (only) should source() this file in the source package directory
# and then run
# buildVennDiagrams()


buildVennDiagrams <- function(rebuild=FALSE) {
# using functions defined in AWFE.Rnw


dataDir <- "./data"
dataFile <- file.path(dataDir,"VennDiagrams.rda")
if (!file.exists(dataFile)) {
	VennDiagrams <- list()
	save(VennDiagrams,file=dataFile,compress="xz")
}

SetBoundaries <- list()
SetBoundaries[["AWFE"]] <- makeAWFESets(10,type="AWFE",hmax=.58) # only works up to n=8
SetBoundaries[["AWFEscale"]] <- makeAWFESets(7,type="AWFEscale",hmax=.6) # can't make it work out to n=9
SetBoundaries[["cog"]] <- makeAWFESets(9,hmax=.6,type="cog")
SetBoundaries[["battle"]] <- makeAWFESets(9,type="battle")

saveTD <- function(type,rebuild=FALSE,Setlist ,dataFile) {
	cat(sprintf("Checking type %s\n",type))
	loaded <- load(dataFile)
	stopifnot("VennDiagrams" %in% loaded)
	if (is.null(VennDiagrams[[type]])) { VennDiagrams[[type]] <- list() }
	if (rebuild | length(VennDiagrams[[type]])==0) {
		VennDiagrams[[type]][[1]]  <- makeAWFE(n=1,Setlist=Setlist)
	}
	if (rebuild) { nstart <- 2 } else { 
		if (length(Setlist) <= length(VennDiagrams[[type]])) {
			return()
		} else {
			nstart <- length(VennDiagrams[[type]])+1
		}
	}
	if (nstart>length(Setlist)) { stop("Don't have enough sets") }
	for (n in nstart:length(Setlist)) { 
		TD <- makeAWFE(n,AWFEnminus1=VennDiagrams[[type]][[n-1]],Setlist=Setlist)	
		.validateDrawing(TD)
		VennDiagrams[[type]][[n]] <- TD
		save(VennDiagrams,file=dataFile,compress="xz")
	}
}



saveTD("battle",Setlist =SetBoundaries[["battle"]],dataFile=dataFile,rebuild=rebuild)
saveTD(type="AWFEscale",Setlist =SetBoundaries[["AWFEscale"]],dataFile=dataFile,rebuild=rebuild)
saveTD("AWFE",Setlist =SetBoundaries[["AWFE"]][1:10],dataFile=dataFile,rebuild=rebuild)

	loaded <- load(dataFile)
	stopifnot("VennDiagrams" %in% loaded)
	vde <- VennDiagrams[["ellipses"]]
	if (is.null(vde)) {
		 vde  <- list()
	}
	if (rebuild | length(vde)<4) {
		vde[[4]] <- NULL
		E4 <- make.E4(dx=0.05)
		.validateDrawing(E4)
		vde[[4]] <- E4
		VennDiagrams[["ellipses"]] <- vde
		save(VennDiagrams,file=dataFile,compress="xz")
	}
}


