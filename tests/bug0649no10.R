library(Vennerable)
setList <- strsplit(month.name,split="")
names(setList) <- month.name
VN3 <- VennFromSets( setList[1:3])
V2 <- VN3[,c("January","February"),]
V2.no10 <- V2
Weights(V2.no10)["10"] <- 0

#debug(joinEdgesInDrawing)
C2 <- compute.Venn(V2.no10,doWeights=FALSE,doEuler=TRUE,type="circles")
plot(C2)

V2.no01 <- V2
Weights(V2.no01)["01"] <- 0
C2 <- compute.Venn(V2.no01,doWeights=FALSE,doEuler=TRUE,type="circles")

