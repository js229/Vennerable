library(Vennerable)
setList <- strsplit(month.name,split="")
names(setList) <- month.name
Vmonth3 <- VennFromSets( setList[1:3])
Vmonth2 <- Vmonth3[,c("January","February"),]
Vmonth2.no11 <- Vmonth2
Weights(Vmonth2.no11)["11"] <- 0
V <- Vmonth2.no11
#undebug(remove.nonintersectionpoints)
#undebug(joinEdgesInDrawing)
C2 <- compute.Venn(V,doWeights=FALSE,doEuler=TRUE,type="circles")
