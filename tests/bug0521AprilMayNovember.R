library(Vennerable)
setList <- strsplit(month.name,split="")
names(setList) <- month.name
Vempty2 <- VennFromSets( setList[c(4,5,11)])
TAMN <- compute.Venn(Vempty2)
.validateDrawing(TAMN )
#grid.newpage();plot(TAMN )
 
