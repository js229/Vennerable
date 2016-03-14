library(Vennerable)
setList <- strsplit(month.name,split="")
names(setList) <- month.name
Vempty2 <- VennFromSets( setList[c(4,5,6)])
TJAM <- compute.Venn(Vempty2)
Vennerable:::.validateDrawing(TJAM)

 
