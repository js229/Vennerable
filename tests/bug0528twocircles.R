library(Vennerable)
centre.xy <- c(0,-2)
VDC1b <- Vennerable:::newTissueFromCircle(centre.xy,radius=2,Set=1,nodes=4)
VDC2b <- Vennerable:::newTissueFromCircle(centre.xy+c(0,3),radius=1,Set=2)
TN2b <- (Vennerable:::addSetToDrawing(VDC1b,VDC2b))
TN2b
Vennerable:::.validateDrawing(TN2b)
 
 
 

 
