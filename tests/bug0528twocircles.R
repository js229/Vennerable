library(Vennerable)
centre.xy <- c(0,-2)
VDC1b <- newTissueFromCircle(centre.xy,radius=2,Set=1,nodes=4)
VDC2b <- newTissueFromCircle(centre.xy+c(0,3),radius=1,Set=2)
TN2b <- (addSetToDrawing(VDC1b,VDC2b))
TN2b
.validateDrawing(TN2b)
 
 
 

 
