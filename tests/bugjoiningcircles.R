library(Vennerable)
VD2 <- compute.Venn(Venn(n=2))
VD3 <- newTissueFromCircle (centre.xy =c(2,0), radius=.6,Set=3)
VD23 <- VD2
VD23@faceList <- c(VD2@faceList,VD3@faceList)
VD23@edgeList <- c(VD2@edgeList,VD3@edgeList)
VD23@setList <- c(VD2@setList,VD3@setList)

library(grid)
grid.newpage()
pushViewport(plotViewport(c(1,1,1,1)))
makevp.eqsc(c(-2,3),c(-2,2));grid.xaxis();grid.yaxis()
PlotSetBoundaries(VD23)

drawing <- VD23
innerFaceName <- "1"
.create.edge.joining.faces(drawing,"DarkMatter","1")
