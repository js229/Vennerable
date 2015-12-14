\name{VennSetSetLabels}
\Rdversion{1.1}
\alias{VennSetSetLabels}
\alias{VennGetSetLabels}
\alias{VennSetFaceLabels}
\alias{VennGetFaceLabels}
\alias{VennSetUniverseRange}
\alias{VennGetUniverseRange}
\alias{VennGetUniverseRange,VennDrawing-method}
\alias{VennSetUniverseRange,VennDrawing-method}
\title{
Set and get annotation labels and graphical styles
}
\description{
Embedded within a \code{VennDiagram} object are instructions for displaying annotation
for the sets and their faces. Retrieve them so they can be edited
and re-embedded for subsequent display.
}
\usage{
VennGetSetLabels(object)
VennSetSetLabels(object, SetLabels)
VennGetFaceLabels(object)
VennSetFaceLabels(object, FaceLabels)
VennGetUniverseRange(object)
VennSetUniverseRange(object, universe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{VennDrawing}
}
  \item{SetLabels}{
A data.frame. The format of this data.frame may change. 
}
  \item{FaceLabels}{
A data.frame. The format of this data.frame may change. 
}
  \item{universe}{
An \eqn{nx2} matrix defining \eqn{x} and \eqn{y} coordinates of the bounding universe. 
Currently assumed to have 2 rows by much of the code, 
corresponding to the lower left and upper right corners of a bounding rectangle}
}
\author{
Jonathan Swinton (jonathan@swintons.net)
}
\examples{
data(StemCell)
Vstem <- Venn(StemCell)
Vstem3 <- Vstem[,c("OCT4","SOX2","NANOG")]
Cstem3 <- compute.Venn(Vstem3,doWeights=TRUE)
plot(Cstem3)
# don't like the default position of 'SOX2'?
SetLabels <- VennGetSetLabels(Cstem3)
SetLabels[SetLabels$Label=="SOX2","x"] <- 12
Cstem3 <- VennSetSetLabels(Cstem3,SetLabels)
grid::grid.newpage()
plot(Cstem3)
}
\keyword{ graphs }
