\name{PlotVennGeometry}
\Rdversion{1.1}
\alias{PlotVennGeometry}
\title{
Display computed Venn diagrams
}
\description{
This is the function called by the \code{plot} generic for objects of class \code{VennDrawing}.
}
\usage{
PlotVennGeometry(C3, gpList, show = list(FaceText = "weight"))
}
\arguments{
  \item{C3}{
An object of class \code{VennDrawing}, as produced by \code{compute.Venn}
}
  \item{gpList}{
A list of lists of graphical parameters. The \code{VennThemes} function produces a suitable list; see its man page for more details.
}
  \item{show}{
A list of options controlling what is shown. Defaults are \code{universe = TRUE, Sets = TRUE, SetLabels = TRUE, 
        DarkMatter = FALSE, Faces = TRUE, FaceText = "weight"}. If \code{universe = TRUE} the bounding region, currently always a rectangle,
of the diagram is shown. If \code{Sets = TRUE}, the boundaries of the Sets are shown. If \code{SetLabels = TRUE} the names of
the Sets are plotted according to the instructions in \code{VennGetSetLabels(C3)}. If \code{Faces = TRUE}
individual faces are filled. \code{FaceText} is a character vector which may contain any of 
\code{c("weight","signature","sets","elements")}. Currently \code{DarkMatter} must be \code{FALSE}.
}
}
\details{
The algorithm for displaying face annotation is pretty crude. The annotation corresponding to each member of the
\code{FaceText} vector is displayed on a single line. \code{weight} and \code{signature} display the associated 
weight and signature for the face. \code{sets} produces a concatenated version of the names of the sets in the face,
while \code{elements} does the same for the elements; both of these only really work for single character set names or elements.
}
\author{
Jonathan Swinton (jonathan@swintons.net)
}
\examples{
library(Vennerable)
setList <- strsplit(month.name,split="")
names(setList) <- month.name
C3 <- compute.Venn(Venn(setList[1:3]))
PlotVennGeometry(C3,show=list(FaceText="signature"),gp=VennThemes(C3,colourAlgorithm="signature"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs}
