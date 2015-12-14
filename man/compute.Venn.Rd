\name{compute.Venn}
\Rdversion{1.1}
\alias{compute.Venn}
\title{
Calculate the geometry of a Venn or Euler diagram.
}
\description{
Given an object of class \code{Venn}, compute an appropriate diagram.
}
\usage{
compute.Venn(V, doWeights = TRUE, doEuler = FALSE, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
An object of class \code{Venn}
}
  \item{doWeights}{
Logical
}
  \item{doEuler}{
Logical
}
  \item{type}{
One of \code{"circles","squares","triangles","ellipses","AWFE","ChowRuskey","AWFE"}.
}
}
\details{
For more details see the package vignette with \code{vignette("Venn")}
}
\value{
An object of class \code{VennDrawing}
}
\author{
Jonathan Swinton (jonathan@swintons.net)
}
\examples{
setList <- strsplit(month.name,split="")
names(setList) <- month.name
compute.Venn(Venn(setList[1:3]))
}
\keyword{ graphs }
