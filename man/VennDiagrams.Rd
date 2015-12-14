\name{VennDiagrams}
\Rdversion{1.1}
\alias{VennDiagrams}
\docType{data}
\title{
Predrawn Venn diagrams.
}
\description{
Precomputed data on the geometrical structure of various types of Venn diagrams.
}
\usage{data(VennDiagrams)}
\format{A list with elements named \code{battle},\code{AWFEscale},\code{AWFE},
and \code{ellipses}. Each element is itself a list in which the nth element, if not \code{NULL}, is an object of 
class \code{TissueDrawing} representing an example of an unweighted Venn diagram on n sets of the given type.
}
\details{
This is really a device to cache the diagrams and is not for direct use. Use eg \code{compute.Venn(Venn(n=...),type=...} which returns a plottable object of class \code{VennDrawing} instead.
}
\source{
Created by the \code{buildVennDiagrams()} function within \code{Vennerable}
}
\examples{
data(VennDiagrams)
TD <- VennDiagrams[["ellipses"]][[4]]
class(TD)
}
\keyword{datasets}
