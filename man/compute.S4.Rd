\name{compute.S4}
\Rdversion{1.1}
\alias{compute.S4}
\title{Venn diagrams on 4 sets}
\description{
Given a \code{Venn} object on 4 sets, returns an  \code{VennDrawing} object rendering it 
in the Lewis Carroll style.
}
\usage{
compute.S4(V, doWeights = FALSE, s = 0.25, likeSquares = TRUE)
}
\arguments{
  \item{V}{
A \code{Venn} object on 4 sets.}
  \item{doWeights}{
Must be \code{FALSE}. For weighted diagrams, use \code{compute.CR}.
}
  \item{s}{
If \eqn{s=0} then set boundaries may be superimposed; for \eqn{s>0} the larger two rectangles
are offset by \eqn{s}.
}
  \item{likeSquares}{
If \code{TRUE}, make the individual faces as square as possible; if \code{FALSE} remove additional corners and use triangular edges.
}
}
\details{
The placement of the face names is not ideal for some combinations with \eqn{s>0}.
}
\value{
An object of class \code{VennDrawing}
}
\author{
Jonathan Swinton (jonathan@swintons.net)
}
\examples{
compute.S4(Venn(n=4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
