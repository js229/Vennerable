\name{Venn}
\Rdversion{1.1}
\alias{Venn}
\alias{Weights}
\alias{Weights<-}
\title{Construct intersections of sets}
\description{
Given a collection of Sets, \code{Venn} will compute all possible 
combinations of intersections and return an object of class \code{Venn}, storing the combinations as well
as the number of elements in each intersection (the 'weights'). 
It will also accept a specification of the weights directly,
or simply of the number of sets desired.
}
\usage{
Venn(Sets, Weight, SetNames, numberOfSets)
Weights(object)
Weights(object) <- value
}
\arguments{
  \item{Sets}{A list, each element of which can be a list defining the members of a Set}
  \item{Weight}{A named vector of weights to associate with each possible combination of Set intersections.
	Ignored if \code{Sets} is present.}
  \item{SetNames}{A character vector of names for each Set. Ignored if \code{Sets} present and has names}
  \item{numberOfSets}{An integer for the number of Sets. Ignored if  \code{Sets} or \code{SetNames} present.}
  \item{object}{An object of class \code{Venn}}
  \item{value}{An object like \code{Weight}}
}
\value{An object of class \code{Venn}}
\author{Jonathan Swinton (jonathan@swintons.net)}
\examples{
Venn(n=3)
data(StemCell)
w <- Venn(Sets=StemCell[1:2])
Weights(w)
# nb order of intersection subsets not guaranteed
Weights(w) <- 1:4
Venn(SetNames=letters[1:5])

}
\keyword{ graphs}
