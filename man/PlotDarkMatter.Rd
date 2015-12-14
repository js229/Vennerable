\name{PlotDarkMatter}
\Rdversion{1.1}
\alias{PlotDarkMatter}
\title{
Plots dark matter
}
\description{
Fills the area of a \code{VennDrawing} universe that is not occupied by any set.
}
\usage{
PlotDarkMatter(VD)
}
\arguments{
  \item{VD}{
An object of class \code{VennDrawing}
}
}
\details{
This works by filling the entire universe with grey (not currently changeable) and then plotting
then filling the inside of the dark matter boundary with white.
}
\value{
Executed for its side effecrs
}
\author{
Jonathan Swinton (jonathan@swintons.net)
}
\examples{
VD <- compute.Venn(Venn(n=2))
PlotDarkMatter(VD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
