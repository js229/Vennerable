\name{StemCell}
\Rdversion{1.1}
\alias{StemCell}
\docType{data}
\title{
Stem cell data
}
\description{
Data on the genes which bind to each of three transcription factors (and a fourth which may also be a transcription factor)
}
\usage{data(StemCell)}
\format{
A named list with four elements named \code{OCT4}, \code{SOX2}, \code{NANOG} and \code{E2F4}.
Each element is a character vector containing the names of 600-1700 associated genes.
}
\source{
Taken from supplementary information in Boyer et al, "Core transcriptional regulatory circuitry in human embryonic stem cells",
\emph{Cell}, 2005, 122:947--956, \code{10.1016/j.cell.2005.08.020},\url{http://www.ncbi.nlm.nih.gov/pubmed/16153702}.
}

\examples{
library(Vennerable)
data(StemCell)
plot(Venn(Sets=StemCell),type="squares",doWeights=FALSE)
}
\keyword{datasets}
