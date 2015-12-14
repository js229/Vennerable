\name{VennDrawing-class}
\Rdversion{1.1}
\docType{class}
\alias{VennDrawing-class}
\alias{Areas,VennDrawing-method}
\alias{IntersectionMidpoints,VennDrawing-method}
\alias{plot,VennDrawing,missing-method}
\alias{PlotUniverse,VennDrawing-method}
\alias{show,VennDrawing-method}
\alias{UniverseRange,VennDrawing-method}
\alias{VisibleRange,TissueDrawing-method}
\alias{VisibleRange-methods}         
\alias{VennSetSetLabels-methods}
\alias{VennSetSetLabels,VennDrawing-method}
\alias{VennGetSetLabels,VennDrawing-method}
\alias{VennSetFaceLabels,VennDrawing-method}
\alias{VennGetFaceLabels,VennDrawing-method}


\title{Class "VennDrawing" }
\description{A geometrical representation of a Venn diagram suitable for display}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VennDrawing",V,TD)}
where \code{V} is an object of class \code{Venn} and \code{TD} is an object of class \code{TissueDrawing} which
has representations for all the Sets and Faces corresponding to the sets and intersections of the \code{V} object.

As well as bringing together the geometric and set-theoretic properties of the Venn diagram,
this class stores information about the universe in which to display the diagram and provides
a number of methods for plotting the final figure.
}
\section{Slots}{
	 \describe{
    \item{\code{universe}:}{A \code{matrix} with the \eqn{x} and \eqn{y} coordinates of the bounding universe. Currently assumed to have 2 rows by much of the code, corresponding to the lower left and upper right corners of a bounding rectangle}
    \item{\code{SetLabels}:}{A data.frame encoding the positions of labels for the Sets}
      \item{\code{FaceLabels}:}{A data.frame encoding the positions of labels for the Faces}
  	\item{\code{setList}:}{See \code{\link{TissueDrawing-class}}  }
    \item{\code{nodeList}:}{See \code{\link{TissueDrawing-class}} }
    \item{\code{recentChanges}:}{See \code{\link{TissueDrawing-class}}}
    \item{\code{edgeList}:}{See \code{\link{TissueDrawing-class}} }
    \item{\code{faceList}:}{See \code{\link{TissueDrawing-class}}}
    \item{\code{faceSignature}:}{See \code{\link{TissueDrawing-class}} }
    \item{\code{IndicatorWeight}:}{See \code{\link{Venn-class}}}
    \item{\code{IntersectionSets}:}{See \code{\link{Venn-class}} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{TissueDrawing}"}, directly.
Class \code{"\linkS4class{Venn}"}, directly.
Class \code{"\linkS4class{TDEdgeList}"}, by class "TissueDrawing", distance 2.
Class \code{"\linkS4class{TDFaceList}"}, by class "TissueDrawing", distance 2.
}
\section{Methods}{
  \describe{
    \item{Areas}{\code{signature(object = "VennDrawing")}: Returns the area of each named face }
    \item{IntersectionMidpoints}{\code{signature(object = "VennDrawing")}: ... }
    \item{plot}{\code{signature(x = "VennDrawing", y = "missing")}: ... }
    \item{PlotUniverse}{\code{signature(object = "VennDrawing")}: ... }
    \item{show}{\code{signature(object = "VennDrawing")}: ... }
    \item{UniverseRange}{\code{signature(object = "VennDrawing")}: ... }
    \item{VisibleRange}{\code{signature(object = "VennDrawing")}: ... }
	 }
}
\author{ Jonathan Swinton (jonathan@swintons.net) }


\examples{
showClass("VennDrawing")
}
\keyword{classes}
