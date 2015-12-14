\name{TissueDrawing-class}
\Rdversion{1.1}
\docType{class}
\alias{TissueDrawing-class}
\alias{CircleDrawing-class}
\alias{SquareDrawing-class}

\alias{PlotFaces,TissueDrawing-method}
\alias{PlotNodes,TissueDrawing-method}
\alias{PlotSetBoundaries,TissueDrawing-method}
\alias{show,TissueDrawing-method}

\title{Class "TissueDrawing" }
\description{A class used internally to implement polygon-polgyon clipping
}
\section{Objects from the Class}{
The class describes diagrams which are constructed by superposing simple polygons (ie those with no holes or self-intersections), called Sets. 
If the Sets are filled with semi-transparent colours, the resulting diagram will have different
regions or Faces depending on the intersections between the sets. The boundaries of these regions will not necessarily
be simple. (For example one circle may have two other circles nested inside it).
A simple model for this process is superpsosing Sets made out of tissue paper, hence the name of the class, which was nearly called
a MAGDrawing after the perspex and lightbox exhibit in the children's gallery of the Manchester Art Gallery.
which also acted as a prototype.

Points where the Sets intersect are always identified and named within the object.
Sets which partially share a common boundary are defined not to intersect on that common boundary,
so there are only a finite number of intersection points. Each Face contains at least one
named point so additional named points will be introduced for each Face that has no other intersections.
In addition, named points may be introduced to make some of the hidden-line alogorithms easier.
These additional named points can be removed (sometimes) with the function \code{remove.nonintersectionpoints}
which leaves the remaining points as the nodes of the graph defined by the Set intersections.

Set and Face boundaries are made up of named edges, each of which starts and finishes at the named points.
Each edge can be either a polygonal segment or a sector of a circle. All of the Faces
in a \code{TissueDrawing} are connected to each other (ie the graph is connected), if necessary
by adding invisible edges which are not part of the boundary of any set.



Objects can be created by calls to the functions \code{newTissueFromPolygon()} or 
\code{newTissueFromCircle}, and combined with the function \code{addSetToDrawing}).


}
\section{Slots}{
	 \describe{
    \item{\code{setList}:}{Named \code{"list"}, with names given by the name of each Set. Each entry is a character vector of edge names }
    \item{\code{nodeList}:}{Named \code{"list"}, with a 1x2 numeric matrix of $(x,y)$ coordinates for each named point. }
    \item{\code{recentChanges}:}{Nasty bodge slot required by internal algorithms }
    \item{\code{edgeList}:}{A \code{TDEdgeList}, a named list of \code{TDedgeDrawn} objects defining the geometry of each edge }
    \item{\code{faceList}:}{A \code{TDFaceList}, a named \code{"list"}, with one element for each correspondingly named face, listing the edges of that face in clockwise order }
    \item{\code{faceSignature}:}{Named  \code{"list"}, with names the same as \code{faceList}, giving the signature for the face.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{TDEdgeList}"}, directly.
Class \code{"\linkS4class{TDFaceList}"}, directly.
}
\section{Methods}{
  \describe{
    \item{PlotFaces}{\code{signature(drawing = "TissueDrawing")}: Displays the faces of the drawing on the current graphics device }
    \item{PlotNodes}{\code{signature(drawing = "TissueDrawing")}: Displays the named points of the drawing }
    \item{PlotSetBoundaries}{\code{signature(drawing = "TissueDrawing")}: Displays the Set boundaries of the drawing }
    \item{show}{\code{signature(object = "TissueDrawing")}: Prints a summary of the drawing }
	 }
}
\author{Jonathan Swinton, jonathan@swintons.net}
\note{Adding even two overlapping simple polygons together can result in multiple
intersections regions. Finding the boundaries of these intersections,
keeping track of which have come from which combination of sets, and, in particular, coping with 
faces embedded inside outher faces is the problem of polygon-polygon clipping. 
This class is designed to solve that problem for the combinations of Sets used
in the \code{Vennerable} package. It is not designed as a general clipping tool, and 
can be very slow, particular for clipping of general polygons. For example, the four ellipse
plot is currently implemented as a combination of many short-segment polygons and several minutes
to run. 

Care is taken to keep track of the outside boundary of all of the sets, or equivalently the inside
boundary of the 'dark matter', the Face defined as the region where none of the Sets intersect.
Nevertheless this portion of the code is currently quite buggy and this can lead to errors in some computed diagrams.

 }

\examples{
showClass("TissueDrawing")
}
\keyword{classes}