\name{plot.czek_matrix}
\alias{plot.czek_matrix}
\title{Produce a Czekanowski's Diagram}
\usage{
\method{plot}{czek_matrix}(x, values = NULL, type = "symbols",
  plot_pch = NULL, plot_cex = 1.5, label.cex = 0.6,
  plot_title = "Czekanowski's diagram", legend = FALSE, axis = TRUE,
  ...)
}
\arguments{
\item{x}{a matrix with class czek_matrix.}

\item{values}{specifies the color or the size of the symbols in the graph. The standard setting is a grey scale for a color graph and a vector with the values 2,1,0.5,0.25 and 0 for a graph with symbols.}

\item{type}{specifies if the graph should use color or symbols. The standard setting is symbols.}

\item{plot_pch}{specifies which symbols the graph should use. The standard setting is 19, which is a black circle.}

\item{plot_cex}{specifies the size of the cells in a color graph. The standard setting is 1.5.}

\item{label.cex}{specifies the size of the labels for the objects. The standard setting is 0.6.}

\item{plot_title}{specifies the main title in the graph.}

\item{legend}{specifies if a legend should be included or not. The standard setting is that the legend will not be included.}

\item{axis}{specifies if the labels for the objects should be included. The standard setting is that the labels are included.}

\item{...}{specifies further parameters that can be passed on to the plot function.}
}
\description{
This is a function that can produce Czekanowski's Diagram.
}

\author{
Albin Vasterlund

Maintainer: Krzysztof Bartoszek <krzbar@protonmail.ch>
}

\references{
K. Bartoszek and A. Vasterlund (2020).
"Old Techniques for New Times": the RMaCzek package
for producing Czekanowski's diagrams.
Biometrical Letters 57(2):89-118.

J. Czekanowski (1909).
Zur Differentialdiagnose der Neandertalgruppe.
Korespondentblatt der Deutschen Gesellschaft fur Anthropologie, 
Ethnologie und Urgeschichte, XL(6/7):44-47,

A. Soltysiak and P. Jaskulski (1999). 
Czekanowski's diagram. a method of multidimensional clustering. 
New Techniques for Old Times. CAA 98. Computer Applications and
Quantitative Methods in Archaeology. Proceedings of the 26th Conference, Barcelona,
March 1998, number 757 in BAR International Series, pages 175-184, Oxford.
 
A. Vasterlund (2019). 
Czekanowski's Diagram: Implementing and exploring Czekanowski's Diagram with different seriation methods
Master thesis, Linkoping University
}


\examples{
# Set data ####
x<-czek_matrix(mtcars)

# Standard plot ############
plot(x)
plot.czek_matrix(x)

# Specify values ############
plot(x,values=c(1.5,1,0.75,0.25,0 ))
plot(x,values=grDevices::colorRampPalette(c("black","red","white"))(5))

# Specify type ############
plot(x,type = "symbols")
plot(x,type = "col")

# Specify plot_pch ############
plot(x,plot_pch = 15)

# Specify plot_cex ############
plot(x,type="col",plot_cex = 1)

# Specify plot_cex ############
plot(x,label.cex = 0.45)

# Specify the main title ############
plot(x,plot_title = "Czekanowski's Diagram of mtcars")

# Add legend ############
plot(x,legend = TRUE)

# Remove axis name ############
plot(x,axis = FALSE)

# Change additional settings to the plot function ############
plot(x,col.main="blue",font.main=9,cex.main=2)
}

\keyword{hplot}
