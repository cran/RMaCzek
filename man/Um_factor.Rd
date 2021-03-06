\name{Um_factor}
\alias{Um_factor}
\title{Calculate the Um factor}
\usage{
Um_factor(
  distMatrix,
  order = NULL,
  matrix_conversion_coefficient = 1,
  inverse_um = TRUE
)
}
\arguments{
\item{distMatrix}{a 'dist' object, matrix of distances between observations.}

\item{order}{a vector, if NULL, then the value of the factor is calculate for the distance matrix as is, otherwise the rows and columns are reordered according to the vector order.}

\item{matrix_conversion_coefficient}{numeric, value to be added to the distances, so that a division by 0 error is not thrown.}

\item{inverse_um}{logical, if TRUE, then the negative is returned. Default TRUE as the function is called in the genetic algorithm maximization procedures.}
}
\value{
The function returns a numeric value equalling the Um_factor.
}
\description{
The function calculates the Um factor associated with an ordering of the rows and columns of a distance matrix. Lower values indicate a better grouping of similar objects. This was the original objective function proposed in the MaCzek program for producing Czekanowski's Diagram.
}

\references{
K. Bartoszek and A. Vasterlund (2020).
"Old Techniques for New Times": the RMaCzek package
for producing Czekanowski's diagrams.
Biometrical Letters 57(2):89-118.

A. Soltysiak and P. Jaskulski (1999). 
Czekanowski's diagram. a method of multidimensional clustering. 
New Techniques for Old Times. CAA 98. Computer Applications and
Quantitative Methods in Archaeology. Proceedings of the 26th Conference, Barcelona,
March 1998, number 757 in BAR International Series, pages 175-184, Oxford.
}



\examples{
# Set data ####
x<-mtcars

mD<-stats::dist(scale(x))
mCz<-czek_matrix(x)
Um_factor(mD)
Um_factor(mD,order=attr(mCz,"order"))
}
