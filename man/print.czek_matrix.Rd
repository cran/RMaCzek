\name{print.czek_matrix}
\alias{print.czek_matrix}
\title{Print a Czekanowski's Diagram}
\usage{
\method{print}{czek_matrix}(x, print_raw = FALSE, ...)
}
\arguments{
\item{x}{a matrix with class czek_matrix.}

\item{print_raw}{logical, if TRUE print out raw, as if the object was a matrix, 
    in particular this prints out the matrix itself, if FALSE (default) print out a summary.
    Furthermore, with print_raw=TRUE the attributes "levels", "partition_boundaries" and
    "n_classes" defining the diagram will be printed out.
    }

\item{...}{specifies further parameters that can be passed on to the print function.}
}
\description{
This is a function that prints out information on a Czekanowski's Diagram.
}

\author{
Krzysztof Bartoszek <krzbar@protonmail.ch>
}


\examples{
czkm<-czek_matrix(mtcars)
#'# Standard print ############
print(czkm)
print.czek_matrix(czkm)
# Print out the raw object ############
print(czkm,print_raw=TRUE)
print.czek_matrix(czkm,print_raw=TRUE)
}

\keyword{print}
