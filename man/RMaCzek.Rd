\docType{package}
\name{RMaCzek}
\alias{RMaCzek}
\alias{RMaCzek-package}
\title{RMaCzek: A package that can produce Czekanowski's diagram}
\description{
This package produces Czekanowski's diagram.

This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
Please understand that there may still be bugs and errors. Use it at your own risk. 
We take no responsibility for any errors or omissions in this package or for any misfortune 
that may befall you or others as a result of its use. Please send comments and report 
bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
}
\details{
\tabular{ll}{
Package: \tab RMaCzek\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2022-03-18\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
This package produces Czekanowski's diagram.
}
\section{The packages functions}{

\describe{
  \item{czek_matrix}{A function that returns a distance matrix where the distances are divided into classes. The return from the function is expected to be passed into the plot function. }
  \item{plot.czek_matrix}{A function that returns Czekanowski's Diagram.}
}
}

\author{
Albin Vasterlund

Piotr Jaskulski

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
Ethnologie und Urgeschichte, XL(6/7):44-47.

A. Soltysiak and P. Jaskulski (1999). 
Czekanowski's diagram. a method of multidimensional clustering. 
New Techniques for Old Times. CAA 98. Computer Applications and
Quantitative Methods in Archaeology. Proceedings of the 26th Conference, Barcelona,
March 1998, number 757 in BAR International Series, pages 175-184, Oxford.
 
A. Vasterlund (2019). 
Czekanowski's Diagram: Implementing and exploring Czekanowski's Diagram with different seriation methods
Master thesis, Linkoping University.
}

\keyword{hplot}

\examples{
# Set data ####
x<-czek_matrix(mtcars)


# Standard plot ############
plot(x)
plot.czek_matrix(x)
}
