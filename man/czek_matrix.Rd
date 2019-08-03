\name{czek_matrix}
\alias{czek_matrix}
\title{Preprocess data to produce a Czekanowski's Diagram.}
\usage{
czek_matrix(x, order = "OLO", n_classes = 5, interval_breaks = NULL,
  monitor = FALSE, distfun = dist, scale_data = TRUE, ...)
}
\arguments{
\item{x}{a numeric matrix, data frame or a 'dist' object.}

\item{order}{specifies which seriation method should be applied. The standard setting is the seriation method OLO.}

\item{n_classes}{specifies how many classes the distances should be divided into. The standard setting is 5 classes.}

\item{interval_breaks}{specifies the partition boundaries for the distances. As a standard setting, each class represents an equal amount of distances.}

\item{monitor}{specifies if the distribution of the distances should be visualized. The standard setting is that the distribution will not be visualized. TRUE and "cumulativ_plot" is available.}

\item{distfun}{specifies which distance function should be used. Standard setting is the dist function which uses the Euclidean distance.}

\item{scale_data}{specifies if the data set should be scaled. The standard setting is that the data will be scaled.}

\item{...}{specifies further parameters that can be passed on to the seriate function in the seriation package.}
}
\value{
The function returns a matrix with class czek_matrix. The return from the function is expected to be passed to the plot function.
}
\description{
This is a function that divided the values inside a distance matrix into classes. The output can be used in the plot function to produce a Czekanowski's Diagram.
}


\author{
Albin Vasterlund

Maintainer: Krzysztof Bartoszek <krzbar@protonmail.ch>
}

\references{
J. Czekanowski. (1909).
Zur Differentialdiagnose der Neandertalgruppe.
Korespondentblatt der Deutschen Gesellschaft fur Anthropologie, 
Ethnologie und Urgeschichte, XL(6/7):44-47,

A Soltysiak and P. Jaskulski (1999). 
Czekanowski's diagram. a method of multidimensional clustering. 
New Techniques for Old Times. CAA 98. Computer Applications and
Quantitative Methods in Archaeology. Proceedings of the 26th Conference, Barcelona,
March 1998, number 757 in BAR International Series, pages 175-184, Oxford.
 
Vasterlund, A. (2019). 
Czekanowski's Diagram: Implementing and exploring Czekanowski's Diagram with different seriation methods
Master thesis, Linkoping University
}

\examples{
# Set data ####
x<-mtcars


# Different type of input that give same result ############
czek_matrix(x)
czek_matrix(stats::dist(scale(x)))

## below a number of other options are shown
## but they take too long to run
\donttest{
# Change seriation method ############
#seriation::show_seriation_methods("dist")
czek_matrix(x,order = "GW")
czek_matrix(x,order = "ga")
czek_matrix(x,order = sample(1:nrow(x)))


# Change number of classes ############
czek_matrix(x,n_classes = 3)


# Change the partition boundaries ############
czek_matrix(x,interval_breaks = c(0.1,0.4,0.5)) #10\%, 40\% and 50\%
czek_matrix(x,interval_breaks = c(0,1,4,6,8.48)) #[0,1] (1,4] (4,6] (6,8.48]
czek_matrix(x,interval_breaks = "equal_width_between_classes") 
#[0,1.7] (1.7,3.39]  (3.39,5.09] (5.09,6.78] (6.78,8.48]


# Change number of classes ############
czek_matrix(x,monitor = TRUE)
czek_matrix(x,monitor = "cumulativ_plot")


# Change distance function ############
czek_matrix(x,distfun = function(x) stats::dist(x,method = "manhattan"))


# Change dont scale the data ############
czek_matrix(x,scale_data = FALSE)
czek_matrix(stats::dist(x))


# Change additional settings to the seriation method ############
czek_matrix(x,order="ga",control=list(popSize=200,
                                     suggestions=c("SPIN_STS","QAP_2SUM")))
}
}
\keyword{hplot}

