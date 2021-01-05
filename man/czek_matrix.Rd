\name{czek_matrix}
\alias{czek_matrix}
\title{Preprocess data to produce a Czekanowski's Diagram.}
\usage{
czek_matrix(x, order = "OLO", n_classes = 5, interval_breaks = NULL,
  monitor = FALSE, distfun = dist, scale_data = TRUE, focal_obj=NULL,
  as_dist=FALSE, original_diagram=FALSE, column_order_stat_grouping=NULL, 
  dist_args=list(), ...)
}
\arguments{
\item{x}{a numeric matrix, data frame or a 'dist' object.}

\item{order}{specifies which seriation method should be applied. The standard setting is the seriation method OLO. If NA or NULL, then no seriation is done and the original ordering is saved. The user may provide their own ordering, through a number vector of indices. Also in this case no rearrangement will be done.}

\item{n_classes}{specifies how many classes the distances should be divided into. The standard setting is 5 classes.}

\item{interval_breaks}{specifies the partition boundaries for the distances. As a standard setting, each class represents an equal amount of distances.  If the interval, breaks are positive and sum up to 1, then it is assumed that they specify percentages of the distances in each interval. Otherwise if provided as a numeric vector not summing up to 1, they specify the exact boundaries for the symbols representing distance groups.}

\item{monitor}{specifies if the distribution of the distances should be visualized. The standard setting is that the distribution will not be visualized. TRUE and "cumulativ_plot" is available.}

\item{distfun}{specifies which distance function should be used. Standard setting is the dist function which uses the Euclidean distance.  The first argument of the function has to be the matrix or data frame containing the data.}

\item{scale_data}{specifies if the data set should be scaled. The standard setting is that the data will be scaled.}

\item{focal_obj}{numbers or names of objects (rows if x is a dataset and not 'dist' object) that are not to take part in the reordering procedure. These observations will be placed as last rows and columns of the output matrix. See Details.}

\item{as_dist}{if TRUE, then the distance matrix of x is returned, with object ordering, instead of the matrix with the levels assigned in place of the original distances. The output will, be of class czek_matrix_dist, if FALSE, then of class czek_matrix.}

\item{original_diagram}{if TRUE, then the returned matrix corresponds as close as possible to the original method proposed by Czekanowski (1909). The levels are column specific and not matrix specific. See Details.}

\item{column_order_stat_grouping}{if original_diagram is TRUE, then here one can pass the partition boundaries for the ranking in each column.}

\item{dist_args}{specifies further parameters that can be passed on to the distance function.}

\item{...}{specifies further parameters that can be passed on to the seriate function in the seriation package.}
}

\details{
In his original paper Czekanowski (1909) did not have as the output a symmetric matrix
where each distance was assigned a level (symbol) depending in which numeric interval it was in.
Instead having the desired ordering, the following procedure was applied to each column.
The three smallest distances (in each column) obtain level (symbol) 1, the fourth smallest
level (symbol) 2, fifth smallest level (symbol) 3, sixth smallest level (symbol) 4 and
all the bigger distances the fifth symbol which was originally just a blank cell in the output
matrix. Here, we give the user more flexibility. In column_order_stat_grouping one may specify
how the order statistics should be grouped in each column. See Example.
The user may also choose some observations not to influence the ordering procedure.
This could be useful if e.g. a single observation is meant to be assigned to a cluster
and for some reason the clusters (that are to be read of from the ordering) should 
not be influenced by this observation. One can pass such observations using the
focal_obs parameter.

A hopefully useful property is that the ordering inside the czek_matrix (and hence
of the diagram when one calls plot) can be manually changed. One merely manipulates
the order attribute as desired. However in such a case one should remember that
the Path_length, criterion_value and Um attributes will have incorrect values
and should be corrected (see Examples).
}
\value{
The function returns a matrix with class czek_matrix. The return from the function is expected 
to be passed to the plot function. If as_dist is passed as TRUE, 
then a czek_matrix_dist object is returned and this is not suitable for the plotting.  
As an attribute of the output the optimized criterion value is returned. However, this is a 
guess based on seriation::seriate()'s and seriation::criterion()'s manuals. 
If something else was optimized, e.g. due to user's parameters, then this will be wrong. 
If unable to guess, then NA saved in the attribute.
}
\description{
This is a function that divided the values inside a distance matrix into classes. The output can be used in the plot function to produce a Czekanowski's Diagram.
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


# Create matrix as originally described by Czekanowski (1909), with each column
# assigned levels according to how the order statistics of the  distances in it
# are grouped. The grouping below is the one used by Czekanowski (1909).
czek_matrix(x,original_diagram=TRUE,column_order_stat_grouping=c(3,4,5,6))

# Create matrix with two focal object that will not influence seriation
czek_matrix(x,focal_obj=c("Merc 280","Merc 450SL"))
# Same results but with object indices
czek_res<-czek_matrix(x,focal_obj=c(10,13))

## we now place the two objects in a new place
czek_res_neworder<-manual_reorder(czek_res,c(1:10,31,11:20,32,21:30),
    orig_data=x)

## the same can be alternatively done by hand
attr(czek_res,"order")<-attr(czek_res,"order")[c(1:10,31,11:20,32,21:30)]
## and then correct the values of the different criteria so that they
## are consistent with the new ordering
attr(czek_res,"Path_length")<-seriation::criterion(stats::dist(scale(x)),
    order=seriation::ser_permutation(attr(czek_res, "order")),method="Path_length")
## Here we need to know what criterion was used for the seriation procedure
## If the seriation package was used, then see the manual for seriation::seriate()
## seriation::criterion().
## If the genetic algorithm shipped with RMaCzek was used, then it was the Um factor.
attr(czek_res,"criterion_value")<-seriation::criterion(stats::dist(scale(x)),
    order=seriation::ser_permutation(attr(czek_res, "order")),method="Path_length")
attr(czek_res,"Um")<-RMaCzek::Um_factor(stats::dist(scale(x)),
    order= attr(czek_res, "order"),inverse_um=FALSE) 
}
}
\keyword{hplot}


