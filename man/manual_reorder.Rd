\name{manual_reorder}
\alias{manual_reorder}
\title{Manually reorder Czekanowski's Diagram}
\usage{
manual_reorder(x, v_neworder, ...)
}
\arguments{
\item{x}{ a matrix with class czek_matrix, czek_matrix_dist or data matrix/data.frame or dist object.}

\item{v_neworder}{a numeric vector with the new ordering.}

\item{...}{specifies further parameters that can be passed on to the czek_matrix function or will be present in derived functions. See details and example code.}
}
\description{
The function returns a Czekanowski's Diagram with the new order and recalculated factors
}

\details{
This is an S3 generic function. The RMaCzek package defines five methods for
it:
\itemize{
\item{manual_reorder.czek_matrix}{In this case the input is of czek_matrix class and the user has to specify a further parameter orig_data. This is as the czek_matrix object does not hold information on the distances and hence they will need to be recalculated. The user should pass all the parameters that were passed to the original call to czek_matrix(), except order and as_dist.}
\item{manual_reorder.czek_matrix_dist}{In this case the input is of czek_matrix_dist class and nothing additionally needs to be provided.}
\item{manual_reorder.data.frame}{In this case the input is a data.frame.  The user should pass all the parameters that they would have passed to the original call to czek_matrix(), except order.}
\item{manual_reorder.matrix}{In this case the input is a matrix of measurements (not distance matrix).  The user should pass all the parameters that they would have passed to the original call to czek_matrix(), except order.}
\item{manual_reorder.dist}{In this case the input is a dist object (e.g. distance matrix after calling as.dist()).  The user should pass all the parameters that they would have passed to the original call to czek_matrix(), except order.}
}
}

\author{
Krzysztof Bartoszek <krzbar@protonmail.ch>
}


\examples{
#'# Set data ####
x<-mtcars
# Calculate Czekanowski's diagram
czkm<-czek_matrix(x)
czkm_dist<-czek_matrix(x,as_dist=TRUE)
# new ordering
neworder<-attr(czkm,"order")
neworder[1:2]<-neworder[2:1]
# reorder the diagram
#if the output was Czekanowski's diagram without the distances
#remembered, then the original data has to be passed so that 
#factors can be recalculated.
new_czkm<-manual_reorder(czkm,v_neworder=neworder,orig_data=x)
new_czkm_dist<-manual_reorder(czkm_dist,v_neworder=neworder)
#we can also pass the original data directly
new_czkm<-manual_reorder(x,v_neworder=neworder)
#and this is equivalent to calling
czkm<-czek_matrix(x,order=neworder)
#up to the value of the "criterion_value" attribute
#which in the second case can be lost, as no information is passed
#on which one was originally used, while in the first case it might
#be impossible to recalculate-only criteria values from seriate are supported
#if a user has a custom seriation function, then they need to recalculate this
#value themselves
}

\keyword{dplot}
