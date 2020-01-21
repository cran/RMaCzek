\name{internet_availability}
\alias{internet_availability}
\docType{data}
\title{Distances between Internet availability in some counties of the Silesia Voivodeship, Poland}
\description{
This is a dataset that describes Internet availability in 36 counties of the Silesia Voivodeship, Poland
to school children. It is derived from counts the number of students per comptuer with Internet
access at various levels of school education (Warzecha 2015).}
\usage{skulls_distances}
\format{
  The format is a list with three fields. The first field, \code{internet_availability_distances},
  is a 36 x 36 matrix. The second field, \code{full_county_names}, is a vector of length 36.
  The third field, \code{MaCzek_order_Warzecha}, is a vector of length 36.
}
\details{
A number of row/column names in the distances matrix are abbreviated in order to fit into 
the diagram. The full names of the counties is in the field \code{full_county_names}.
The field \code{MaCzek_order_Warzecha} contains the ordering found by Warzecha (2015)
using the MaCzek program, for comparison purposes.
}
\source{
The data is obtained personally from Katarzyna Warzecha and is included in the
package with her permission.
}
\references{
Warzecha, K. (2015);
The use of quantitative methods in research on selected behavioral addictions of young people. 
In: Studia Ekonomiczne, 247, 121--139.
}
\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")

internet_availability
internet_availability_distances<-as.dist(internet_availability$internet_availability_distances)

##  Warzecha (2015)'s original ordering by MaCzek
czek_matrix_internet_availability_sym<-czek_matrix(internet_availability_distances,
    original_diagram=FALSE,order=internet_availability$MaCzek_order_Warzecha)
plot(czek_matrix_internet_availability_sym,plot_title="",label.cex=0.5)

## Try to find a better ordering using the OLO method.
## We find the same clusters of counties as MaCzek did.
czek_matrix_internet_availability_OLO_sym<-czek_matrix(internet_availability_distances,
    original_diagram=FALSE,order="OLO")
plot(czek_matrix_internet_availability_OLO_sym,plot_title="",label.cex=0.5)
## Try to find a better ordering using the QAP_2SUM method.
## Even though Um is lower, the ordering does not seem to capture the
## clusters, suggesting that Hamiltonian path length minimization
## is a better strategy than Um minimization.
czek_matrix_internet_availability_qap2sum_sym<-czek_matrix(internet_availability_distances,
    original_diagram=FALSE,order="QAP_2SUM")
plot(czek_matrix_internet_availability_qap2sum_sym,plot_title="",label.cex=0.5)

## We try to see if anything more can be found from Czekanowski's original asymmetric
## diagram, but it seems that the plots are much less informative.
## We first try to make the grouping of similar objects in each column more
## related to the found clusters through the symmetric diagram.
column_order_stat_grouping<-c(8,10,12,16)
czek_matrix_internet_availability<-czek_matrix(internet_availability_distances,
    original_diagram=TRUE,order=internet_availability$MaCzek_order_Warzecha,
    column_order_stat_grouping=column_order_stat_grouping)
plot(czek_matrix_internet_availability,plot_title="",label.cex=0.5)
czek_matrix_internet_availability_OLO<-czek_matrix(internet_availability_distances,
    original_diagram=TRUE,order=attr(czek_matrix_internet_availability_OLO_sym,"order"),
    column_order_stat_grouping=column_order_stat_grouping)
plot(czek_matrix_internet_availability_OLO,plot_title="",label.cex=0.5)
czek_matrix_internet_availability_qap2sum<-czek_matrix(internet_availability_distances,
    original_diagram=TRUE,order=attr(czek_matrix_internet_availability_qap2sum_sym,"order"),
    column_order_stat_grouping=column_order_stat_grouping)
plot(czek_matrix_internet_availability_qap2sum,plot_title="",label.cex=0.5)

RNGversion(as.character(getRversion()))
}
\keyword{datasets}
