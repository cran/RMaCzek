\name{urns}
\alias{urns}
\docType{data}
\title{Measurements of urns from cremational graves}
\description{
Measurements of urns from cremational graves excavated at Paprotki
Kolonia 12 in Poland (Table p. 181, Soltysiak and Jaskulski 1999).
}
\usage{urns}
\format{
  The format is a 15 x 9 matrix.
}
\details{
Urn "gr.52.1" (row 9) has too many missing values and probably should 
be removed for further analyses.

The column names correspond to height (WYS), rim diameter (SW), maximal diameter (MWB), 
bottom diameter (SD), average wall thickness (GS), average bottom thickness (GD) and 
three indices describing proportions of the vessel (W-A, W-B, W-D).
}
\source{
The matrix is taken from Soltysiak and Jaskulski (1999).
}
\references{
Soltysiak, A. and Jaskulski, P. (1999);
"Czekanowski's Diagram. A Method of Multidimensional Clustering." 
New Techniques for Old Times. CAA 98. Computer Applications and Quantitative Methods in 
Archaeology. Proceedings of the 26th Conference, Barcelona, March 1998, 
number 757 in BAR International Series, pp. 175-184. Oxford.
}
\examples{
urns
urns_use<-urns[-9,] ## removed as too many missing values, observation "gr.52.1"

## proposed order by Soltysiak and Jaskulski (1999) from the MaCzek program
urns_use_proposedorder<-c(1,3,5,7,2,4,8,11,10,12,6,9,14,13) 
czkm_urns<-czek_matrix(urns_use_proposedorder,original_diagram=FALSE,
    order=urns_use_proposedorder,scale_data=TRUE)
plot(czkm_urns,plot_title="",label.cex=0.9)

## try to find a better ordering
czkm_urns_OLO<-czek_matrix(urns_use,original_diagram=FALSE,order="OLO",scale_data=TRUE)
plot(czkm_urns_OLO,plot_title="",label.cex=0.9)

## Construct Czekanowski's original non-symmetric diagram
czkm_urns_OLO_2<-czek_matrix(urns_use,original_diagram=TRUE,order="OLO",scale_data=TRUE)
plot(czkm_urns_OLO_2,plot_title="",label.cex=0.9)
}
\keyword{datasets}
