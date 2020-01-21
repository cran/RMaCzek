\name{seals_similarities}
\alias{seals_similarities}
\docType{data}
\title{Similarities between Akkadian cylinder seals depicting the Serpent God}
\description{
A data martix containing the similarities between 37 seals
depicting the Serpent God (Table 10, Soltysiak 2000).
}
\usage{skulls_distances}
\format{
  The format is a 37 x 37 matrix.
}
\details{
In the original similarities table, there are no values on the diagonal 
(i.e. self-similarity). Hence for furtherer analyses, something has to
be filled in, e.g. 100.
}
\source{
The similarity matrix is taken from Soltysiak (2000).
}
\references{
Soltysiak, A. (2000);
"Przedstawienie Boga-Weza na pieczeciach cylindrycznych z okresu
akadyjskiego: analiza ikonograficzna." Studia i Materialy Archeologiczne, 10, 189-214.
(Depiction of the Serpent God on cylinder seals from the Akkadian era: and iconogrphic
analysis, in Polish).
}
\examples{
seals_similarities
## there are only NAs on the diagonal (self-similarity)
diag(seals_similarities)<-100
## need to change it to distance matrix for analysis by czek_matrix
seals_distances<-as.dist(100-seals_similarities)
int_breaks<-c(0,40,60,80,100)

## Soltysiak (2000)'s original ordering
czkm_seals<-czek_matrix(seals_distances,original_diagram=FALSE,order=NA,
    n_classes=4,interval_breaks=int_breaks)
plot(czkm_seals,plot_title="",label.cex=0.7)

## try to find a better ordering
czkm_seals_OLO<-czek_matrix(seals_distances,original_diagram=FALSE,order="OLO",
    n_classes=4,interval_breaks=int_breaks)    
plot(czkm_seals_OLO,plot_title="",label.cex=0.7)

## Construct Czekanowski's original non-symmetric diagram
czkm_seals_OLO_2<-czek_matrix(seals_distances,original_diagram=TRUE,order="OLO")
plot(czkm_seals_OLO_2,plot_title="",label.cex=0.7)
}
\keyword{datasets}
