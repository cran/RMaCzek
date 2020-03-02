\name{skulls_distances}
\alias{skulls_distances}
\docType{data}
\title{Distances between archaic human skulls}
\description{
A data matrix containing the distances between 13 skulls
of archaic humans (Table II, Czekanowski, 1909).
}
\usage{skulls_distances}
\format{
  The format is a 13 x 13 matrix.
}
\details{
In the original paper (Table II, Czekanowski, 1909) there are minor 
typographic error. Firstly,
distance(Neandertal,Galey Hill)=10.54
while
distance(Galey Hill,Neandertal)=10.504 . 
However, if one analysis the source of the data (Stolyhwa 1908),
one will see that it should be 10.504 .
Furthermore, the third row is called "Krapina G", while it should be 
"Krapina C".
}
\source{
The distance matrix is taken from Czekanowski (1909), while the original
skull measurements can be found in Stolyhwa (1908).
}
\references{
Czekanowski, J. (1909); 
Zur Differentialdiagnose der Neandertalgruppe. 
In: Korespondentblatt der Deutschen Gesellschaft fur Anthropologie, Ethnologie und Urgeschichte, 
XL(6/7), 44--47.

Stolyhwa, K. (1908); 
Czaszka z Nowosiolki jako dowod istnienia w okresie historycznym ksztaltow pokrewnych z Homo primigenius.
In: Rozprawy Wydzialu matematyczno--przyrodnicznego Akademii Umiejetnosci, 
XLVIII(B), 1--27. (The skull from Nowosiolka as proof of existence during the era of history shapes common 
with Homo primigenius, in  Polish).
}
\examples{
## correcting typo rowname
rownames(skulls_distances)[3]<-"Krapina C"
## d(Neandertal,Galey Hill)!=d(Galey Hill,Neandertal)
## in the original paper so we need to correct, 
## e.g. as based on Stolyhwa (1908)'s data
skulls_distances
sym_skulls_distances<-skulls_distances
sym_skulls_distances[5,9]<-10.504
## To obtain Czekanowski (1909)'s original diagram (with different symbols)
czek_matrix_skulls<-czek_matrix(as.dist(sym_skulls_distances),order=NA,original_diagram=TRUE)
plot(czek_matrix_skulls,plot_title="",label.cex=0.5)
## Obtain a symmetric version of the diagram
czek_matrix_skulls_sym<-czek_matrix(as.dist(sym_skulls_distances),order=NA,original_diagram=FALSE)
plot(czek_matrix_skulls_sym,plot_title="",label.cex=0.5)
## Try to find a better ordering
czek_matrix_skulls_OLO<-czek_matrix(as.dist(sym_skulls_distances),order="OLO",original_diagram=TRUE)
plot(czek_matrix_skulls_OLO,plot_title="",label.cex=0.5)
## Better ordering with original symmetric diagram
czek_matrix_skulls_OLO_sym<-czek_matrix(as.dist(sym_skulls_distances),order="OLO",
    original_diagram=FALSE)
plot(czek_matrix_skulls_OLO_sym,plot_title="",label.cex=0.5)
}
\keyword{datasets}
