\name{read_maczek_file}
\alias{read_maczek_file}
\title{Read in MaCzek data.}
\usage{
read_maczek_file( filepath )
}
\arguments{
\item{filepath}{ the mdt MaCzek file to be read in.}
}
\value{
The function returns a data frame with the data.
}
\description{
The function reads in an mdt file that contains the data prepared for the MaCzek program
(MaCzek 3.3 - \url{http://www.antropologia.uw.edu.pl/MaCzek/maczek.html}).
The example data file (slabosz.mdt), concerning skull measurements from Slaboszewo,
comes from Piontek (1981).
}

\author{
Piotr Jaskulski

Maintainer: Krzysztof Bartoszek <krzbar@protonmail.ch>
}


\references{
J. Piontek J (1981).
Biologiczna charakterystyka sredniowiecznej populacji ze Slaboszewa, woj. bydgoskie.
Zrodla do badan biologii i historii populacji slowianskich, pages 39-83, Poznan UAM.
(in Polish, Biological characteristics of the medieval population of Slaboszewo,
Bydgoszcz Voivodeship. Materials for studying the biology and history of
Slavic populations.)
}



\examples{
filepath <- system.file("extdata", "slabosz.mdt", package="RMaCzek")
slabosz <- read_maczek_file(filepath)
res <- czek_matrix(slabosz)
}
