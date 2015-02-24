\name{phenability-package}
\alias{phenability-package}
\alias{phenability}
\docType{package}
\title{Nonparametric Stability Analysis}
\description{
An alternative to carrying out phenotypic adaptability and stability analysis, taking into account nonparametric statistics. Can be used as a robust approach, and less sensitive to departures from common genotypic, environmental, and GxE effects data assumptions (e.g., normal distribution of errors).
}
\details{
\tabular{ll}{
Package: \tab phenability\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-02-24\cr
License: \tab GPU-3\cr
}
The package can be applied to MET (multienvironment trials) data containing the phenotypic mean of genotype x environment (dataframe with genotypes by rows and environments into columns) and performs nonparametric analysis according to for criteria: Huehn (1979), and Nassar and Huehn (1987); Kang's (1988); Fox et al. (1990) and; Thennarasu's (1995).
}

\arguments{
  \item{function}{ nahu, kang, fox, or thsu}
  \item{data}{ A numeric data.frame with fenotypic means values of MET, and genotypes and environments by rows and columns, respectively.}
  \item{interaction}{do a interaction plot or not (FALSE = no interaction plot output, TRUE = draws a interaction plot between genotypes and environments).}
}

\author{
Leonardo Castelo Branco

Maintainer: Leonardo C. Branco <leonardocb@usp.br>
}
\examples{
require(phenability)
data(bean.data)
bean.data
x <- bean.data

fox(x)
fox(x, interaction=TRUE)

kang(x)
kang(x, interaction=TRUE)

thsu(x)
thsu(x, interaction=TRUE)

nahu(x)
nahu(x, interaction=TRUE)

}

\references{
Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R. Cormier.
1990. Yield and adaptation of hexaploid spring triticale. Euphytica
47:57.\cr

Huehn, V.M. 1979. Beitrage zur erfassung der phanotypischen
stabilitat. EDV Med. Biol. 10:112.\cr

Kang, M.S. 1988. A ranksum method for selecting high yielding,
stable corn genotypes. Cereal Res. Comm. 16:113.\cr

Lin, C.S., and M.R. Binns. 1994. Concepts and methods for analyzing
regional trial data for cultivar and location selection. Plant Breed.
Rev. 12:271.\cr

Lin, C.S., M.R. Binns, and L.P. Lefkovitch. 1986. Stability analysis:
Where do we stand? Crop Sci. 26:894.\cr

Nassar, R., and M. Huehn. 1987. Studies on estimation of phenotypic
stability: Tests of significance for nonparametric measures of
phenotypic stability. Biometrics 43:45.\cr

Shukla, G.K. 1972. Some aspects of partitioning genotype environmental
components of variability. Heredity 28:237.\cr

Thennarasu, K. 1995. On certain nonparametric procedures for
studying genotype x environment interactions and yield stability.
Ph.D. thesis. P.J. School, IARI, New Delhi, India.


}


\keyword{phenability}
