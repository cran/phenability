\name{kang}
\alias{kang}
\title{
Kang Function
}
\description{
This function performs a stability analysis based on the Kang (1988) criteria. Kang nonparametric stability (ranksum) uses both "trait single value" and stability variance (Shukla, 1972), and the genotype with the lowest ranksum is commonly the most favorable one.
}
\usage{
kang(data)
}
\arguments{
  \item{kang}{kang function}
  \item{data}{A numeric data.frame with fenotypic means values of MET, and genotypes and       environments by rows and columns, respectively.}
}


\value{
Is returned a numerical "ranksum" value by genotype.
}
\references{
Kang, M.S. 1988. A ranksum method for selecting high yielding,
stable corn genotypes. Cereal Res. Comm. 16:113.\cr

Shukla, G.K. 1972. Some aspects of partitioning genotype environmental
components of variability. Heredity 28:237.
}
\author{
Leonardo Castelo Branco
}


\examples{
require("phenability")
data(bean.data)
bean.data
x <- bean.data

kang(x)

}
\keyword{ kang }
\keyword{ phenability }
