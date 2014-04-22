\name{thsu}
\alias{thsu}
\title{
Thennarasu Function
}
\description{
This function performs a stability analysis according to Thennarasu (1995) methodology. therefore the statistics N1, N2, N3, N4 are presented in output. 
}
\usage{
thsu(data)
}
\arguments{
  \item{thsu}{thsu function}
  \item{data}{ A numeric data.frame with fenotypic means values of MET, and genotypes and       environments by rows and columns, respectively.}
}


\value{
Is returned four numerical statistics values (N1, N2, N3, and N4) by genotype.
}
\references{
Thennarasu, K. 1995. On certain nonparametric procedures for
studying genotype x environment interactions and yield stability.
Ph.D. thesis. P.J. School, IARI, New Delhi, India.
}
\author{
Leonardo Castelo Branco
}

\examples{
require("phenability")
data(bean.data)
bean.data
x <- bean.data

thsu(x)
}
\keyword{ ~Thennarasu }
\keyword{ ~phenability }
