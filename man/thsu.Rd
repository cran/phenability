\name{thsu}
\alias{thsu}
\title{
Thennarasu Function
}
\description{
Performs a stability analysis according to Thennarasu (1995) methodology. Therefore the statistics N1, N2, N3, N4 are presented in output. 
}
\usage{
thsu(data, interaction=FALSE)
}
\arguments{
  \item{thsu}{thsu function}
  \item{data}{ A numeric data.frame with fenotypic means values of MET, and genotypes and       environments by rows and columns, respectively.}
  \item{interaction}{do a interaction plot or not (FALSE = no interaction plot output, TRUE = draws a interaction plot between genotypes and environments).}
}


\value{
Is returned a list with four numerical statistics values (N1, N2, N3, and N4) by genotype, the ranking for the response variable in each environment, and the Spearman correlation matrix between environments.
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

thsu(x, interaction = TRUE)
}
\keyword{ ~Thennarasu }
\keyword{ ~phenability }
