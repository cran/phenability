\name{nahu}
\alias{nahu}
\title{
Nassar and Huehn Function
}
\description{
Performs a stability analysis based on Nassar and Huehn (1987) and, Huehn (1979) statistics. the four nonparametric measures of phenotypic stability are: S1 (mean of the absolute rank differences of a genotype over the n environments), S2 (variance among the ranks over the k environments), S3 (sum of the absolute deviations), and S6 (relative sum of squares of rank for each genotype).
}
\usage{
nahu(data, interaction=FALSE)
}

\arguments{
  \item{nahu}{nahu function}
  \item{data}{ A numeric data.frame with fenotypic means values of MET, and genotypes and environments by rows and columns, respectively.}
  \item{interaction}{do a interaction plot or not (FALSE = no interaction plot output, TRUE = draws a interaction plot between genotypes and environments).}
}


\value{
Is returned a list with four numerical statistics values (S1, S2, S3, and S6) by genotype, the ranking for the response variable in each environment, and the Spearman correlation matrix between environments.
}
\references{
Nassar, R., and M. Huehn. 1987. Studies on estimation of phenotypic
stability: Tests of significance for nonparametric measures of
phenotypic stability. Biometrics 43:45.\cr

Huehn, V.M. 1979. Beitrage zur erfassung der phanotypischen
stabilitat. EDV Med. Biol. 10:112.
}
\author{
Leonardo Castelo Branco
}

\examples{
require("phenability")
data(bean.data)
bean.data
x <- bean.data

nahu(x, interaction = TRUE)
}
\keyword{ huehn }
\keyword{ phenability }
