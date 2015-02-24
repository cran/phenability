\name{fox}
\alias{fox}
\title{
Fox Function
}
\description{
Performs a stability analysis based on the criteria of Fox et al. (1990), using the statistical "TOP third" only. In Fox function, a stratified ranking of the genotypes at each environment separately is done. The proportion of sites at which the genotype occurred in the top third are expressed in TOP output.
}
\usage{
fox(data, interaction=FALSE)
}
\arguments{
  \item{fox}{fox function}
  \item{data}{A numeric data.frame with fenotypic means values of MET, and genotypes and       environments by rows and columns, respectively.}
  \item{interaction}{do a interaction plot or not (FALSE = no interaction plot output, TRUE = draws a interaction plot between genotypes and environments).}
}


\value{
Is returned a list with numerical TOP ranking value by genotype, the ranking for the response variable in each environment, and the Spearman correlation matrix between environments.
}
\references{
Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R. Cormier.
1990. Yield and adaptation of hexaploid spring triticale. Euphytica
47:57.
}
\author{
Leonardo Castelo Branco
}

\examples{
require("phenability")
data(bean.data)
bean.data
x <- bean.data

fox(x, interaction=TRUE)

}

\keyword{ Fox }
\keyword{ phenability }
