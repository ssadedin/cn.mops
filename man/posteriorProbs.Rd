\name{posteriorProbs}
\alias{posteriorProbs}
\title{This generic function returns the posterior probabilities of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.
The posterior probabilities are represented
as a three dimensional array, where the three dimensions are segment,
copy number and individual.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{posteriorProbs} returns a three dimensional array.
}
\description{
  This generic function returns the posterior probabilities
  of a CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}. The posterior
  probabilities are represented as a three dimensional
  array, where the three dimensions are segment, copy
  number and individual.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
posteriorProbs(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

