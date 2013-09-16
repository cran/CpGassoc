\name{design}
\alias{design}
\alias{design.data.frame}
\alias{design.NULL}

\title{
Create full and reduced design matrices for the cpg.assoc function.
}
\description{
Designed to be used by \code{cpg.assoc} and \code{cpg.perm}. Creates a full and reduced design matrices.
}
\usage{
design(covariates, indep, chip.id, random)
}
\arguments{
  \item{covariates}{
A data frame consisting of the covariates of interest. covariates can also be a matrix if it is a model matrix minus the intercept column. 
It can also be a vector if there is only one covariate of interest.
If no covariates must be specified as \code{NULL}.
}
  \item{indep}{
A vector containing the main variable of interest.  \code{cpg.assoc} will evaluate the association between indep and the beta values.
}
  \item{chip.id}{
An optional vector containing chip or batch identities. If specified, \code{chip.id} will be included as a factor in the model.
}
  \item{random}{
Is the model going to be a mixed effects. If so, \code{chip.id} will not be included in the design matrices.
}
}

\value{
Returns a list containing the full and reduced design matrices.
  \item{full }{The full design matrix.}
  \item{reduced }{The reduced design matrix.}

}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}
\note{
The \code{design} function  is designed to be used exclusively by the \code{cpg.assoc} and \code{cpg.perm} functions.
}



\seealso{
\code{\link{cpg.assoc}}
\code{\link{cpg.perm}}
\code{\link{cpg.work}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.combine}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
}
\examples{
data(samplecpg,samplepheno,package="CpGassoc")
#Example where there are covariates:
covar<-data.frame(samplepheno$weight,samplepheno$Distance)
test<-design(covar,samplepheno$SBP,samplepheno$chip,FALSE)
dim(test$full)
dim(test$reduced)
test$reduced[1:5,1:5]
test$full[1:5,1:5]
#When no covariates or chip.id:
test2<-design(NULL,samplepheno$SBP,NULL,FALSE)
dim(test2$full)
dim(test2$reduced)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
