\name{Sample Data CpGassoc}
\alias{samplecpg}
\alias{samplepheno}
\alias{annotation}
\docType{data}
\title{
Sample data from \pkg{CpGassoc}
}
\description{
  \itemize{
  \item{samplecpg}{Matrix containing sample (fake) methylation data for 258 individuals over 1228 CpG sites}
  \item{samplepheno}{Matrix with phenotype info for the 258 individuals}
  \item{annotation}{Matrix with annotation information for the 1228 CpG sites}
  }
}
\usage{data(samplecpg)}
\format{
 \itemize{
  \item{samplecpg}{Matrix, dimensions 1228 x 258}
  \item{samplepheno}{Matrix, dimensions 258 x 6. Phenotype info on \code{"Dose"}, \code{"SBP"}, \code{"Distance"}, \code{"weight"}, \code{"Disease"}, and \code{"chip"}}
  \item{annotation}{Matrix, dimensions 1228 x 3. Header names are \code{"TargetID"}, \code{"CHR"}, and \code{"MAPINFO"}}
  }
}

\seealso{
\code{\link{cpg.perm}}
\code{\link{cpg.assoc}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}

}

\examples{
##See help pages for other functions for usage of these datasets
}
\keyword{datasets}
