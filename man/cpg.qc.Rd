\name{cpg.qc}
\alias{cpg.qc}

\title{
Performs quality control on Illumina data.
}
\description{
\code{cpg.qc} is designed to perform quality control on Illumina data prior to analysis. 
In addition to the matrix of beta values, this function requires as input matrices of Signal A, Signal B, and detection p-values.
It will remove samples that have low intensity (mean signal intensity less than half of the overall median or 2000). 
It can also set to NA datapoints with detection p-values exceeding a user-specified cutoff, and can remove samples or 
sites that have a missing rate above a user-specified value. 
Finally, users can opt to compute beta values as M/(U+M) or M/(U+M+100).
}
\usage{
cpg.qc(beta.orig,siga,sigb,pval,p.cutoff=.001,cpg.miss=NULL,sample.miss=NULL,
constant100=FALSE,sig.return=FALSE)
}

\arguments{
   \item{beta.orig}{
The original beta values matrix obtained from GenomeStudio.
}
  \item{siga}{
The unmethylated signals matrix obtained from GenomeStudio.
}
  \item{sigb}{
The methylated signals matrix obtained from GenomeStudio.
}
  \item{pval}{
A matrix of detection p-values obtained from GenomeStudio. 
 pval should have the same dimension as the beta values and signals: one row for each site and one column for each individual.
}
  \item{p.cutoff}{
The user-specified cutoff for detection p-values (default=.001).
}
  \item{cpg.miss}{
Optional cutoff value. If specified, cpg.qc will remove cpg sites where the proportion of missing values exceeds this cutoff.
}
  \item{sample.miss}{
Optional cutoff value. If specified, cpg.qc will remove samples where the proportion of missing values exceeds this cutoff.
}
  \item{constant100}{
Logical. If true, the new beta values will be calculated as M/(U+M+100); if false (default) they will be calculated as M/(U+M).
}
  \item{sig.return}{
Logical. If true, cpg.qc returns the qced signal data as well.
}
}
\details{
It is important that all the matrices or data frames listed above (\option{pval}, \option{siga}, \option{sigb}, \option{beta.orig}) 
are ordered similarly with respect to samples and CpG sites.
}
\value{
\code{cpg.qc} returns a new matrix of beta values that has been subjected to the specified quality control filters.  This matrix can be input directly into \code{cpg.assoc}.
}

\author{
Barfield, R.; Conneely, K.; Kilaru,V.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}

\seealso{
\code{\link{cpg.work}}
\code{\link{cpg.perm}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.combine}}
\code{\link{cpg.assoc}}

}


