\name{cpg.perm}
\alias{cpg.perm}

\title{
Perform a Permutation Test of the Association Between Methylation and a Phenotype of Interest
}
\description{
Calls \code{cpg.assoc} to get the observed P-values from the study and then performs a user-specified number of permutations to calculate an emperical p-value.
In addition to the same test statistics computed by \code{cpg.assoc}, \code{cpg.perm} will compute the permutation p-values for the observed p-value, the number of Holm significant sites, and the number of FDR significant sites.
}
\usage{
cpg.perm(beta.values, indep, covariates = NULL, nperm, data = NULL, seed = NULL, logit.transform = FALSE, chip.id = NULL, subset = NULL, random = FALSE, fdr.cutoff = 0.05, fdr.method = "BH",large.data=TRUE)
}

\arguments{
  \item{beta.values}{
A vector, matrix, or data frame containing the beta values of interest (1 row per CpG site, 1 column per individual).
}
  \item{indep}{
A vector containing the main variable of interest.  \code{cpg.assoc} will evaluate the association between indep and the beta values. 
}
  \item{covariates}{
A data frame consisting of the covariates of interest. covariates can also be a matrix if it is a model matrix minus the intercept column. 
It can also be a vector if there is only one covariate of interest.  Can also be a formula(e.g. \code{~cov1+cov2}).
}
  \item{nperm}{
The number of permutations to be performed.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from the environment from which \code{cpg.perm} is called.
}
  \item{seed}{
The required seed for random number generation. If not input, will use R's internal seed.
}
  \item{logit.transform}{
logical. If \code{TRUE}, the logit transform of the beta values log(beta.val/(1-beta.val)) will be used. Any values equal to zero or one will be set to the 
next smallest or  next largest value respectively; values <0 or >1 will be set to NA.
}
  \item{chip.id}{
An optional vector containing the chip information. If specified, chip id will be included as a factor in the model.
}
  \item{subset}{
An optional logical vector specifying a subset of observations to be used in the fitting process.
}
  \item{random}{
logical. If \code{TRUE}, the \option{chip.id} will be processed as a random effect, and a random intercept model will be fitted.
}
  \item{fdr.cutoff}{
The threshold at which to compare the FDR values. The default setting is .05. Any FDR values less than .05 will be considered significant.
}
  \item{fdr.method}{
Character. Method used to calculate False Discovery Rate. Can be any of the methods listed in \code{\link{p.adjust}} or "qvalue" 
for John Storey's qvalue method (required to have \pkg{qvalue} package installed). The default method is "BH" for the Benjamini & Hochberg method.
}
  \item{large.data}{
Logical.  Enables analyses of large datasets.  When \code{large.data=TRUE}, \code{cpg.assoc} avoids memory problems by performing the analysis in chunks. 
}
}

\value{
The item returned will be of class \code{"cpg.perm"}. It will contain all of the values of class cpg (\code{\link{cpg.assoc}}) and a few more:

\item{permutation.matrix}{
A matrix consisting of the minimum observed P-value, the number of Holm significant CpG sites, and the number of FDR significant sites for each permutation.
}
\item{gc.permutation.matrix}{
Similar to the permutation.matrix only in relation to the genomic control adjusted p-values.
}
\item{perm.p.values}{
A data frame consisting of the permutation P-values, and the number of permutations performed.
  }
\item{perm.tstat}{
If one hundred or more permutations were performed and indep is a continuous variable, consists of the quantile .025 and .975 of observed t-statistcs for each permutation, ordered from smallest to largest.
perm.tstat is used by \code{plot.cpg.perm} to compute the confidence intervals for the QQ plot of t-statistics. Otherwise \code{NULL}.
}
\item{perm.pval}{
If one hundred or more permutations were performed, consists of the observed p-values for each permutation, ordered from smallest to largest. perm.pval is usd by \code{plot.cpg.perm} to compute
the confidence intervals for the QQ plot of the p-values. Otherwise NULL.
}
}

\author{
Barfield, R.; Conneely, K.; Kilaru,V.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}


\seealso{
\code{\link{cpg.assoc}}
\code{\link{cpg.work}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.combine}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
\code{\link{sort.cpg.perm}}
\code{\link{sort.cpg}}
\code{\link{cpg.qc}}
}
\examples{
##Loading the data
data(samplecpg,samplepheno,package="CpGassoc")
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
#Performing a permutation 10 times
Testperm<-cpg.perm(samplecpg[1:200,],samplepheno$weight,seed=2314,nperm=10,large.data=FALSE)
Testperm
#All the contents of CpGassoc are included in the output from Testperm

#summary function works on objects of class cpg.perm
summary(Testperm)

}
