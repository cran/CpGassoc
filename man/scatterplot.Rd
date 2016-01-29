\name{scatterplot}
\alias{scatterplot}

\title{
Plot beta values of individual CpG sites against the independent variable.
}
\description{
Plot beta values of individual CpG sites against the independent variable. Can create scatterplots and boxplots. If scatterplots the intercept will be adjusted
for any covariates that were included in the model.
}
\usage{
scatterplot(x, cpg.rank = NULL, cpg.name = NULL, save.plot = NULL, file.type="pdf", 
eps.size = c(5, 5), popup.pdf = FALSE, beta.values = NULL,main.title=NULL, ...)
}

\arguments{
  \item{x}{
Object of class \code{"cpg"} or \code{"cpg.perm"}.
}
  \item{cpg.rank}{
A vector listing the rank of sites to be plotted. The rank is based on the ordered p-values.
}
  \item{cpg.name}{
A character vector containing the names of CpG sites to be plotted against the phenotype of interest. This option is ignored if \option{cpg.rank} is specified.
}
  \item{save.plot}{
Prefix of the filename for the plot(s) to be saved to. If specified, plot filenames will be created by appending this prefix to  
either cpg.rank or cpg.name.  If not specified, plot will not be saved.
}
  \item{file.type}{
Type of file to be saved. Can either be \code{"pdf"} or \code{"eps"}. Selecting \code{file.type="eps"} will
result in publication quality editable postscript files that can be opened by Adobe Illustrator or Photoshop.
  }
  \item{eps.size}{
Vector indicating the size of .eps file (if creating one). Correponds to horrizontal and height.
}
  \item{popup.pdf}{
\code{TRUE} or \code{FALSE}. If creating a pdf file, this indicates if the plot should appear in a popup window as well. If running in a 
cluster-like environment, best to leave \code{FALSE}.
}
  \item{beta.values}{
If the object has been renamed (i.e. \code{x$info$betainfo} is no longer in \code{ls(.GlobalEnv)}) then specify the new object here.
}
  \item{main.title}{
Main title to be put on the graph. If \code{NULL} one based on the analysis will be used
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters.
}
}
\details{
An unlimited number of CpG sites can be selected for plotting by specifying either \option{cpg.rank} or \option{cpg.name}, as shown in the Examples below.  Note that only one of these options is needed; if  both are entered, \option{cpg.rank} will be used.
}
\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}




\seealso{
\code{\link{cpg.assoc}}
\code{\link{plot.cpg}}
\code{\link{manhattan}}
\code{\link{cpg.perm}}
\code{\link{plot.cpg.perm}}
}
\examples{
#Load the data:
data(samplecpg,samplepheno,package="CpGassoc")
###NOTE: If you are dealing with large data, do not specify large.data=FALSE.
###The default option is true.
##This will involve partitioning up the data and performing more gc() to clear up space
test<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
##Using rank, will plot the top three sites in order of significance:
scatterplot(test,c(1:3))
##Using name, specify three sites:
scatterplot(test,cpg.name=c("CpG1182","CpG1000","CpG42"))

##Plotting something that is categorical in nature:
test2<-cpg.assoc(samplecpg[1:200,],factor(samplepheno$Disease),large.data=FALSE)
scatterplot(test2,c(2),beta.values=samplecpg[1:200,])
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
