\name{manhattan.reflect}
\alias{manhattan.reflect}

\title{
Create a Reflective Manhattan plot
}
\description{
This function will produce a reflective manhattan plot for the observed P-values from an object of class \code{"cpg"} or \code{"cpg.perm"}.
The original analysis needs to be performed on a continuous variable (need T-statistics).
}
\usage{
manhattan.reflect(x, cpgname, chr, pos, save.plot = NULL, file.type="pdf",
popup.pdf = FALSE, eps.size = c(15, 5), main.title = NULL, cpg.labels = NULL,
chr.list = NULL, color.list = NULL, fdr.cutoff=NULL, point.size=NULL, ...)
}

\arguments{
  \item{x}{
Object of class \code{"cpg"} or \code{"cpg.perm"}.
}
  \item{cpgname}{
A vector consisting of the labels for each CpG site.
}
  \item{chr}{
A vector consisting of the chromosome number for each CpG site.
}
  \item{pos}{
The map position of each CpG site within its chromosome.
}
  \item{save.plot}{
Name of the file for the plot to be saved to. If not specified, plot will not be saved.
}
  \item{file.type}{
Type of file to be saved. Can either be \code{"pdf"} or \code{"eps"}. Selecting \code{file.type="eps"} will
result in publication quality editable postscript files that can be opened by Adobe Illustrator or Photoshop.
  }
  \item{popup.pdf}{
\code{TRUE} or \code{FALSE}. If creating a pdf file, this indicates if the plot should appear in a popup window as well. If running in a
cluster-like environment, best to leave \code{FALSE}.
}
  \item{eps.size}{
Vector indicating the size of .eps file (if creating one). Corresponds to horizontal and height.
}
  \item{main.title}{
Main title to be put on the graph. If \code{NULL} one based on the analysis will be used.
}
  \item{cpg.labels}{
A character scalar of either \code{"FDR"} or \code{"HOLM"} which will label the significant sites on the manhattan plot.
}
  \item{chr.list}{
A vector listing the chromosomes to be plotted (all available chromosomes are plotted by default). The X and Y chromosomes
can be denoted by 23 and 24.
}
  \item{color.list}{
A vector of custom colors to be used for each chromosomes in the manhattan plot.
}
  \item{fdr.cutoff}{
A numeric scalar between 0 and 1 to indicate what to consider FDR significant. Defaults to \code{NULL}.
}
  \item{point.size}{
The size of the points in the manhattan plot, if NULL, default to our default, where significant probes have different sizes.  
}  
 \item{\dots}{
Arguments to be passed to methods, such as graphical parameters.
}
}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <barfieldrichard8@gmail.com>
}
\note{
\option{cpgname}, \option{chr}, and \option{pos} must be sorted in the same order, so that the first cpgname[1] corresponds to chr[1] and pos[1], and so on.
}



\seealso{
\code{\link{cpg.perm}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.assoc}}
\code{\link{plot.cpg.perm}}
\code{\link{manhattan}}
}
\examples{
#Doing a Manhattan plot. First load the data:
data(samplecpg,samplepheno,annotation,package="CpGassoc")

examplemanhat<-cpg.assoc(samplecpg,samplepheno$Disease,large.data=FALSE)

manhattan.reflect(examplemanhat,annotation$TargetID,annotation$CHR,annotation$MAPINFO)


}
