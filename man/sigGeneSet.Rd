\name{sigGeneSet}
\Rdversion{1.1}
\alias{sigGeneSet}

\title{
Significant gene set from GAGE analysis
}
\description{
  This function sorts and counts signcant gene sets based on q- or p-value
  cutoff.
}
\usage{
sigGeneSet(setp, cutoff = 0.1, dualSig = (0:2)[2], qpval = c("q.BH",
"P.erlang")[1])
}

\arguments{
  \item{setp}{
    the result object returned by \code{gage} function, either a numeric
    matrix or a list of two such matrices. Check \code{gage} help
    information for details.
}
  \item{cutoff}{
    numeric, q- or p-value cutoff, between 0 and 1. Default 0.1 (for
    q-value). When p-value is used, recommended cutoff value is 0.001 for data with more than 2
    replicates per condition or 0.01 for les sample sizes.
}
  \item{dualSig}{
    integer, switch argument controlling how dual-signficant gene sets
    should be treated. 0: discard such gene sets from
    the final significant gene set list; 1: keep such gene sets in the
    more significant direction and remove them from the less significant
    direction; 2: keep such gene sets in the lists for both
    directions. default to 1. Dual-signficant means a gene set is called
    significant simultaneously in both 1-direction tests (up- and
    down-regulated). Check the details for more information.
}
  \item{qpval}{
    character, specifies the column name used for gene set selection,
    i.e. what type of q- or p-value to use in gene set
    selection. Default to 
    be "q.BH" (q-value using BH procedure). "P.erlang" is the unadjusted
    global p-value and may be used as selection criterion
    sometimes. 
}
}
\details{
    Dual signficant gene sets are rare cases, but do occurred in large
    clinic datasets. Gene sets are signficantly up-regulated in a subset
    of experiments, but down-regulated in another subset. Note that dual-signficant gene sets are not the
    same as gene sets called signficant in 2-directional tests, although
    they are related.
}
\value{
  \code{sigGeneSet} function returns a data matrix of the same structure
  as \code{gage} result matrix. Check \code{gage} help information for
  details. This data matrix combined both up-regulated and
  down-regulated gene set data for 1-directional tests.
}
\references{
  Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE:
  Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC
  Bioinformatics 2009, 10:161
}
\author{
  Weijun Luo <luo@cshl.edu>
}

\seealso{
  \code{\link{gage}} the main function for GAGE analysis;
  \code{\link{esset.grp}} non-redundant signcant gene set list;
  \code{\link{essGene}} essential member genes in a gene set;
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
data(kegg.gs)

#kegg test for 1-directional changes
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis)
#kegg test for 2-directional changes
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis, same.dir = FALSE)
gse16873.kegg.sig<-sigGeneSet(gse16873.kegg.p, dualSig=1)
str(gse16873.kegg.sig)
gse16873.kegg.2d.sig<-sigGeneSet(gse16873.kegg.2d.p)
str(gse16873.kegg.2d.sig)
}

\keyword{htest}
\keyword{multivariate}
\keyword{manip}
