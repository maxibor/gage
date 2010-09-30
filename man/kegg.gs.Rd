\name{kegg.gs}
\Rdversion{1.1}
\alias{kegg.gs}
\alias{go.gs}
\alias{carta.gs}
\docType{data}
\title{
  Common gene set data collections
}
\description{
  The gene set data collections derived from KEGG, GO
  and BioCarta databases.
}
\usage{
data(kegg.gs)
data(go.gs)
data(carta.gs)
}
\format{
 kegg.gs is a named list of 205 elements. Each element is a character
 vector of member gene Entrez IDs for a single KEGG pathway. Type
 \code{head(kegg.gs, 3)} for the first 3 gene sets or pathways.

 go.gs is a named list of 1000 elements in this package. It is
 a truncated list in this package. The ful list of go.gs has
 ~10000 elements and is provided with an experimental data package
 gageData. Each element is a character
 vector of member gene Entrez IDs for a single Gene Ontology term. Type
 \code{head(go.gs, 3)} for the first 3 gene sets or GO terms.

 carta.gs is a named list of 259 elements. Each element is a character
 vector of member gene Entrez IDs for a single BioCarta pathway. Type
 \code{head(carta.gs, 3)} for the first 3 gene sets or pathways.
}
\details{
These gene set data were compiled using Entrez Gene IDs, gene set names
and mapping information from multiple Bioconductor packages, including:
org.Hs.eg.db, kegg.db, go.db and cMAP. Please check the corresponding
packages for more information.

We only provide gene set data for human with this package. For other
species, please check the experiment data package list of Bioconductor
website or use the Bioconductor package GSEABase to build the users' own
gene set collections.
}
\source{
  Data from multiple Bioconductor packages, including: org.Hs.eg.db,
  kegg.db, go.db and cMAP.
}
\references{
  Entrez Gene <URL: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene>
  KEGG pathways <URL: ftp://ftp.genome.ad.jp/pub/kegg/pathways>
  Gene Ontology <URL: http://www.geneontology.org/>
  cMAP <URL: http://cmap.nci.nih.gov/PW>
  }

\examples{
#load expression and gene set data
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)

data(kegg.gs)
data(go.gs)

#make sure the gene IDs are the same for expression data and gene set
#data
lapply(kegg.gs[1:6],head)
lapply(go.gs[1:6],head)
head(rownames(gse16873))

#GAGE analysis
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis)
gse16873.go.p <- gage(gse16873, gsets = go.gs,
    ref = hn, samp = dcis)
}

\keyword{datasets}
