% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_GeneAccessors.R
\name{proteinSequence}
\alias{proteinSequence}
\alias{proteinSequence,ProteinCodingGene-method}
\title{Accessor function for proteinSequence}
\usage{
proteinSequence(gene)

\S4method{proteinSequence}{ProteinCodingGene}(gene)
}
\arguments{
\item{gene}{An instance of the ProteinCodingGene class}
}
\value{
The proteinSequence of the gene.

The Gene object with updated proteinSequence
}
\description{
Accessor function for proteinSequence
}
\examples{
library(GenomicRanges)
library(Biostrings)
pcg_example <- new("ProteinCodingGene",
  geneID = "ENSG00000139618",
  symbol = "BRCA2",
  name = "Breast cancer 2",
  description = "Gene implicated in breast cancer.",
  structure = GRanges(
    seqnames = Rle("17"),
    ranges = IRanges(start = 43044295, end = 43125482),
    strand = Rle("+")
  ),
  proteinID = "P51587",
  proteinSequence = AAString("MARKSLEMSIR")
)
proteinSequence(pcg_example)
}
