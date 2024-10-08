% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_GeneAccessors.R
\name{snRNASequence}
\alias{snRNASequence}
\alias{snRNASequence,SmallNuclearRNAGene-method}
\title{Accessor function for snRNASequence}
\usage{
snRNASequence(gene)

\S4method{snRNASequence}{SmallNuclearRNAGene}(gene)
}
\arguments{
\item{gene}{An instance of the SmallNuclearRNAGene class}
}
\value{
The snRNASequence of the gene.

The Gene object with updated snRNASequence
}
\description{
Accessor function for snRNASequence
}
\examples{
library(GenomicRanges)
library(Biostrings)
snRNA_example <- new("SmallNuclearRNAGene",
  geneID = "ENSG00000234745",
  symbol = "snRNP-U1",
  name = "Small nuclear ribonucleoprotein U1",
  description = "snRNA component of the spliceosome",
  structure = GRanges(
    seqnames = Rle("Y"),
    ranges = IRanges(start = 26357, end = 29428),
    strand = Rle("-")
  ),
  snRNAID = "SNR0001",
  snRNASequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
)
snRNASequence(snRNA_example)
}
