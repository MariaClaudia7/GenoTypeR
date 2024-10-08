% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_GeneClasses.R
\docType{class}
\name{PiwiInteractingRNAGene-class}
\alias{PiwiInteractingRNAGene-class}
\title{PiwiInteractingRNAGene}
\value{
Returns an object of class \code{PiwiInteractingRNAGene}
}
\description{
Piwi interacting RNA Gene Class
}
\details{
S4 class representing piwi-interacting RNA gene.
Inherits all slots from \code{\link{Gene}}.
}
\section{Slots}{

\describe{
\item{\code{piRNAID}}{Unique piRNA identifier}

\item{\code{piSequence}}{Nucleotide piRNA sequence}
}}

\examples{
library(GenomicRanges)
library(Biostrings)
piRNA_example <- new("PiwiInteractingRNAGene",
  geneID = "ENSG00000241599",
  symbol = "PIWI001",
  name = "Piwi-interacting RNA 1",
  description = "Piwi-interacting RNA involved in gene regulation",
  structure = GRanges(
    seqnames = Rle("x"),
    ranges = IRanges(start = 1547627, end = 1577799),
    strand = Rle("+")
  ),
  piRNAID = "PIR0001",
  piSequence = RNAString("UGAGGUAGUAGGUU")
)
print(piRNA_example)
}
