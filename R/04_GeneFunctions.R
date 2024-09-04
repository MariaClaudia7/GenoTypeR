#' @importFrom methods setMethod setGeneric
#' @importFrom Biostrings AAString RNAString
NULL

#' Product length of a gene
#'
#' @param gene The gene object from which the product length will be calculated.
#'
#' @return Returns the length of the gene product.
#' @export
setGeneric("lengthProduct", function(gene) standardGeneric("lengthProduct"))

#' Product length for Protein Coding Gene
#'
#' @param gene A ProteinCodingGene object
#'
#' @return Integer representing the length of the protein sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' pcg_example <- new("ProteinCodingGene",
#'   geneID = "ENSG00000139618",
#'   symbol = "BRCA2",
#'   name = "Breast cancer 2",
#'   description = "Gene implicated in breast cancer.",
#'   structure = GRanges(
#'     seqnames = Rle("17"),
#'     ranges = IRanges(start = 43044295, end = 43125482),
#'     strand = Rle("+")
#'   ),
#'   proteinID = "P51587",
#'   proteinSequence = AAString("MVLSPADKTNVK")
#' )
#' pcg_length <- lengthProduct(pcg_example)
#' print(pcg_length)
setMethod("lengthProduct", signature = "ProteinCodingGene", function(gene) {
    nchar(as.character(gene@proteinSequence))
})

#' Length of the RNA product for LongNonCodingRNAGene.
#'
#' @param gene A `LongNonCodingRNAGene` object
#'
#' @return Integer representing the length of the RNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' lnc_example <- new("LongNonCodingRNAGene",
#'   geneID = "ENSG00000228630",
#'   symbol = "HOTAIR",
#'   name = "HOX transcript antisense RNA",
#'   description = "LncRNA implicated in gene silencing.",
#'   structure = GRanges(
#'     seqnames = Rle("12"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("-")
#'   ),
#'   lncRNAID = "NR_00316",
#'   RNASequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' lnc_length <- lengthProduct(lnc_example)
#' print(lnc_length)
setMethod("lengthProduct", signature = "LongNonCodingRNAGene", function(gene) {
    nchar(as.character(gene@RNASequence))
})

#' Length of the piRNA sequence for Piwi-interactingRNAGene.
#'
#' @param gene A PiwiInteractingRNAGene object
#'
#' @return Integer representing the length of the piRNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' piRNA_example <- new("PiwiInteractingRNAGene",
#'   geneID = "ENSG00000241599",
#'   symbol = "PIWI001",
#'   name = "Piwi-interacting RNA 1",
#'   description = "Piwi-interacting RNA involved in gene regulation",
#'   structure = GRanges(
#'     seqnames = Rle("x"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("+")
#'   ),
#'   piRNAID = "PIR0001",
#'   piSequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' piRNA_length <- lengthProduct(piRNA_example)
#' print(piRNA_length)
setMethod("lengthProduct",
          signature = "PiwiInteractingRNAGene", function(gene) {
              nchar(as.character(gene@piSequence))
          })

#' Length of the RNA sequence for MicroRNAGene.
#'
#' @param gene A MicroRNAGene object
#'
#' @return Integer representing the length of the RNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' miRNA_example <- new("MicroRNAGene",
#'   geneID = "ENSG00000207734",
#'   symbol = "MIR21",
#'   name = "MicroRNA 21",
#'   description = "miRNA involved in cancer regulation and expression.",
#'   structure = GRanges(
#'     seqnames = Rle("17"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("+")
#'   ),
#'   microRNAID = "MI000077",
#'   seedSequence = RNAString("UAGCUUAUCAGAC")
#' )
#' miRNA_length <- lengthProduct(miRNA_example)
#' print(miRNA_length)
setMethod("lengthProduct", signature = "MicroRNAGene", function(gene) {
    nchar(as.character(gene@seedSequence))
})

#' Length of the rRNA sequence for RibosomialRNAGene.
#'
#' @param gene A RibosomialRNAGene`object
#'
#' @return Integer representing the length of the rRNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' rRNA_example <- new("RibosomialRNAGene",
#'   geneID = "ENSG00000279457",
#'   symbol = "rRNA-5s",
#'   name = "5s ribosomial RNA",
#'   description = "rRNA component of the ribosome 5S subunit",
#'   structure = GRanges(
#'     seqnames = Rle("1"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("+")
#'   ),
#'   rRNAID = "RRNA0005",
#'   rRNASequence = RNAString("GCUCAGUUGGGAGAG")
#' )
#' rRNA_length <- lengthProduct(rRNA_example)
#' print(rRNA_length)
setMethod("lengthProduct", signature = "RibosomialRNAGene", function(gene) {
    nchar(as.character(gene@rRNASequence))
})

#' Length of the tRNA sequence for TransferRNAGene.
#'
#' @param gene A `TransferRNAGene` object
#'
#' @return Integer representing the length of the tRNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' tRNA_example <- new("TransferRNAGene",
#'   geneID = "ENSG00000201145",
#'   symbol = "tRNA-Lys",
#'   name = "Transfer RNA lysine",
#'   description = "tRNA involved in protein synthesis",
#'   structure = GRanges(
#'     seqnames = Rle("Y"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("-")
#'   ),
#'   tRNAID = "TRNA023",
#'   tRNASequence = RNAString("GCCUUUAGCUCAGUUGGGAGAG")
#' )
#' tRNA_length <- lengthProduct(tRNA_example)
#' print(tRNA_length)
setMethod("lengthProduct", signature = "TransferRNAGene", function(gene) {
    nchar(as.character(gene@tRNASequence))
})

#' Length of the snRNA sequence for SmallNuclearRNAGene.
#'
#' @param gene A SmallNuclearRNAGene object
#'
#' @return Integer representing the length of the snRNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' snRNA_example <- new("SmallNuclearRNAGene",
#'   geneID = "ENSG00000234745",
#'   symbol = "snRNP-U1",
#'   name = "Small nuclear ribonucleoprotein U1",
#'   description = "snRNA component of the spliceosome",
#'   structure = GRanges(
#'     seqnames = Rle("Y"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("-")
#'   ),
#'   snRNAID = "SNR0001",
#'   snRNASequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' snRNA_length <- lengthProduct(snRNA_example)
#' print(snRNA_length)
setMethod("lengthProduct", signature = "SmallNuclearRNAGene", function(gene) {
    nchar(as.character(gene@snRNASequence))
})

#' Length of the snRNA sequence for Signal recognition particle RNAGene.
#'
#' @param gene A SRPRNAGene object
#'
#' @return Integer representing the length of the SRP RNA sequence.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' SRP_example <- new("SRPRNAGene",
#'   geneID = "ENSG00000248527",
#'   symbol = "SRP-RNA",
#'   name = "Signal recognition particle RNA",
#'   description = "SRP RNA involved in targeting secretory proteins",
#'   structure = GRanges(
#'     seqnames = Rle("5"),
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("-")
#'   ),
#'   SRPRNAID = "SNR0001",
#'   SRPSequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' SRP_length <- lengthProduct(SRP_example)
#' print(SRP_length)
setMethod("lengthProduct", signature = "SRPRNAGene", function(gene) {
    nchar(as.character(gene@SRPSequence))
})
