#' @importFrom methods setClass setGeneric setMethod
#' @importFrom GenomicRanges GRanges
NULL

#' Accessor function for geneID
#'
#' @param gene An instance of the Gene class
#'
#' @return the geneID of the gene
#' @export
setGeneric("geneID", function(gene) standardGeneric("geneID"))

#' @rdname geneID
#' @param gene An instance of the Gene class
#'
#' @return The Gene object with updated geneID
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' geneID(pcg_example)
setMethod("geneID", "Gene", function(gene) gene@geneID)


#' Accessor function for symbol
#'
#' @param gene An instance of the Gene class
#'
#' @return The symbol of the gene.
#' @export
setGeneric("symbol", function(gene) standardGeneric("symbol"))

#' @rdname symbol
#' @param gene An instance of the Gene class
#'
#' @return The Gene object with updated symbol
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' symbol(pcg_example)
setMethod("symbol", "Gene", function(gene) gene@symbol)


#' Accessor function for name
#'
#' @param gene An instance of the Gene class.
#'
#' @return The name of the gene.
#' @export
setGeneric("name", function(gene) standardGeneric("name"))

#' @rdname name
#' @param gene An instance of the Gene class.
#'
#' @return The Gene object with updated name
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' name(pcg_example)
setMethod("name", "Gene", function(gene) gene@name)


#' Accessor function for description
#'
#' @param gene An instance of the Gene class
#'
#' @return The description of the gene.
#' @export
setGeneric("description", function(gene) standardGeneric("description"))

#' @rdname description
#' @param gene An instance of the Gene class
#'
#' @return The Gene object with updated description
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' description(pcg_example)
setMethod("description", "Gene", function(gene) gene@description)

#' Accessor function for structure
#'
#' @param gene An instance of the Gene class.
#' @return The structure of the gene.
#' @export
setGeneric("structure", function(gene) standardGeneric("structure"))

#' @rdname structure
#' @param gene An instance of the Gene class
#'
#' @return The Gene object with updated structure
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' structure(pcg_example)
setMethod("structure", "Gene", function(gene) gene@structure)


#' Accessor function for proteinID
#'
#' @param gene An instance of the ProteinCodingGene class
#'
#' @return The proteinID of the gene.
#' @export
setGeneric("proteinID", function(gene) standardGeneric("proteinID"))

#' @rdname proteinID
#' @param gene An instance of the ProteinCodingGene class
#'
#' @return The Gene object with updated proteinID
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' proteinID(pcg_example)
setMethod("proteinID", "ProteinCodingGene", function(gene) gene@proteinID)

#' Accessor function for proteinSequence
#'
#' @param gene An instance of the ProteinCodingGene class
#'
#' @return The proteinSequence of the gene.
#' @export
setGeneric("proteinSequence", function(gene) standardGeneric("proteinSequence"))

#' @rdname proteinSequence
#' @param gene An instance of the ProteinCodingGene class
#'
#' @return The Gene object with updated proteinSequence
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
#'   proteinSequence = AAString("MARKSLEMSIR")
#' )
#' proteinSequence(pcg_example)
setMethod("proteinSequence",
          "ProteinCodingGene", function(gene) gene@proteinSequence)

#' Accessor function for lncRNAID
#'
#' @param gene An instance of the LongNonCodingRNAGene class
#'
#' @return The lncRNAID of the gene.
#' @export
setGeneric("lncRNAID", function(gene) standardGeneric("lncRNAID"))

#' @rdname lncRNAID
#' @param gene An instance of the LongNonCodingRNAGene class
#'
#' @return The Gene object with updated lncRNAID
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
#'     ranges = IRanges(start = 53947302, end = 53948402),
#'     strand = Rle("-")
#'   ),
#'   lncRNAID = "NR_00316",
#'   RNASequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' lncRNAID(lnc_example)
setMethod("lncRNAID", "LongNonCodingRNAGene", function(gene) gene@lncRNAID)

#' Accessor function for RNASequence
#'
#' @param gene An instance of the LongNonCodingRNAGene class
#'
#' @return The RNASequence of the gene.
#' @export
setGeneric("RNASequence", function(gene) standardGeneric("RNASequence"))

#' @rdname RNASequence
#' @param gene An instance of the LongNonCodingRNAGene class
#'
#' @return The Gene object with updated RNASequence
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
#'     ranges = IRanges(start = 53947302, end = 53948402),
#'     strand = Rle("-")
#'   ),
#'   lncRNAID = "NR_00316",
#'   RNASequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' RNASequence(lnc_example)
setMethod("RNASequence", 
          "LongNonCodingRNAGene", function(gene) gene@RNASequence)

#' Accessor function for piRNAID
#'
#' @param gene An instance of the PiwiInteractingRNAGene class
#'
#' @return The piRNAID of the gene.
#' @export
setGeneric("piRNAID", function(gene) standardGeneric("piRNAID"))

#' @rdname piRNAID
#' @param gene An instance of the PiwiInteractingRNAGene class
#'
#' @return The Gene object with updated piRNAID
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
#'     ranges = IRanges(start = 1547627, end = 1577799),
#'     strand = Rle("+")
#'   ),
#'   piRNAID = "PIR0001",
#'   piSequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' piRNAID(piRNA_example)
setMethod("piRNAID", "PiwiInteractingRNAGene", function(gene) {
    gene@piRNAID
})

#' Accessor function for piSequence
#'
#' @param gene An instance of the PiwiInteractingRNAGene class
#'
#' @return The piSequence of the gene.
#' @export
setGeneric("piSequence", function(gene) standardGeneric("piSequence"))

#' @rdname piSequence
#' @param gene An instance of the PiwiInteractingRNAGene class
#'
#' @return The Gene object with updated piSequence
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
#'     ranges = IRanges(start = 1547627, end = 1577799),
#'     strand = Rle("+")
#'   ),
#'   piRNAID = "PIR0001",
#'   piSequence = RNAString("UGAGGUAGUAGGUU")
#' )
#' piSequence(piRNA_example)
setMethod("piSequence", "PiwiInteractingRNAGene", function(gene) gene@piSequence)

#' Accessor function for microRNAID
#'
#' @param gene An instance of the MicroRNAGene class
#'
#' @return The microRNAID of the gene.
#' @export
setGeneric("microRNAID", function(gene) standardGeneric("microRNAID"))

#' @rdname microRNAID
#' @param gene An instance of the MicroRNAGene class
#'
#' @return The Gene object with updated microRNAID
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
#'     ranges = IRanges(start = 29570, end = 29691),
#'     strand = Rle("+")
#'   ),
#'   microRNAID = "MI000077",
#'   seedSequence = RNAString("UAGCUUAUCAGAC")
#' )
#' microRNAID(miRNA_example)
setMethod("microRNAID", "MicroRNAGene", function(gene) gene@microRNAID)


#' Accessor function for seedSequence
#'
#' @param gene An instance of the MicroRNAGene class
#'
#' @return The seedSequence of the gene.
#' @export
setGeneric("seedSequence", function(gene) standardGeneric("seedSequence"))

#' @rdname seedSequence
#' @param gene An instance of the MicroRNAGene class
#'
#' @return The Gene object with updated seedSequence
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
#'     ranges = IRanges(start = 29570, end = 29691),
#'     strand = Rle("+")
#'   ),
#'   microRNAID = "MI000077",
#'   seedSequence = RNAString("UAGCUUAUCAGAC")
#' )
#' seedSequence(miRNA_example)
setMethod("seedSequence", "MicroRNAGene", function(gene) gene@seedSequence)

#' Accessor function for rRNAID
#'
#' @param gene An instance of the RibosomialRNAGene class
#'
#' @return The rRNAID of the gene.
#' @export
setGeneric("rRNAID", function(gene) standardGeneric("rRNAID"))

#' @rdname rRNAID
#' @param gene An instance of the RibosomialRNAGene class
#'
#' @return The Gene object with updated rRNAID
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
#'     ranges = IRanges(start = 3413204, end = 3413561),
#'     strand = Rle("+")
#'   ),
#'   rRNAID = "RRNA0005",
#'   rRNASequence = RNAString("GCUCAGUUGGGAGAG")
#' )
#' rRNAID(rRNA_example)
setMethod("rRNAID", "RibosomialRNAGene", function(gene) gene@rRNAID)

#' Accessor function for rRNASequence
#'
#' @param gene An instance of the RibosomialRNAGene class
#'
#' @return The rRNASequence of the gene.
#' @export
setGeneric("rRNASequence", function(gene) standardGeneric("rRNASequence"))

#' @rdname rRNASequence
#' @param gene An instance of the RibosomialRNAGene class
#'
#' @return The Gene object with updated rRNASequence
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
#'     ranges = IRanges(start = 3413204, end = 3413561),
#'     strand = Rle("+")
#'   ),
#'   rRNAID = "RRNA0005",
#'   rRNASequence = RNAString("GCUCAGUUGGGAGAG")
#' )
#' rRNASequence(rRNA_example)
setMethod("rRNASequence", "RibosomialRNAGene", function(gene) gene@rRNASequence)

#' Accessor function for tRNAID
#'
#' @param gene An instance of the TransferRNAGene class
#'
#' @return The tRNAID of the gene.
#' @export
setGeneric("tRNAID", function(gene) standardGeneric("tRNAID"))

#' @rdname tRNAID
#' @param gene An instance of the TransferRNAGene class
#'
#' @return The Gene object with updated tRNAID
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
#'     ranges = IRanges(start = 26357, end = 29428),
#'     strand = Rle("-")
#'   ),
#'   tRNAID = "TRNA023",
#'   tRNASequence = RNAString("GCCUUUAGCUCAGUUGGGAGAG")
#' )
#' tRNAID(tRNA_example)
setMethod("tRNAID", "TransferRNAGene", function(gene) gene@tRNAID)

#' Accessor function for tRNASequence
#'
#' @param gene An instance of the TransferRNAGene class
#'
#' @return The tRNASequence of the gene.
#' @export
setGeneric("tRNASequence", function(gene) standardGeneric("tRNASequence"))

#' @rdname tRNASequence
#' @param gene An instance of the TransferRNAGene class
#'
#' @return The Gene object with updated tRNASequence
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
#'     ranges = IRanges(start = 26357, end = 29428),
#'     strand = Rle("-")
#'   ),
#'   tRNAID = "TRNA023",
#'   tRNASequence = RNAString("GCCUUUAGCUCAGUUGGGAGAG")
#' )
#' tRNASequence(tRNA_example)
setMethod("tRNASequence", "TransferRNAGene", function(gene) gene@tRNASequence)

#' Accessor function for snRNAID
#'
#' @param gene An instance of the SmallNuclearRNAGene class
#'
#' @return the snRNAID of the gene.
#' @export
setGeneric("snRNAID", function(gene) standardGeneric("snRNAID"))

#' @rdname snRNAID
#' @param gene An instance of the SmallNuclearRNAGene class
#'
#' @return The Gene object with updated snRNAID
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
#'     ranges = IRanges(start = 26357, end = 29428),
#'     strand = Rle("-")
#'   ),
#'   snRNAID = "SNR0001",
#'   snRNASequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' snRNAID(snRNA_example)
setMethod("snRNAID", "SmallNuclearRNAGene", function(gene) gene@snRNAID)

#' Accessor function for snRNASequence
#'
#' @param gene An instance of the SmallNuclearRNAGene class
#'
#' @return The snRNASequence of the gene.
#' @export
setGeneric("snRNASequence", function(gene) standardGeneric("snRNASequence"))

#' @rdname snRNASequence
#' @param gene An instance of the SmallNuclearRNAGene class
#'
#' @return The Gene object with updated snRNASequence
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
#'     ranges = IRanges(start = 26357, end = 29428),
#'     strand = Rle("-")
#'   ),
#'   snRNAID = "SNR0001",
#'   snRNASequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' snRNASequence(snRNA_example)
setMethod("snRNASequence", "SmallNuclearRNAGene", function(gene) gene@snRNASequence)

#' Accessor function for SRPRNAID
#'
#' @param gene An instance of the SRPRNAGene class
#'
#' @return The SRPRNAID of the gene.
#' @export
setGeneric("SRPRNAID", function(gene) standardGeneric("SRPRNAID"))

#' @rdname SRPRNAID
#' @param gene An instance of the SRPRNAGene class
#'
#' @return The Gene object with updated SRPRNAID
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
#'     ranges = IRanges(start = 112345, end = 112423),
#'     strand = Rle("-")
#'   ),
#'   SRPRNAID = "SNR0001",
#'   SRPSequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' SRPRNAID(SRP_example)
setMethod("SRPRNAID", "SRPRNAGene", function(gene) gene@SRPRNAID)

#' Accessor function for SRPSequence
#'
#' @param gene An instance of the SRPRNAGene class
#'
#' @return The SRPSequence of the gene.
#' @export
setGeneric("SRPSequence", function(gene) standardGeneric("SRPSequence"))

#' @rdname SRPSequence
#' @param gene An instance of the SRPRNAGene class
#'
#' @return The Gene object with updated SRPSequence
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
#'     ranges = IRanges(start = 112345, end = 112423),
#'     strand = Rle("-")
#'   ),
#'   SRPRNAID = "SNR0001",
#'   SRPSequence = RNAString("GCCUUUAGCUCAAGCGGCUUAGAG")
#' )
#' SRPSequence(SRP_example)
setMethod("SRPSequence", "SRPRNAGene", function(gene) gene@SRPSequence)
