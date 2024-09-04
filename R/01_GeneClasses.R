#' @importFrom GenomicRanges GRanges
#' @importFrom methods setClass
#' @importFrom IRanges IRanges
NULL

#' Gene Class
#'
#' This virtual class represents the fundamental attributes for a gene.
#'
#' @title Gene
#'
#' @slot geneID Unique gene identifier
#' @slot symbol gene symbol (e.g., HUGO symbol).
#' @slot name full gene name
#' @slot description Brief gene description 
#' @slot structure Genomic structure as GRanges
#' 
#' @return No objects are returned because gene class is virtual
setClass(
    "Gene",
    slots = list(
        geneID = "character",
        symbol = "character",
        name = "character",
        description = "character",
        structure = "GRanges"
    ),
    contains = "VIRTUAL",
)

#' Protein-Coding Gene Class
#'
#' S4 class representing a protein-coding gene.
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title ProteinCodingGene
#'
#' @slot proteinID Unique protein identifier
#' @slot proteinSequence Amino acid sequence (AAString)
#' 
#' @return Returns an object of class \code{ProteinCodingGene}
#' 
#' @export
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' pcg_example <- new("ProteinCodingGene",
#' geneID = "ENSG00000139618",
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
#' print(pcg_example)
setClass(
    "ProteinCodingGene",
    slots = list(
        proteinID = "character",
        proteinSequence = "AAString"
    ),
    contains = "Gene"
)

#' Long non-coding RNA Gene Class
#'
#' S4 class representing a long non-coding gene.
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title Long non-coding Gene class
#'
#' @slot lncRNAID Unique lncRNA identifier
#' @slot RNASequence Nucleotide sequence of lncRNA.
#' 
#' @return Returns an object of class \code{LongNonCodingRNAGene}
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
#' print(lnc_example)
setClass(
    "LongNonCodingRNAGene",
    slots = list(
        lncRNAID = "character",
        RNASequence = "RNAString"
    ),
    contains = "Gene"
)

#' Piwi interacting RNA Gene Class
#'
#' S4 class representing piwi-interacting RNA gene.
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title PiwiInteractingRNAGene
#'
#' @slot piRNAID Unique piRNA identifier
#' @slot piSequence Nucleotide piRNA sequence
#' 
#' @return Returns an object of class \code{PiwiInteractingRNAGene}
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
#' print(piRNA_example)
setClass(
    "PiwiInteractingRNAGene",
    slots = list(
        piRNAID = "character",
        piSequence = "RNAString"
    ),
    contains = "Gene"
)

#' Micro RNA Gene Class
#'
#' S4 class representing microRNA gene.
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title MicroRNAGene
#'
#' @slot microRNAID Unique microRNA identifier
#' @slot seedSequence Nucleotide microRNA sequence
#' 
#' @return Returns an object of class \code{MicroRNAGene}
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
#' print(miRNA_example)
setClass(
    "MicroRNAGene",
    slots = list(
        microRNAID = "character",
        seedSequence = "RNAString"
    ),
    contains = "Gene"
)

#' Ribosomial RNA Gene Class
#'
#' S4 class representing microRNA gene.
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title RibosomialRNAGene
#'
#' @slot rRNAID Unique microRNA identifier
#' @slot rRNASequence Nucleotide microRNA sequence
#' 
#' @return Returns an object of class \code{RibosomialRNAGene}
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
#' print(rRNA_example)
setClass(
    "RibosomialRNAGene",
    slots = list(
        rRNAID = "character",
        rRNASequence = "RNAString"
    ),
    contains = "Gene"
)

#' Transfer RNA Gene Class
#'
#' This S4 class represents transferRNA gene
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title TransferRNAGene
#'
#' @slot tRNAID Unique tRNA identifier
#' @slot tRNASequence Nucleotide tRNA sequence
#' 
#' @return Returns an object of class \code{TransferRNAGene}
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
#' print(tRNA_example)
setClass(
    "TransferRNAGene",
    slots = list(
        tRNAID = "character",
        tRNASequence = "RNAString"
    ),
    contains = "Gene"
)

#' Small Nuclear RNA Gene Class
#'
#' This S4 class represents smallnuclearRNA gene
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title SmallNuclearRNAGene
#'
#' @slot snRNAID Unique snRNA identifier
#' @slot snRNASequence Nucleotide snRNA sequence
#' 
#' @return Returns an object of class \code{SmallNuclearRNAGene}
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
#' print(snRNA_example)
setClass(
    "SmallNuclearRNAGene",
    slots = list(
        snRNAID = "character",
        snRNASequence = "RNAString"
    ),
    contains = "Gene"
)

#' Signal Recognition Particle RNA Gene Class
#'
#' This S4 class represents smallnuclearRNA gene
#' Inherits all slots from \code{\link{Gene}}.
#'
#' @title SRPRNAGene
#'
#' @slot SRPRNAID Unique snRNA identifier
#' @slot SRPSequence Nucleotide snRNA sequence
#' 
#' @return Returns an object of class \code{SRPRNAGene}
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
#' print(SRP_example)
setClass(
    "SRPRNAGene",
    slots = list(
        SRPRNAID = "character",
        SRPSequence = "RNAString"
    ),
    contains = "Gene"
)
