#' @importFrom methods setClass new validObject
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings AAString RNAString
#' @importFrom S4Vectors Rle
NULL

#### Gene class ####

#' Create a Gene object
#'
#' This function initializes a new Gene object with specified attributes.
#'
#' @aliases Gene
#' @aliases createGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#'
#' @return A new Gene object.
#'
#' @export
Gene <- function(geneID, symbol, name, description, structure) {
    if (!inherits(structure, "GRanges")) {
        stop("Structure must be a GRanges object.")
    }
    new("Gene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure)
}

#### Protein Coding Gene class ####

#' Create a ProteinCodingGene object
#'
#' Constructor function for creating an instance of ProteinCodingGene class.
#'
#' @aliases ProteinCodingGeneGene
#' @aliases createProteinCodingGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param proteinID Unique identifier for the protein.
#' @param proteinSequence Amino acid sequence of the protein.
#'
#' @return A new ProteinCodingGene object.
#'
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
#'     ranges = IRanges(start = 500, end = 600),
#'     strand = Rle("+")
#'   ),
#'   proteinID = "P51587",
#'   proteinSequence = AAString("MVLSPADKTNVK")
#' )
#' pcg_example
#'
#' @export
ProteinCodingGene <- function(geneID, symbol, name, description, 
                              structure, proteinID, proteinSequence) {
    if (!is.character(proteinID) || !nzchar(proteinID)) {
        stop("proteinID must be a non-empty character string.")
    }
    if (!is.character(proteinSequence) || !nzchar(proteinSequence)) {
        stop("proteinSequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an AAString
    if (!inherits(proteinSequence, "AAString")) {
        tryCatch(
            {
                proteinSequence <- AAString(proteinSequence)
            },
            error = function(e) {
                stop("proteinSequence must be a valid protein sequence.")
            }
        )
    }
    new("ProteinCodingGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        proteinID = proteinID, proteinSequence = proteinSequence)
}

#### Long Non-coding RNA Gene class ####

#' Create a LongNonCodingRNAGene object
#'
#' Constructor function for creating an instance of LongNonCodingRNAGene class.
#'
#' @aliases LongNonCodingRNAGene
#' @aliases createLongNonCodingRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param lncRNAID Unique identifier for lncRNA.
#' @param RNASequence An RNAString containing the nucleotide sequence.
#'
#' @return A new LongNonCodingRNAGene object.
#'
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
#' print(lnc_example)
#'
#' @export
LongNonCodingRNAGene <- function(geneID, symbol, name, description, 
                                 structure, lncRNAID, RNASequence) {
    if (!is.character(lncRNAID) || !nzchar(lncRNAID)) {
        stop("lncRNAID must be a non-empty character string.")
    }
    if (!is.character(RNASequence) || !nzchar(RNASequence)) {
        stop("RNASequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(RNASequence, "RNAString")) {
        tryCatch(
            {
                RNASequence <- RNAString(RNASequence)
            },
            error = function(e) {
                stop("RNASequence must be a valid RNA sequence.")
            }
        )
    }
    new("LongNonCodingRNAGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        lncRNAID = lncRNAID, RNASequence = RNASequence)
}

#### Piwi interacting RNA Gene Class ####

#' Create a PiwiInteractingRNAGene object
#'
#' Constructor function for creating an instance of PiwiInteractingRNAGene class
#'
#' @aliases PiwiInteractingRNAGene
#' @aliases createPiwiInteractingRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param piRNAID Unique identifier for piRNA.
#' @param piSequence An RNAString containing the nucleotide sequence.
#'
#' @return A new PiwiInteractingRNAGene object.
#'
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
#' print(piRNA_example)
#'
#' @export
PiwiInteractingRNAGene <- function(geneID, symbol, name, description, 
                                   structure, piRNAID, piSequence) {
    if (!is.character(piRNAID) || !nzchar(piRNAID)) {
        stop("piRNAID must be a non-empty character string.")
    }
    if (!is.character(piSequence) || !nzchar(piSequence)) {
        stop("piSequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(piSequence, "RNAString")) {
        tryCatch(
            {
                piSequence <- RNAString(piSequence)
            },
            error = function(e) {
                stop("piSequence must be a valid RNA sequence.")
            }
        )
    }
    new("PiwiInteractingRNAGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        piRNAID = piRNAID, piSequence = piSequence)
}

#### Micro RNA Gene Class ####

#' Create a MicroRNAGene object
#'
#' Constructor function for creating an instance of MicroRNAGene class.
#'
#' @aliases MicroRNAGene
#' @aliases createMicroRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param microRNAID Unique identifier for microRNA.
#' @param seedSequence An RNAString containing the nucleotide sequence.
#'
#' @return A new MicroRNAGene object.
#'
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
#' print(miRNA_example)
#'
#' @export
MicroRNAGene <- function(geneID, symbol, name, description, 
                         structure, microRNAID, seedSequence) {
    if (!is.character(microRNAID) || !nzchar(microRNAID)) {
        stop("microRNAID must be a non-empty character string.")
    }
    if (!is.character(seedSequence) || !nzchar(seedSequence)) {
        stop("seedSequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(seedSequence, "RNAString")) {
        tryCatch(
            {
                seedSequence <- RNAString(seedSequence)
            },
            error = function(e) {
                stop("seedSequence must be a valid RNA sequence.")
            }
        )
    }
    new("MicroRNAGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        microRNAID = microRNAID, seedSequence = seedSequence)
}

#### Ribosomial RNA Gene class ####

#' Create a RibosomialRNAGene object
#'
#' Constructor function for creating an instance of RibosomialRNAGene class.
#'
#' @aliases RibosomialRNAGene
#' @aliases createRibosomialRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param rRNAID Unique identifier for rRNA.
#' @param rRNASequence An RNAString containing the nucleotide sequence.
#'
#' @return A new RibosomialRNAGene object.
#'
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
#' print(rRNA_example)
#'
#' @export
RibosomialRNAGene <- function(geneID, symbol, name, description, 
                              structure, rRNAID, rRNASequence) {
    if (!is.character(rRNAID) || !nzchar(rRNAID)) {
        stop("rRNAID must be a non-empty character string.")
    }
    if (!is.character(rRNASequence) || !nzchar(rRNASequence)) {
        stop("rRNASequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(rRNASequence, "RNAString")) {
        tryCatch(
            {
                rRNASequence <- RNAString(rRNASequence)
            },
            error = function(e) {
                stop("rRNASequence must be a valid RNA sequence.")
            }
        )
    }
    new("RibosomialRNAGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        rRNAID = rRNAID, rRNASequence = rRNASequence)
}


#### Transfer RNA Gene Class ####

#' Create a TransferRNAGene object
#'
#' Constructor function for creating an instance of TransferRNAGene class.
#'
#' @aliases TransferRNAGene
#' @aliases createTransferRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param tRNAID Unique identifier for tRNA.
#' @param tRNASequence An RNAString containing the nucleotide sequence.
#'
#' @return A new TransferRNAGene object.
#'
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
#' print(tRNA_example)
#'
#' @export
TransferRNAGene <- function(geneID, symbol, name, description, 
                            structure, tRNAID, tRNASequence) {
    if (!is.character(tRNAID) || !nzchar(tRNAID)) {
        stop("tRNAID must be a non-empty character string.")
    }
    if (!is.character(tRNASequence) || !nzchar(tRNASequence)) {
        stop("tRNASequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(tRNASequence, "RNAString")) {
        tryCatch(
            {
                tRNASequence <- RNAString(tRNASequence)
            },
            error = function(e) {
                stop("tRNASequence must be a valid RNA sequence.")
            }
        )
    }
    new("TransferRNAGene", geneID = geneID, symbol = symbol, name = name, 
        description = description, structure = structure, 
        tRNAID = tRNAID, tRNASequence = tRNASequence)
}

#### Small Nuclear RNA Gene Class ####

#' Create a SmallNuclearRNAGene object
#'
#' Constructor function for creating an instance of SmallNuclearRNAGene class.
#'
#' @aliases SmallNuclearRNAGene
#' @aliases createSmallNuclearRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param snRNAID Unique identifier for snRNA.
#' @param snRNASequence An RNAString containing the nucleotide sequence.
#'
#' @return A new SmallNuclearRNAGene object.
#'
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
#' print(snRNA_example)
#'
#' @export
SmallNuclearRNAGene <- function(geneID, symbol, name, description,
                                structure, snRNAID, snRNASequence) {
    if (!is.character(snRNAID) || !nzchar(snRNAID)) {
        stop("snRNAID must be a non-empty character string.")
    }
    if (!is.character(snRNASequence) || !nzchar(snRNASequence)) {
        stop("snRNASequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(snRNASequence, "RNAString")) {
        tryCatch(
            {
                snRNASequence <- RNAString(snRNASequence)
            },
            error = function(e) {
                stop("snRNASequence must be a valid RNA sequence.")
            }
        )
    }
    new("SmallNuclearRNAGene", geneID = geneID, symbol = symbol, name = name,
        description = description, structure = structure,
        snRNAID = snRNAID, snRNASequence = snRNASequence)
}

#### Signal Recognition Particle RNA Gene Class ####

#' Create a SRPRNAGene object
#'
#' Constructor function for creating an instance of SRPRNAGene class.
#'
#' @aliases SRPRNAGene
#' @aliases createSRPRNAGene
#'
#' @param geneID Unique identifier for the gene.
#' @param symbol Standard symbol associated with the gene.
#' @param name Full name of the gene.
#' @param description Detailed description of the gene.
#' @param structure A GRanges object of the genomic structure of the gene.
#' @param SRPRNAID Unique identifier for SRPRNA.
#' @param SRPSequence An RNAString containing the nucleotide sequence.
#'
#' @return A new SRPRNAGene object.
#'
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
#' print(SRP_example)
#'
#' @export
SRPRNAGene <- function(geneID, symbol, name, description,
                       structure, SRPRNAID, SRPSequence) {
    if (!is.character(SRPRNAID) || !nzchar(SRPRNAID)) {
        stop("SRPRNAID must be a non-empty character string.")
    }
    if (!is.character(SRPSequence) || !nzchar(SRPSequence)) {
        stop("SRPSequence must be a non-empty character string.")
    }
    # Attempt to convert the character string to an RNAString
    if (!inherits(SRPSequence, "RNAString")) {
        tryCatch(
            {
                SRPSequence <- RNAString(SRPSequence)
            },
            error = function(e) {
                stop("SRPSequence must be a valid RNA sequence.")
            }
        )
    }
    new("SRPRNAGene", geneID = geneID, symbol = symbol, name = name,
        description = description, structure = structure,
        SRPRNAID = SRPRNAID, SRPSequence = SRPSequence)
}
