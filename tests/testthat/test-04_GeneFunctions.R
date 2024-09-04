library(testthat)
library(GenoTypeR)
library(GenomicRanges)
library(Biostrings)

### ProteinCodingGene Tests

test_that("ProteinCodingGene class functionality", {
    gr <- GRanges(seqnames = Rle("chr1"),
                  ranges = IRanges(start = 100, end = 200), strand = Rle("+"))
    pcg <- new("ProteinCodingGene",
               geneID = "PCG123", symbol = "PCG",
               name = "ProteinCodingGene",
               description = "A protein-coding gene",
               structure = gr, proteinID = "P123",
               proteinSequence = AAString("MKTEVDQ")
    )
    
    p_length <- lengthProduct(pcg)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("MKTEVDQ")))
})

### LongNonCodingRNAGene Tests

test_that("LongNonCodingRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chr1"),
                  ranges = IRanges(start = 100, end = 200), strand = Rle("+"))
    lnc <- new("LongNonCodingRNAGene",
               geneID = "LNC123", symbol = "LNCG",
               name = "LongNonCodingRNAGene",
               description = "A long non-coding RNA gene",
               structure = gr, lncRNAID = "P123",
               RNASequence = RNAString("AUCG")
    )
    
    p_length <- lengthProduct(lnc)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("AUCG")))
})

### PiwiInteractingRNAGene Tests

test_that("PiwiInteractingRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chr1"),
                  ranges = IRanges(start = 100, end = 200), strand = Rle("+"))
    pig <- new("PiwiInteractingRNAGene",
               geneID = "PI123", symbol = "PIG",
               name = "PiwiInteractingRNAGene",
               description = "A piwi interacting RNA gene",
               structure = gr, piRNAID = "P123",
               piSequence = RNAString("AUCG")
    )
    
    p_length <- lengthProduct(pig)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("AUCG")))
})

### MicroRNAGene Tests

test_that("MicroRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chrX"),
                  ranges = IRanges(start = 1000, end = 1100), strand = Rle("-"))
    micro_gene <- new("MicroRNAGene",
                      geneID = "MG001", symbol = "miR-21",
                      name = "Micro RNA 21",
                      description = "A Micro RNA gene",
                      structure = gr, microRNAID = "MIR021",
                      seedSequence = RNAString("CUUACGAUAG")
    )
    
    p_length <- lengthProduct(micro_gene)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("CUUACGAUAG")))
})

### TransferRNAGene Tests

test_that("TransferRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chrY"),
                  ranges = IRanges(start = 500, end = 600), strand = Rle("+"))
    trna_gene <- new("TransferRNAGene",
                     geneID = "TG001", symbol = "tRNA-Gly",
                     name = "Transfer RNA Glycine",
                     description = "tRNA produce the glycine amino acid",
                     structure = gr, tRNAID = "TRNA001", 
                     tRNASequence = RNAString("GCCUAGCUAG")
    )
    
    p_length <- lengthProduct(trna_gene)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("GCCUAGCUAG")))
})

### SmallNuclearRNAGene Tests

test_that("SmallNuclearRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chr2"),
                  ranges = IRanges(start = 200, end = 300), strand = Rle("-"))
    snrna_gene <- new("SmallNuclearRNAGene",
                      geneID = "SN001", symbol = "snRNA-U6",
                      name = "Small Nuclear RNA U6",
                      description = "snRNA involved in splicing",
                      structure = gr, snRNAID = "SNRNA001",
                      snRNASequence = RNAString("AUAUCGAU")
    )
    
    p_length <- lengthProduct(snrna_gene)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("AUAUCGAU")))
})

### SRPRNAGene Tests

test_that("SRPRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chr22"),
                  ranges = IRanges(start = 2200, end = 2300), strand = Rle("+"))
    srp_gene <- new("SRPRNAGene",
                    geneID = "SRP001", symbol = "SRP RNA",
                    name = "Signal Recognition Particle RNA",
                    description = "RNA component of the SRP",
                    structure = gr, SRPRNAID = "SRPNA001",
                    SRPSequence = RNAString("GUCAGUCAGU")
    )
    
    p_length <- lengthProduct(srp_gene)
    expect_true(is(p_length, "integer"))
    expect_equal(p_length, nchar(as.character("GUCAGUCAGU")))
})
