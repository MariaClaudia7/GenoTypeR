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
    
    expect_equal(geneID(pcg), "PCG123")
    expect_equal(symbol(pcg), "PCG")
    expect_equal(name(pcg), "ProteinCodingGene")
    expect_equal(description(pcg), "A protein-coding gene")
    expect_true(is(structure(pcg), "GRanges"))
    expect_equal(proteinID(pcg), "P123")
    expect_true(is(proteinSequence(pcg), "AAString"))
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
    
    expect_equal(geneID(lnc), "LNC123")
    expect_equal(symbol(lnc), "LNCG")
    expect_equal(name(lnc), "LongNonCodingRNAGene")
    expect_equal(description(lnc), "A long non-coding RNA gene")
    expect_true(is(structure(lnc), "GRanges"))
    expect_equal(lncRNAID(lnc), "P123")
    expect_true(is(RNASequence(lnc), "RNAString"))
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
    
    expect_equal(geneID(pig), "PI123")
    expect_equal(symbol(pig), "PIG")
    expect_equal(name(pig), "PiwiInteractingRNAGene")
    expect_equal(description(pig), "A piwi interacting RNA gene")
    expect_true(is(structure(pig), "GRanges"))
    expect_equal(piRNAID(pig), "P123")
    expect_true(is(piSequence(pig), "RNAString"))
})

### MicroRNAGene Tests

test_that("MicroRNAGene class functionality", {
    gr <- GRanges(seqnames = Rle("chrX"),
                  ranges = IRanges(start = 1000, end = 1100), strand = Rle("-"))
    micro_gene <- new("MicroRNAGene",
                      geneID = "MG001", symbol = "miR-21", name = "Micro RNA 21",
                      description = "A Micro RNA gene", 
                      structure = gr, microRNAID = "MIR021", 
                      seedSequence = RNAString("CUUACGAUAG")
    )
    
    expect_equal(geneID(micro_gene), "MG001")
    expect_equal(symbol(micro_gene), "miR-21")
    expect_equal(name(micro_gene), "Micro RNA 21")
    expect_equal(description(micro_gene), "A Micro RNA gene")
    expect_true(is(structure(micro_gene), "GRanges"))
    expect_equal(microRNAID(micro_gene), "MIR021")
    expect_true(is(seedSequence(micro_gene), "RNAString"))
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
    
    expect_equal(geneID(trna_gene), "TG001")
    expect_equal(symbol(trna_gene), "tRNA-Gly")
    expect_equal(name(trna_gene), "Transfer RNA Glycine")
    expect_equal(description(trna_gene), "tRNA produce the glycine amino acid")
    expect_true(is(structure(trna_gene), "GRanges"))
    expect_equal(tRNAID(trna_gene), "TRNA001")
    expect_true(is(tRNASequence(trna_gene), "RNAString"))
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
    
    expect_equal(geneID(snrna_gene), "SN001")
    expect_equal(symbol(snrna_gene), "snRNA-U6")
    expect_equal(name(snrna_gene), "Small Nuclear RNA U6")
    expect_equal(description(snrna_gene), "snRNA involved in splicing")
    expect_true(is(structure(snrna_gene), "GRanges"))
    expect_equal(snRNAID(snrna_gene), "SNRNA001")
    expect_true(is(snRNASequence(snrna_gene), "RNAString"))
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
    
    expect_equal(geneID(srp_gene), "SRP001")
    expect_equal(symbol(srp_gene), "SRP RNA")
    expect_equal(name(srp_gene), "Signal Recognition Particle RNA")
    expect_equal(description(srp_gene), "RNA component of the SRP")
    expect_true(is(structure(srp_gene), "GRanges"))
    expect_equal(SRPRNAID(srp_gene), "SRPNA001")
    expect_true(is(SRPSequence(srp_gene), "RNAString"))
})
