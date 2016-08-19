#' Load Kolodziejczyk Counts
#'
#' Load counts of Kolodziejczyk data produced by featureCounts as a scater
#' SCESet object
#'
#' @return SCESet object
#' @examples
#' kolod <- loadKolodCounts()
loadKolodCounts <- function() {
    # Read data
    dir <- "/group/bioi1/shared/public_data/Kolodziejczyk-mESCs-scRNAseq/"
    data <- readr::read_tsv(file.path(dir, "counts/counts.txt"), comment = "#")
    sdrf <- readr::read_tsv(file.path(dir, "metadata/E-MTAB-2600.sdrf.txt"))
    star <- readr::read_csv(file.path(dir, "stats/star_stats.txt"))
    fcount <- readr::read_tsv(file.path(dir, "counts/counts.txt.summary"))
    is.dup <- duplicated(sdrf$"Comment[ENA_RUN]")
    sdrf <- sdrf[!is.dup, ]

    # Extract annotation and counts
    annot <- data[, 1:6]
    colnames(annot)[1] <- "Gene"
    counts <- as.matrix(data[, 7:ncol(data)])
    storage.mode(counts) <- "double"

    # Extract sample names
    samples <- colnames(counts)
    samples <- basename(samples)
    samples <- gsub("\\.Aligned\\.out\\.bam", "", samples)
    colnames(counts) <- samples

    # Remove version numbers from gene names
    genes <- strsplit(annot$Gene, ".", fixed = TRUE)
    genes.ercc <- grep(pattern = "ERCC-", genes, value = TRUE)
    genes <- genes[1:(length(genes) - length(genes.ercc))]
    genes <- matrix(unlist(genes), ncol = 2, byrow = TRUE)
    genes <- c(genes[, 1], genes.ercc)
    rownames(counts) <- genes

    # Set rownames and reorder SDRF
    sdrf <- as.data.frame(sdrf)
    rownames(sdrf) <- sdrf$"Comment[ENA_RUN]"
    sdrf <- sdrf[samples, ]

    # Extract conditions
    conds <- sdrf$"FactorValue[growth condition]"

    # Extract plates
    plates <- sdrf$"Source Name"
    plates <- gsub("_..$","", (gsub("_.$","", plates)))

    # Find bulk samples
    is.bulk <- sdrf$"Material Type" == "cells"

    # Alignment stats
    input <- star$InputNum
    names(input) <- star$Sample
    input <- input[samples]
    mapped <- star$MappedNum
    names(mapped) <- star$Sample
    mapped <- mapped[samples]

    # Counting stats
    fcount <- as.data.frame(fcount)
    colnames(fcount) <- c("Status", samples)
    rownames(fcount) <- fcount$Status
    fcount <- t(fcount[, -1])
    assigned <- fcount[, "Assigned"]
    assigned <- assigned[samples]

    # Phenotype data
    pheno <- data.frame(CellID = samples, Condition = conds, Plate = plates,
                        Input = input, Mapped = mapped,
                        pctMapped = (mapped / input) * 100, Assigned = assigned,
                        pctAssigned = (assigned / input) * 100,
                        pctMappedAssigned = (assigned / mapped) * 100)

    # Produce SCESet object
    kolod <- scater::newSCESet(countData = counts)
    Biobase::fData(kolod) <- as.data.frame(annot)
    Biobase::fData(kolod)$Gene <- genes
    rownames(Biobase::fData(kolod)) <- genes
    Biobase::pData(kolod) <- pheno

    # Remove bulk samples
    kolod <- kolod[, !is.bulk]

    # Calculate QC metrics
    kolod <- scater::calculateQCMetrics(kolod, feature_controls = genes.ercc)

    return(kolod)
}
