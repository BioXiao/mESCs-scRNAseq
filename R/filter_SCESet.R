#' Filter Kolodziejczyk counts
#'
#' Take a SCESet object containing the Kolodziejczyk mESC counts and filter it
#' as described in `Filtering-mESCs-01.Rmd`.
#'
#' @param kolod SCESet object containing counts
#' @param genes Logical. Whether to filter high dropout genes.
#' @param cells Logical. Whether to filter bad cells based on dropout and
#'        library size.
#' @param plates Logical. Whether to filter known bad plates.
#' @param ercc Logical. Whether to remove ERCC genes.
#' @param condition String. Condition to keep. Options are: "all",
#'        "alternative 2i media + LIF", "serum + LIF" or
#'        "standard 2i media + LIF".
#'
#' @return Filtered SCSSet object.
#' @examples
#' kolod <- loadKolodCounts()
#' kolod.filter <- filterKolodCounts(kolod)
filterKolodCounts <- function(kolod, genes = TRUE, cells = TRUE,
                              plates = TRUE, ercc = FALSE,
                              condition = c("all", "alternative 2i media + LIF",
                                            "serum + LIF",
                                            "standard 2i media + LIF")) {

    condition <- match.arg(condition)

    gene.names <- Biobase::featureNames(kolod)
    genes.ercc <- grep(pattern = "ERCC-", gene.names, value = TRUE)

    # Keep cells with low dropout and decent library size
    if (cells) {
        # Select only endogenous genes
        kolod.end <- kolod[!Biobase::fData(kolod)$is_feature_control, ]
        kolod.end <- scater::calculateQCMetrics(kolod.end)
        # Keep cells with dropout < 80%
        low.drop <- Biobase::pData(kolod.end)$pct_dropout < 80
        # Keep cells with sqrt(lib.size) > 750
        hi.lib <- sqrt(Biobase::pData(kolod.end)$total_counts) > 750
        # Filter SCESet
        kolod <- kolod[, low.drop & hi.lib]
    }

    # Remove known bad plates
    if (plates) {
        bad.plates <- c("2i_3","a2i_3","serum_3")
        is.bad <- Biobase::pData(kolod)$Plate %in% bad.plates
        kolod <- kolod[, !is.bad]
    }

    # Remove genes with dropout > 90%
    if (genes) {
        kolod <- scater::calculateQCMetrics(kolod)
        lo.drop <- Biobase::fData(kolod)$pct_dropout <= 90
        kolod <- kolod[lo.drop, ]
    }

    # Select condition
    if (condition != "all") {
        kolod <- kolod[, Biobase::pData(kolod)$Condition == condition]
    }

    # Remove spike-ins
    if (ercc) {
        kolod <- kolod[!(Biobase::fData(kolod)$Gene %in% genes.ercc), ]
    }

    # Make sure metrics are correct
    kolod <- scater::calculateQCMetrics(kolod, feature_controls = genes.ercc)

    return(kolod)
}
