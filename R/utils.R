#' bindByQuantiles
#'
#' @param vec 
#' @param q_low 
#' @param q_high 
#'
#' @return ...
#'
#' @export

bindByQuantiles <- function(vec, q_low = 0, q_high = 0.99) {
	mi <- quantile(vec, q_low)
	ma <- quantile(vec[vec != 0], q_high)
	bound_vec <- vec
	bound_vec[bound_vec > ma] <- ma
	bound_vec[bound_vec < mi] <- mi
	return(bound_vec)
}

#' Checking gene names in SCE
#'
#' @param sce 
#' @param genes 
#'
#' @import SingleCellExperiment
#' 
#' @return ...

.checkGenes <- function(sce, genes) {
	absent <- genes[!genes %in% rownames(sce)]
	assertthat::assert_that(
		length(absent) == 0, 
		msg = message(paste(absent, collapse = ', '), ": not in the SCE. Aborting.")
	)
}

#' Checking colData names in SCE
#'
#' @param sce 
#' @param col 
#'
#' @import SingleCellExperiment
#' 
#' @return ...

.checkColData <- function(sce, col) {
	absent <- col[!col %in% colnames(colData(sce))]
	assertthat::assert_that(
		length(absent) == 0, 
		msg = message(paste(absent, collapse = ', '), ": column not in SCE colData. Aborting.")
	)
}

#' Checking embedding SCE
#'
#' @param sce 
#' @param embedding 
#'
#' @import SingleCellExperiment
#' 
#' @return ...

.checkEmbedding <- function(sce, embedding) {
	assertthat::assert_that(
		embedding %in% names(reducedDims(sce)), 
		msg = message(embedding, ": Embedding not in SCE. Aborting.")
	)
}


