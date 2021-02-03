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
#' @return ...

.checkGenes <- function(sce, genes) {
	absent <- genes[!genes %in% rownames(sce)]
	assertthat::assert_that(
		length(absent) == 0, 
		msg = message(paste(absent, collapse = ', '), ": not in the SCE. Aborting.")
	)
}


