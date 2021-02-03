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

