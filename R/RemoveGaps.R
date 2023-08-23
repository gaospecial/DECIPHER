#' Remove Gap Characters in Sequences
#' 
#' Removes gaps ("-" or "." characters) in a set of sequences, either deleting
#' all gaps or only those shared by all sequences in the set.
#' 
#' The \code{removeGaps} argument controls which gaps are removed in
#' \code{myXStringSet}.  Setting \code{removeGaps} to \code{"all"} will remove
#' all gaps in the input sequences, whereas setting \code{removeGaps} to
#' \code{"common"} will remove only gaps that exist in the same position in
#' every sequence.  Therefore, the latter method will leave gaps in place that
#' are not shared by every sequence, requiring that the sequences in
#' \code{myXStringSet} all be the same length (i.e., be aligned).  Setting
#' \code{removeGaps} to \code{"none"} will simply return \code{myXStringSet}
#' unaltered.
#' 
#' @name RemoveGaps
#' @param myXStringSet An \code{AAStringSet}, \code{DNAStringSet}, or
#' \code{RNAStringSet} object containing sequences.
#' @param removeGaps Determines how gaps ("-" or "." characters) are removed in
#' the sequences.  This should be (an unambiguous abbreviation of) one of
#' \code{"none"}, \code{"all"} or \code{"common"}.
#' @param includeMask Logical specifying whether to consider the mask character
#' ("+") as a gap.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @return An \code{XStringSet} of the same type as \code{myXStringSet}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignSeqs}}
#' @examples
#' 
#' dna <- DNAStringSet(c("ACT-G-", "AC--G-"))
#' dna
#' RemoveGaps(dna, "all")
#' RemoveGaps(dna, "common")
#' 
#' @export RemoveGaps
RemoveGaps <- function(myXStringSet,
	removeGaps="all",
	includeMask=FALSE,
	processors=1) {
	
	# error checking
	if (is(myXStringSet, "DNAStringSet")) {
		type <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		type <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		type <- 3L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	if (!is.logical(includeMask))
		stop("includeMask must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (removeGaps==2L) { # all gaps
		ns <- names(myXStringSet)
		myXStringSet <- .Call("removeGaps",
			myXStringSet,
			type,
			includeMask,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	} else if (removeGaps==3L) { # common gaps
		ns <- names(myXStringSet)
		myXStringSet <- .Call("removeCommonGaps",
			myXStringSet,
			type,
			includeMask,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	}
	
	return(myXStringSet)
}
