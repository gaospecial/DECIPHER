RemoveGaps <- function(myXStringSet,
	removeGaps="all",
	processors=1) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
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
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	} else if (removeGaps==3L) { # common gaps
		ns <- names(myXStringSet)
		myXStringSet <- .Call("removeCommonGaps",
			myXStringSet,
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	}
	
	return(myXStringSet)
}
