#' Determine the Number of Terminal Characters
#' 
#' Counts the number of terminal characters for every sequence in an
#' \code{XStringSet}.  Terminal characters are defined as a specific character
#' repeated at the beginning and end of a sequence.
#' 
#' 
#' @name TerminalChar
#' @param myXStringSet An \code{XStringSet} object of sequences.
#' @param char A single character giving the terminal character to count, or an
#' empty character ("") indicating to count both gap ("-") and unknown (".")
#' characters.
#' @return A \code{matrix} containing the results for each sequence in its
#' respective row.  The first column contains the number of leading
#' \code{char}, the second contains the number of trailing \code{char}, and the
#' third contains the total number of characters in-between.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{IdLengths}}
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' dna <- SearchDB(db)
#' t <- TerminalChar(dna)
#' 
#' @export TerminalChar
TerminalChar <- function(myXStringSet,
	char="") {
	
	# error checking
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
		
	
	if (char=="") {
		if (is(myXStringSet, "BStringSet"))
			stop("A single character must be specified with a BStringSet input.")
		gaps <- .Call("gaps",
			myXStringSet,
			ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
			PACKAGE="DECIPHER")
		dimnames(gaps) <- list(names(myXStringSet),
			c("leadingChar","trailingChar","difference"))
	} else {
		numF <- length(myXStringSet)
		maxW <- max(width(myXStringSet))
		
		gaps <- matrix(data=0,
			nrow=numF,
			ncol=3,
			dimnames=list(names(myXStringSet),
				c("leadingChar",
					"trailingChar",
					"difference")))
		myXStringSetRev <- reverse(myXStringSet)
		
		lead <- isMatchingAt(char,
			myXStringSet,
			seq_len(maxW))
		trail <- isMatchingAt(char,
			myXStringSetRev,
			seq_len(maxW))
		
		for (i in seq_len(numF)) {
			# leading characters
			index <- which(!lead[seq_len(width(myXStringSet[i])), i])[1L]
			
			if (is.na(index)) { # sequence is all that character
				gaps[i, 1] <- width(myXStringSet[i])
				gaps[i, 2] <- width(myXStringSet[i])
				gaps[i, 3] <- NA
			} else {
				gaps[i, 1] <- index - 1L
				
				# trailing characters
				index <- which(!trail[, i])[1L]
				gaps[i, 2] <- index - 1L
				
				# difference between leading and trailing gaps
				gaps[i, 3] <- width(myXStringSet[i]) - gaps[i, 2] - gaps[i, 1]
			}
		}
	}
	return(gaps)
}
