#' Expand Ambiguities into All Permutations of a DNAStringSet
#' 
#' Performs the inverse function of \code{ConsensusSequence} by expanding any
#' ambiguities present in sequences.
#' 
#' Ambiguity codes in the \code{IUPAC_CODE_MAP} can be used to represent
#' multiple nucleotides at a single position.  Using these letters, multiple
#' oligonucleotide permutations can be represented with a single ambiguous
#' sequence.  This function expands each sequence in the \code{DNAStringSet}
#' input into all of its permutations.  Note that sequences with many
#' ambiguities can result in a very large number of potential permutations.
#' 
#' @name Disambiguate
#' @param myXStringSet A \code{DNAStringSet} or \code{RNAStringSet} object of
#' sequences.
#' @return A \code{DNAStringSetList} or \code{RNAStringSetList} with one
#' element for each sequence in \code{myXStringSet}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{ConsensusSequence}}
#' @examples
#' 
#' dna <- DNAStringSet(c("ACST", "NNN"))
#' dna_list <- Disambiguate(dna)
#' dna_list[[1]]
#' dna_list[[2]]
#' unlist(dna_list)
#' 
#' rna <- RNAStringSet(c("ACGU", "AGAU")) # 2 permutations
#' rna <- ConsensusSequence(rna) # "ASRU"
#' Disambiguate(rna) # 4 permutations
#' 
#' @export Disambiguate
Disambiguate <- function(myXStringSet) {
	
	# error checking
	if (is(myXStringSet, "DNAStringSet")) {
		type <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		type <- 2L
	} else {
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	}
	
	myXStringSet <- .Call("expandAmbiguities",
		myXStringSet,
		ifelse(type==1,
			"T",
			"U"),
		PACKAGE="DECIPHER")
	
	if (type==1) {
		myXStringSetList <- DNAStringSetList(myXStringSet)
	} else {
		myXStringSetList <- RNAStringSetList(myXStringSet)
	}
	
	return(myXStringSetList)
}
