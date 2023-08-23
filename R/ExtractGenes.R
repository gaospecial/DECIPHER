#' Extract Predicted Genes from a Genome
#' 
#' Extracts predicted genes from the genome used for prediction.
#' 
#' Extracts a set of gene predictions as either DNA, mRNA, or proteins.
#' 
#' @name ExtractGenes
#' @param x An object of class \code{Genes}.
#' @param myDNAStringSet The \code{DNAStringSet} object used in generating
#' \code{x}.
#' @param type The class of sequences to return.  This should be (an
#' unambiguous abbreviation of) one of \code{"AAStringSet"},
#' \code{"DNAStringSet"} (the default), or \code{"RNAStringSet"}.
#' @param \dots Other parameters passed directly to \code{translate}.
#' @return An \code{"AAStringSet"}, \code{"DNAStringSet"}, or
#' \code{"RNAStringSet"} determined by \code{type}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{FindGenes}}, \code{\link{Genes-class}},
#' \code{\link{WriteGenes}}
#' @examples
#' 
#' # import a test genome
#' fas <- system.file("extdata",
#' 	"Chlamydia_trachomatis_NC_000117.fas.gz",
#' 	package="DECIPHER")
#' genome <- readDNAStringSet(fas)
#' 
#' x <- FindGenes(genome)
#' genes <- ExtractGenes(x, genome)
#' proteins <- ExtractGenes(x, genome, type="AAStringSet")
#' 
#' @export ExtractGenes
ExtractGenes <- function(x,
	myDNAStringSet,
	type="DNAStringSet",
	...) {
	
	if (!is(x, "Genes"))
		stop("x must be an object of class 'Genes'.")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet")
	if (length(type)==0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (type==3 && any(1/x[, "Gene"] < 0))
		stop("Only open reading frames in x can be translated when type is 'AAStringSet'.")
	if (!all(attr(x, "widths")==width(myDNAStringSet)))
		stop("Mismatch between x and myDNAStringSet.")
	
	y <- rep(list(IRanges()),
		length(myDNAStringSet))
	t <- tapply(seq_len(nrow(x)),
		x[, "Index"],
		c)
	for (i in seq_along(t)) {
		y[[as.numeric(names(t)[i])]] <- IRanges(x[t[[i]], "Begin"],
			x[t[[i]], "End"])
	}
	y <- IRangesList(y)
	z <- unlist(extractAt(myDNAStringSet, y))
	w <-  which(x[, "Strand"]==1)
	if (length(w) > 0)
		z[w] <- reverseComplement(z[w])
	
	if (type==2L) {
		z <- RNAStringSet(z)
	} else if (type==3L) {
		z <- translate(z,
			genetic.code=attr(x, "geneticCode"),
			...)
	}
	
	names(z) <- NULL
	
	return(z)
}
