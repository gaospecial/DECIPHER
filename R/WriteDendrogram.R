#' Write a Dendrogram to Newick Format
#' 
#' Writes a dendrogram object to a file in Newick (also known as New Hampshire)
#' parenthetic format.
#' 
#' \code{WriteDendrogram} will write a dendrogram object to a \code{file} in
#' standard Newick format.  Note that special characters (commas, square
#' brackets, colons, semi-colons, and parentheses) present in leaf labels will
#' likely cause a broken Newick file unless \code{quoteLabels} is \code{TRUE}
#' (the default).
#' 
#' @name WriteDendrogram
#' @param x An object of class \code{dendrogram}.
#' @param file A connection or a character string naming the file path where
#' the tree should be exported.  If "" (the default), the tree is printed to
#' the standard output connection, the console unless redirected by sink.
#' @param quoteLabels Logical specifying whether to place leaf labels in double
#' quotes.
#' @param convertBlanks Logical specifying whether to convert spaces in leaf
#' labels to underscores.
#' @param internalLabels Logical indicating whether to write any ``edgetext''
#' preceding a node as an internal node label.
#' @param digits The maximum number of digits to print for edge lengths.
#' @param append Logical indicating whether to append to an existing
#' \code{file}.  Only applicable if \code{file} is a character string.  If
#' \code{FALSE} (the default), then the file is overwritten.
#' @return \code{NULL}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{TreeLine}}, \code{\link{ReadDendrogram}}
#' @examples
#' 
#' dists <- matrix(c(0, 10, 20, 10, 0, 5, 20, 5, 0),
#'     nrow=3,
#'     dimnames=list(c("dog", "elephant", "horse")))
#' dend <- TreeLine(myDistMatrix=dists, method="NJ")
#' WriteDendrogram(dend)
#' 
#' @export WriteDendrogram
WriteDendrogram <- function(x,
	file="",
	quoteLabels=TRUE,
	convertBlanks=!quoteLabels,
	internalLabels=TRUE,
	digits=10,
	append=FALSE) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x is not a dendrogram.")
	if (!is.logical(quoteLabels))
		stop("quoteLabels must be a logical.")
	if (!is.logical(convertBlanks))
		stop("convertBlanks must be a logical.")
	if (!is.numeric(digits))
		stop("digits must be a numeric.")
	if (floor(digits)!=digits)
		stop("digits must be a whole number.")
	if (digits < 1)
		stop("digits must be at least 1.")
	if (!is.logical(append))
		stop("append must be a logical.")
	if (!is.logical(internalLabels))
		stop("internalLabels must be a logical.")
	
	if (is.character(file)) {
		if (file == "") {
			file <- stdout()
		} else if (substring(file, 1L, 1L) == "|") {
			file <- pipe(substring(file, 2L), "w")
			on.exit(close(file))
		} else {
			file <- file(file, "w")
			on.exit(close(file))
		}
	}
	
	getLab <- function(LAB) {
		if (is.null(LAB))
			return("")
		lab <- gsub("'", "''", LAB, fixed=TRUE)
		if (convertBlanks)
			lab <- gsub(" ", "_", lab, fixed=TRUE)
		if (quoteLabels)
			lab <- paste('"', lab, '"', sep="")
		return(lab)
	}
	
	.dendrogram2newick <- function(x, height=attr(x, "height"), root=TRUE) {
		if (is.leaf(x)) {
			cat(getLab(attr(x, "label")),
				":",
				round(height - attr(x, "height"),
					digits=digits),
				sep="",
				file=file,
				append=TRUE)
		} else {
			cat("(",
				file=file,
				append=TRUE)
			for (i in seq_along(x)) {
				.dendrogram2newick(x[[i]],
					attr(x, "height"),
					root=FALSE)
				if (i < length(x))
					cat(",",
						file=file,
						append=TRUE)
			}
			if (root) {
				cat(");\n",
					file=file,
					append=TRUE)
			} else {
				cat(")",
					ifelse(internalLabels,
						getLab(attr(x, "edgetext")),
						""),
					":",
					round(height - attr(x, "height"),
						digits=digits),
					sep="",
					file=file,
					append=TRUE)
			}
		}
	}
	
	if (!append) # overwrite the file
		cat("", file=file)
	.dendrogram2newick(x)
	invisible(NULL)
}
