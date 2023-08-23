#' Compression/Decompression of Character Vectors
#' 
#' Compresses character vectors into raw vectors, or decompresses raw vectors
#' into character vectors using a variety of codecs.
#' 
#' \code{Codec} can be used to compress/decompress character vectors using
#' different algorithms.  The \code{"nbit"} and \code{"qbit"} methods are
#' tailored specifically to nucleotides and quality scores, respectively.
#' These two methods will store the data as plain text (\code{"ASCII"} format)
#' when it is incompressible.  In such cases, a second \code{compression}
#' method can be given to use in lieu of plain text.  For example
#' \code{compression = c("nbit", "gzip")} will use \code{"gzip"} compression
#' when \code{"nbit"} compression is inappropriate.
#' 
#' When performing the reverse operation, decompression, the type of
#' \code{compression} is automatically detected based on the unique signature
#' ("magic number") added by each compression algorithm.
#' 
#' @name Codec
#' @param x Either a character vector to be compressed, or a list of raw
#' vectors to be decompressed.
#' @param compression The type of compression algorithm to use when \code{x} is
#' a character vector.  This should be (an unambiguous abbreviation of) one of
#' \code{"nbit"} (for nucleotides), \code{"qbit"} (for quality scores),
#' \code{"gzip"}, \code{"bzip2"}, or \code{"xz"}.  If \code{compression} is
#' \code{"nbit"} or \code{"qbit"} then a second method can be provided for
#' cases when \code{x} is incompressible.  Decompression type is determined
#' automatically.  (See details section below.)
#' @param compressRepeats Logical specifying whether to compress exact repeats
#' and reverse complement repeats in a character vector input (\code{x}). Only
#' applicable when \code{compression} is \code{"nbit"}.  Repeat compression in
#' long DNA sequences generally increases compression by about 2\% while
#' requiring three-fold more compression time.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @return If \code{x} is a character vector to be compressed, the output is a
#' list with one element containing a raw vector per character string.  If
#' \code{x} is a list of raw vectors to be decompressed, then the output is a
#' character vector with one string per list element.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @examples
#' 
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- as.character(readDNAStringSet(fas)) # aligned sequences
#' object.size(dna)
#' 
#' # compression
#' system.time(x <- Codec(dna, compression="nbit"))
#' object.size(x)/sum(nchar(dna)) # bytes per position
#' 
#' system.time(g <- Codec(dna, compression="gzip"))
#' object.size(g)/sum(nchar(dna)) # bytes per position
#' 
#' # decompression
#' system.time(y <- Codec(x))
#' stopifnot(dna==y)
#' 
#' system.time(z <- Codec(g))
#' stopifnot(dna==z)
#' 
#' @export Codec
Codec <- function(x,
	compression=c("nbit", "gzip"),
	compressRepeats=FALSE,
	processors=1) {
	
	# error checking
	if (length(compression)==2) {
		if (compression[1] != "nbit" && compression[1] != "qbit")
			stop("The first element of compression must be 'nbit' or 'qbit' when two elements are provided.")
		if (!(compression[2] %in% c("gzip", "bzip2", "xz")))
			stop("The second element of compression must be  'gzip', 'bzip2', or 'xz' when two elements are provided.")
	} else if (length(compression) != 1 || !is.character(compression)) {
		stop("compression must be a character string.")
	} else if (!(compression %in% c("nbit", "qbit", "gzip", "bzip2", "xz"))) {
		stop("Invalid type of compression.")
	}
	if (typeof(x)=="list") {
		if (length(x)==0)
			return(character())
		if (!all(unlist(lapply(x, is.raw))))
			stop("All elements of x must be raw vectors.")
	} else if (!typeof(x)=="character") {
		stop("Invalid type of x.")
		if (length(x)==0)
			stop("Length of x must be greater than zero.")
	}
	if (!is.logical(compressRepeats))
		stop("compressRepeats must be a logical.")
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
	
	if (typeof(x)=="character") {
		if (compression[1]=="nbit") {
			y <- .Call(compression[1],
				x,
				2L - length(compression),
				compressRepeats,
				processors,
				PACKAGE="DECIPHER")
			if (length(compression)==2) {
				w <- which(lengths(y)==0)
				for (i in w)
					y[[i]] <- memCompress(x[i],
						compression[2])
			}
		} else if (compression[1]=="qbit") {
			y <- .Call(compression[1],
				x,
				2L - length(compression),
				processors,
				PACKAGE="DECIPHER")
			if (length(compression)==2) {
				w <- which(lengths(y)==0)
				for (i in w)
					y[[i]] <- memCompress(x[i],
						compression[2])
			}
		} else {
			y <- lapply(x,
				memCompress,
				compression)
		}
	} else { # x is a list
		y <- .Call("decompress",
			x,
			processors,
			PACKAGE="DECIPHER")
		
		bz_header <- as.raw(c(0x42, 0x5a, 0x68))
		xz_header <- as.raw(c(0xfd, 0x37, 0x7a))
		.decompress <- function(x) {
			# choose decompression type
			if (identical(x[1:3], bz_header)) {
				memDecompress(x,
					type="bzip2",
					asChar=TRUE)
			} else if (identical(x[1:3], xz_header)) {
				memDecompress(x,
					type="xz",
					asChar=TRUE)
			} else { # variable header
				memDecompress(x,
					type="gzip",
					asChar=TRUE)
			}
		}
		
		w <- which(unlist(lapply(y, is.na)))
		for (i in w)
			y[i] <- .decompress(x[[i]])
	}
	
	if (!is.null(names(x)))
		names(y) <- names(x)
	
	return(y)
}
