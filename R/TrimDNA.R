#' Trims DNA Sequences to the High Quality Region Between Patterns
#' 
#' Aids in trimming DNA sequences to the high quality region between a set of
#' patterns that are potentially present on the left and right sides.
#' 
#' After a sequencing run, it is often necessary to trim the resulting
#' sequences to the high quality region located between a set of patterns.
#' \code{TrimDNA} works as follows: first left and right patterns are
#' identified within the sequences if \code{allowInternal} is \code{TRUE} (the
#' default).  If the patterns are not found internally, then a search is
#' conducted at the flanking ends for patterns that partially overlap the
#' sequence.  The region between the \code{leftPatterns} and
#' \code{rightPatterns} is then returned, unless quality information is
#' provided.  Note that the patterns must be in the same orientation as the
#' sequence, which may require using the \code{reverseComplement} of a PCR
#' primer.
#' 
#' If \code{quality} contains quality scores, these are converted to error
#' probabilities and an exponential moving average is applied to smooth the
#' signal.  The longest region between the \code{leftPatterns} and
#' \code{rightPatterns} where the average error probability is below
#' \code{threshold} is then returned, so long as it has an average error rate
#' of at most \code{maxAverageError}.  Note that it is possible to only filter
#' by \code{maxAverageError} by setting \code{threshold} to \code{1}, or
#' vise-versa by setting \code{maxAverageError} to the same value as
#' \code{threshold}.
#' 
#' @name TrimDNA
#' @param myDNAStringSet A \code{DNAStringSet} or
#' \code{QualityScaledDNAStringSet} object containing the sequences to be
#' trimmed.  If \code{"type"} is \code{"sequences"} then the output class will
#' match the class of \code{myXStringSet}.  Note that the qualities of a
#' \code{QualityScaledDNAStringSet} are ignored because the \code{quality}
#' argument must be supplied separately.
#' @param leftPatterns A \code{DNAStringSet} or character vector of patterns to
#' remove from the left side of \code{myDNAStringSet}, or \code{""} to prevent
#' trimming patterns on the left.
#' @param rightPatterns A \code{DNAStringSet} or character vector of patterns
#' to remove from the right side of \code{myDNAStringSet}, or \code{""} to
#' prevent trimming patterns on the right.
#' @param type Character string indicating the type of results desired.  This
#' should be (an abbreviation of) either \code{"ranges"}, \code{"sequences"} or
#' \code{"both"}.
#' @param quality Either \code{NULL} (the default) to skip quality trimming, or
#' a \code{PhredQuality}, \code{SolexaQuality}, or \code{IlluminaQuality}
#' object containing the quality scores corresponding to \code{myDNAStringSet}.
#' @param maxDistance Numeric between zero and one giving the maximum distance
#' of a match from the \code{leftPatterns} and \code{rightPatterns} to initiate
#' trimming. For example, \code{0.1} (the default) would allow up to 10\%
#' mismatches between a pattern and sequence.
#' @param minOverlap Integer specifying the minimum number of nucleotides the
#' \code{leftPatterns} and \code{rightPatterns} must overlap a sequence to
#' initiate trimming.
#' @param allowInternal Logical initiating whether to search for the
#' \code{leftPatterns} and \code{rightPatterns} within \code{myDNAStringSet},
#' or (\code{FALSE} for) only overlapping the ends.
#' @param alpha Numeric between zero and one giving the smoothing parameter for
#' an exponential moving average that is applied to the quality scores before
#' trimming.  Higher values result in less smoothing than lower values.
#' @param threshold Numeric between zero and one specifying the threshold above
#' which to trim the poor quality regions of the sequence.  Higher values allow
#' more sequence to be preserved at the expense of a greater error rate.
#' @param maxAverageError Numeric between zero and \code{threshold} indicating
#' the maximum average error rate of the trimmed region of the sequence.
#' Trimmed sequences with average error rates above \code{maxAverageError} will
#' be rejected.  Note that the expected number of errors in a sequence is equal
#' to the average error rate multiplied by the length of the sequence.
#' @param maxAmbiguities Numeric between zero and one giving the maximum
#' fraction of ambiguous (e.g., \code{"N"}) positions that are tolerated within
#' the trimmed region of the sequence.  Trimmed sequences with a greater
#' fraction of ambiguities than \code{maxAmbiguities} will be rejected.
#' @param minWidth Integer giving the minimum number of nucleotides a pattern
#' must overlap the sequence to initiate trimming.
#' @param verbose Logical indicating whether to display progress.
#' @return \code{TrimDNA} can return two \code{type}s of results:
#' \code{IRanges} that can be used for trimming \code{myDNAStringSet}, or a
#' trimmed \code{DNAStringSet} or \code{QualityScaledDNAStringSet} containing
#' only those sequences over \code{minWidth} nucleotides after trimming.  Note
#' that ambiguity codes (\code{IUPAC_CODE_MAP}) are supported in the
#' \code{leftPatterns} and \code{rightPatterns}, but not in
#' \code{myDNAStringSet} to prevent trivial matches (e.g., runs of N's).
#' 
#' If \code{type} is \code{"ranges"} (the default) the output is an
#' \code{IRanges} object with the start, end, and width of every sequence.
#' This information can be accessed with the corresponding accessor function
#' (see examples below).  Note that the start will be \code{1} and the end will
#' be \code{0} for sequences that were not at least \code{minWidth} nucleotides
#' after trimming.
#' 
#' If \code{type} is \code{"sequences"} then the trimmed sequences are returned
#' that are at least \code{minWidth} nucleotides in length.
#' 
#' If \code{type} is \code{"both"} the output is a list of two components, the
#' first containing the ranges and the second containing the sequences.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{CorrectFrameshifts}}
#' @examples
#' 
#' # simple example of trimming a single sequence
#' dna <- DNAStringSet("AAAAAAAAAATTACTTCCCCCCCCCC")
#' qscores <- PhredQuality("0000000000AAAAAAAAAAAAAAAA")
#' 
#' x <- TrimDNA(dna,
#'             leftPatterns="AAAAAA",
#'             rightPatterns="CCCCCC",
#'             quality=qscores,
#'             minWidth=1,
#'             allowInternal=TRUE,
#'             type="both")
#' 
#' x[[1]]
#' start(x[[1]])
#' end(x[[1]])
#' width(x[[1]])
#' subseq(dna, start(x[[1]]), end(x[[1]]))
#' x[[2]]
#' 
#' # example of trimming a FASTQ file by quality scores
#' fpath <- system.file("extdata",
#' 	"s_1_sequence.txt",
#' 	package="Biostrings")
#' reads <- readQualityScaledDNAStringSet(fpath)
#' trimmed <- TrimDNA(reads,
#' 	leftPatterns="",
#' 	rightPatterns="",
#' 	type="sequences",
#' 	quality=quality(reads))
#' trimmed
#' DNAStringSet(trimmed) # drop the qualities
#' 
#' @export TrimDNA
TrimDNA <- function(myDNAStringSet,
	leftPatterns,
	rightPatterns,
	type="ranges",
	quality=NULL,
	maxDistance=0.1,
	minOverlap=5,
	allowInternal=TRUE,
	alpha=0.1,
	threshold=0.01,
	maxAverageError=threshold,
	maxAmbiguities=0.1,
	minWidth=36,
	verbose=TRUE) {
	
	# error checking:
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (!(is(quality, "NULL") ||
		is(quality, "PhredQuality") ||
		is(quality, "SolexaQuality") ||
		is(quality, "IlluminaQuality")))
		stop("quality must be a PhredQuality, SolexaQuality, or IlluminaQuality object.")
	if (!is.null(quality) &&
		any(width(myDNAStringSet) != width(quality)))
		stop("The widths of myDNAStringSet must match the widths of quality.")
	if (length(leftPatterns)==0)
		stop("leftPatterns must be at least length one.")
	if (is(leftPatterns, "DNAString") ||
		is(leftPatterns, "DNAStringSet")) {
		leftPatterns <- as.character(leftPatterns)
	} else if (!is.character(leftPatterns)) {
		stop("leftPatterns must be a DNAStringSet or character vector.")
	}
	if (length(rightPatterns)==0)
		stop("rightPatterns must be at least length one.")
	if (is(rightPatterns, "DNAString") ||
		is(rightPatterns, "DNAStringSet")) {
		rightPatterns <- as.character(rightPatterns)
	} else if (!is.character(rightPatterns)) {
		stop("rightPatterns must be a DNAStringSet or character vector.")
	}
	if (length(type) != 1L)
		stop("type must be a character string of length one.")
	TYPES <- c("ranges", "sequences", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (!is.numeric(maxDistance))
		stop("maxDistance must be a numeric.")
	if (maxDistance >= 1 || maxDistance < 0)
		stop("maxDistance must be between zero and one.")
	if (!is.numeric(minOverlap))
		stop("minOverlap must be a numeric.")
	if (minOverlap < 1)
		stop("minOverlap must be at least one.")
	if (!is.logical(allowInternal))
		stop("allowInternal must be a logical.")
	if (!is.numeric(alpha))
		stop("alpha must be a numeric.")
	if (alpha > 1 || alpha <= 0)
		stop("alpha must be between zero and one.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold > 1 || threshold <= 0)
		stop("threshold must be between zero and one.")
	if (!is.numeric(maxAverageError))
		stop("maxAverageError must be a numeric.")
	if (maxAverageError > threshold || maxAverageError <= 0)
		stop("maxAverageError must be between zero and threshold.")
	if (!is.numeric(maxAmbiguities))
		stop("maxAmbiguities must be a numeric.")
	if (maxAmbiguities > 1 || maxAmbiguities < 0)
		stop("maxAmbiguities must be between zero and one.")
	if (!is.numeric(minWidth))
		stop("minWidth must be a numeric.")
	if (minWidth < 0)
		stop("minWidth must be at least zero.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
	if (verbose)
		time.1 <- Sys.time()
	
	.findPatternEnd <- function(pattern) {
		left <- numeric(length(myDNAStringSet))
		l <- nchar(pattern)
		
		if (allowInternal) {
			# find interal matches without indels
			m <- vmatchPattern(pattern,
				myDNAStringSet,
				max.mismatch=floor(l*maxDistance),
				fixed="subject")
			s <- endIndex(m)
			x <- which(lengths(s) > 0)
			if (length(x) > 0) {
				left[x] <- sapply(s[x],
					function(x) {
						if (length(x) > 1) {
							max(x)
						} else {
							x
						}
					})
			}
			
			if (verbose) {
				n <- length(x)
				cat(round(n/length(left)*100, 1),
					"% internal, ",
					sep="")
				flush.console()
			}
			
			x <- which(lengths(s)==0)
		} else {
			x <- seq_len(length(myDNAStringSet))
			n <- 0
		}
		
		if (l > minOverlap) {
			s <- substring(pattern, seq_len(l - minOverlap) + 1)
			for (i in seq_along(s)) {
				if (length(x)==0)
					break
				w <- which(isMatchingAt(s[i],
					.subset(myDNAStringSet, x),
					max.mismatch=floor((l - i + 1)*maxDistance),
					with.indels=TRUE,
					fixed="subject"))
				if (length(w) > 0) {
					left[x[w]] <- l - i
					x <- x[-w]
				}
			}
		}
		
		if (verbose) {
			n <- length(left) - length(x) - n
			cat(round(n/length(left)*100, 1),
				"% flanking\n",
				sep="")
			flush.console()
		}
		
		return(left)
	}
	
	# try to find the flanking left pattern
	lefts <- integer(length(myDNAStringSet))
	lefts[] <- 1L
	for (k in which(nchar(leftPatterns) >= minOverlap)) {
		if (verbose) {
			cat("Finding left pattern",
				ifelse(length(leftPatterns) > 1,
					paste(" (#", k, "): ", sep=""),
					": "),
				sep="")
			flush.console()
		}
		left <- .findPatternEnd(leftPatterns[k]) + 1L
		
		w <- which(left > lefts)
		if (length(w) > 0)
			lefts[w] <- left[w]
	}
	
	# reverse the patterns and subjects
	dna <- myDNAStringSet
	myDNAStringSet <- reverse(myDNAStringSet)
	rightPatterns <- sapply(strsplit(rightPatterns,
			"",
			fixed=TRUE),
		function(x) {
			paste(rev(x), collapse="")
		})
	
	# try to find the flanking right pattern
	rights <- width(myDNAStringSet)
	for (k in which(nchar(rightPatterns) >= minOverlap)) {
		if (verbose) {
			cat(ifelse(k == 1L && any(nchar(leftPatterns) >= minOverlap), "\n", ""),
				"Finding right pattern",
				ifelse(length(rightPatterns) > 1,
					paste(" (#", k, "): ", sep=""),
					": "),
				sep="")
			flush.console()
		}
		right <- width(myDNAStringSet) - .findPatternEnd(rightPatterns[k])
		
		w <- which(right < rights)
		if (length(w) > 0)
			rights[w] <- right[w]
	}
	
	if (!is.null(quality)) {
		if (verbose) {
			cat(ifelse(any(nchar(leftPatterns) >= minOverlap) ||
				any(nchar(rightPatterns) >= minOverlap), "\n", ""),
				"Trimming by quality score: ",
				sep="")
			flush.console()
		}
		bounds <- .Call("movAvg",
			quality,
			match(class(quality)[1],
				c("PhredQuality",
					"SolexaQuality",
					"IlluminaQuality")),
			alpha,
			threshold,
			maxAverageError,
			as.integer(lefts),
			as.integer(rights),
			PACKAGE="DECIPHER")
		if (verbose) {
			wl <- which(lefts != bounds[[1]])
			wr <- which(rights != bounds[[2]])
			cat(100*round(length(wl)/length(lefts), 1),
				"% left, ",
				100*round(length(wr)/length(rights), 1),
				"% right\n",
				sep="")
			flush.console()
		}
		lefts <- bounds[[1]]
		rights <- bounds[[2]]
	}
	
	w <- which((rights - lefts + 1L) < minWidth)
	if (length(w) > 0) {
		lefts[w] <- 1L
		rights[w] <- 0L
	}
	
	ns <- names(myDNAStringSet)
	w <- which(lefts <= rights)
	if (length(w) > 0) {
		myDNAStringSet <- dna[w]
		myDNAStringSet <- subseq(myDNAStringSet,
			lefts[w],
			rights[w])
		
		a <- alphabetFrequency(myDNAStringSet, as.prob=TRUE, baseOnly=TRUE)
		x <- which(a[, "other"] > maxAmbiguities)
		if (length(x) > 0) {
			myDNAStringSet <- myDNAStringSet[-x]
			if (type==1 || type==3) {
				lefts[w[x]] <- 1L
				rights[w[x]] <- 0L
			}
		}
	} else {
		myDNAStringSet <- DNAStringSet()
	}
	
	if (type==1 || type==3)
		ranges <- IRanges(start=lefts,
			end=rights,
			names=ns)
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type==1) {
		return(ranges)
	} else if (type==2) {
		return(myDNAStringSet)
	} else { # type==3
		return(list(ranges,
			myDNAStringSet))
	}
}
