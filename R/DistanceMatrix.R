#' Calculate the Distances Between Sequences
#' 
#' Calculates a distance matrix for an \code{XStringSet}.  Each element of the
#' distance matrix corresponds to the dissimilarity between two sequences in
#' the \code{XStringSet}.
#' 
#' The uncorrected (\code{correction = "none"}) distance matrix represents the
#' hamming distance between each of the sequences in \code{myXStringSet}.
#' Ambiguity can be represented using the characters of the
#' \code{IUPAC_CODE_MAP} for \code{DNAStringSet} and \code{RNAStringSet}
#' inputs, or using the \code{AMINO_ACID_CODE} for an \code{AAStringSet} input.
#' For example, the distance between an 'N' and any other nucleotide base is
#' zero.  The letters B (N or D), J (I or L), Z (Q or E), and X (any letter)
#' are degenerate in the \code{AMINO_ACID_CODE}.
#' 
#' If \code{includeTerminalGaps = FALSE} then terminal gaps ("-" or "."
#' characters) are not included in sequence length.  This can be faster since
#' only the positions common to each pair of sequences are compared.  Sequences
#' with no overlapping region in the alignment are given a value of \code{NA},
#' unless \code{includeTerminalGaps = TRUE}, in which case distance is 100\%.
#' Masked characters (\code{"+"}) in either sequence are not considered in
#' distance.
#' 
#' Penalizing gap-to-letter mismatches specifies whether to penalize these
#' special mismatch types and include them in the total length when calculating
#' distance.  Both "-" and "." characters are interpreted as gaps.  The default
#' behavior is to calculate distance as the fraction of positions that differ
#' across the region of the alignment shared by both sequences (not including
#' gap-to-gap matches).
#' 
#' Two correction factors are available, \code{"JC69"} and \code{"F81"}, which
#' are described in \code{\link{MODELS}}.  Both transform raw distance
#' (\eqn{d}) by \deqn{-E * \log \left( 1 - d / E \right)}, where \deqn{E = 1 -
#' \sum_{i \in sym} f_i^2} and (\eqn{f}) is the relative frequency of each
#' symbol (\eqn{sym}).  In the \code{"JC69"} model symbols are assumed to have
#' equal frequency, whereas in the \code{"F81"} model the symbol frequencies
#' are empirically derived from the input \code{myXStringSet}.  Note that gaps
#' are treated as an additional symbol when \code{penalizeGapLetterMatches} is
#' \code{TRUE}.
#' 
#' The elements of the distance matrix can be referenced by \code{dimnames}
#' corresponding to the \code{names} of the \code{XStringSet}.  Additionally,
#' an attribute named "correction" specifying the method of correction used can
#' be accessed using the function \code{attr}.
#' 
#' @name DistanceMatrix
#' @param myXStringSet A \code{DNAStringSet}, \code{RNAStringSet}, or
#' \code{AAStringSet} object of aligned sequences.
#' @param type Character string indicating the type of output desired.  This
#' should be either \code{"matrix"} or \code{"dist"}.  (See value section
#' below.)
#' @param method Character string determining the region in which distance is
#' calculated.  This should be (an unambiguous abbreviation of) one of
#' \code{"overlap"}, \code{"shortest"}, or \code{"longest"}.  The default
#' \code{method} (\code{"overlap"}) calculates distance in the overlapping
#' region between terminal gaps when \code{includeTerminalGaps} is \code{FALSE}
#' and the entire alignment otherwise.  Setting \code{method} to
#' \code{"shortest"} or \code{"longest"} will use the region between the start
#' and end of the shortest or longest sequence, respectively, for each pairwise
#' distance.  The \code{method} is only relevant when
#' \code{includeTerminalGaps} is \code{TRUE}.
#' @param includeTerminalGaps Logical specifying whether or not to include
#' terminal gaps ("-" or "." characters on each end of the sequence) into the
#' calculation of distance.
#' @param penalizeGapLetterMatches Logical specifying whether or not to
#' consider gap-to-letter matches as mismatches.  If \code{FALSE}, then
#' gap-to-letter matches are not included in the total length used to calculate
#' distance.  The default (\code{TRUE}) is to penalize gap-to-letter mismatches
#' the same as letter-to-letter mismatches.  If \code{NA} then gap-to-letter
#' mismatches are only penalized once per insertion or deletion, i.e., when
#' changing to gap-to-letter or letter-to-gap.
#' @param minCoverage Numeric between zero and one giving the minimum fraction
#' of sequence positions (not gap or mask) in the shortest sequence that must
#' be overlapping with the longer sequence in each pair.  Sequences below
#' \code{minCoverage} will be assigned \code{NA} distances.  Note that
#' non-overlapping sequences are always given \code{NA} distances, regardless
#' of \code{minCoverage}, unless \code{includeTerminalGaps} is \code{TRUE}
#' (distance = 100\%).
#' @param correction The substitution model used for distance correction.  This
#' should be (an abbreviation of) either \code{"none"}, \code{"Jukes-Cantor"}
#' (i.e., \code{"JC69"}), or \code{"F81"}.  For \code{"F81"} letter frequencies
#' are derived from \code{myXStringSet}.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @return If \code{type} is \code{"matrix"}, a symmetric matrix where each
#' element is the distance between the sequences referenced by the respective
#' row and column.  The \code{dimnames} of the matrix correspond to the
#' \code{names} of the \code{XStringSet}.
#' 
#' If \code{type} is \code{"dist"}, an object of \code{class} \code{"dist"}
#' that contains one triangle of the distance matrix as a vector.  Since the
#' distance matrix is symmetric, storing only one triangle is more memory
#' efficient.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{TreeLine}}
#' @examples
#' 
#' # example of using the defaults:
#' dna <- DNAStringSet(c("ACTG", "ACCG"))
#' dna
#' DistanceMatrix(dna)
#' 
#' # changing the output type to "dist":
#' d <- DistanceMatrix(dna, type="dist")
#' d
#' length(d) # minimal memory space required
#' m <- as.matrix(d)
#' length(m) # more memory space required
#' 
#' # supplying an AAStringSet
#' aa <- AAStringSet(c("ASYK", "ATYK", "CTWN"))
#' aa
#' DistanceMatrix(aa)
#' 
#' # defaults compare intersection of internal ranges:
#' dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
#' dna
#' d <- DistanceMatrix(dna)
#' d
#' # d[1,2] is 1 base in 4 = 0.25
#' 
#' # compare union of internal positions, without terminal gaps:
#' dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
#' dna
#' d <- DistanceMatrix(dna, includeTerminalGaps=TRUE)
#' d
#' # d[1,2] is now 2 bases in 5 = 0.40
#' 
#' # gap ("-") and unknown (".") characters are interchangeable:
#' dna <- DNAStringSet(c("ANGCT.", ".ACCT-"))
#' dna
#' d <- DistanceMatrix(dna, includeTerminalGaps=TRUE)
#' d
#' # d[1,2] is still 2 bases in 5 = 0.40
#' 
#' # compare different methods for calculating distance:
#' dna <- DNAStringSet(c("--ACTG", "TGAGT-"))
#' dna
#' DistanceMatrix(dna, method="overlap") # 1/3
#' DistanceMatrix(dna, method="overlap",
#'                minCoverage=1) # NA (only 3/4 overlapping)
#' DistanceMatrix(dna, method="shortest",
#'                includeTerminalGaps=FALSE) # 1/3
#' DistanceMatrix(dna, method="shortest",
#'                includeTerminalGaps=TRUE) # 2/4
#' DistanceMatrix(dna, method="shortest",
#'                includeTerminalGaps=TRUE,
#'                penalizeGapLetterMatches=FALSE) # 1/3
#' DistanceMatrix(dna, method="longest",
#'                includeTerminalGaps=FALSE) # 1/3
#' DistanceMatrix(dna, method="longest",
#'                includeTerminalGaps=TRUE) # 3/5
#' DistanceMatrix(dna, method="longest",
#'                includeTerminalGaps=TRUE,
#'                penalizeGapLetterMatches=FALSE) # 1/3
#' 
#' # neither internal nor external gap/gap matches are considered:
#' dna <- DNAStringSet(c("--A-CTA", "-AG-C--"))
#' dna
#' DistanceMatrix(dna) # 1/2
#' DistanceMatrix(dna, includeTerminalGaps=TRUE) # 4/5
#' DistanceMatrix(dna, includeTerminalGaps=TRUE,
#'                method="shortest") # 2/3
#' DistanceMatrix(dna, includeTerminalGaps=TRUE,
#'                method="longest") # 3/4
#' 
#' @export DistanceMatrix
DistanceMatrix <- function(myXStringSet,
	method="overlap",
	type="matrix",
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=TRUE,
	minCoverage=0,
	correction="none",
	processors=1,
	verbose=TRUE) {
	
	# error checking
	TYPES <- c("matrix", "dist")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	METHODS <- c("overlap", "shortest", "longest")
	method <- pmatch(method[1], METHODS)
	if (is.na(method))
		stop("Invalid method.")
	if (method==-1)
		stop("Ambiguous method.")
	CORRECTIONS <- c("none", "Jukes-Cantor", "JC", "JC69", "F81")
	correction <- pmatch(correction, CORRECTIONS)
	if (is.na(correction))
		stop("Invalid distance correction method.")
	if (correction==-1)
		stop("Ambiguous distance correction method.")
	if (correction==3 || correction==4)
		correction <- 2
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	if (!is.numeric(minCoverage))
		stop("maxCoverage must be a numeric.")
	if (minCoverage < 0)
		stop("minCoverage must be at least zero.")
	if (minCoverage > 1)
		stop("minCoverage can be at most one.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
	if (is(myXStringSet, "BStringSet"))
		stop("myXStringSet cannot be a BStringSet.")
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
	
	maxW <- unique(width(myXStringSet))
	if (length(maxW)!=1) {
		warning("\n",
			length(maxW),
			" different sequence lengths.\n",
			"Using shorter length in each comparison.\n")
	}
	numF <- length(myXStringSet)
	if (numF < 2) {
		stop("At least two sequences are required.")
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	# calculate distance correction
	if (correction != 1L) {
		if (is(myXStringSet, "DNAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(DNA_BASES, "-")
			} else {
				alphabet <- DNA_BASES
			}
		} else if (is(myXStringSet, "RNAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(RNA_BASES, "-")
			} else {
				alphabet <- RNA_BASES
			}
		} else if (is(myXStringSet, "AAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(AA_STANDARD, "-")
			} else {
				alphabet <- AA_STANDARD
			}
		}
		if (correction == 5L) {
			E <- letterFrequency(myXStringSet,
				alphabet)
			E <- colSums(E)
			E <- E/sum(E)
			E <- 1 - sum(E^2)
		} else { # JC69
			E <- 1/length(alphabet) # assumes even frequencies
			E <- 1 - length(alphabet)*sum(E^2)
		}
	} else {
		E <- 0
	}
	
	# calculate the distance matrix
	distMatrix <- .Call("distMatrix",
		myXStringSet,
		ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
		includeTerminalGaps,
		penalizeGapLetterMatches,
		TRUE, # full matrix
		type,
		E,
		minCoverage,
		method,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	if (type==1) { # matrix
		dimnames(distMatrix) <- list(names(myXStringSet),
			names(myXStringSet))
	} else { # dist
		attr(distMatrix, "Size") <- length(myXStringSet)
		if (!is.null(names(myXStringSet)))
			attr(distMatrix, "Labels") <- names(myXStringSet)
		attr(distMatrix, "Diag") <- TRUE
		attr(distMatrix, "Upper") <- TRUE
		class(distMatrix) <- "dist"
	}
	attr(distMatrix, "correction") <- CORRECTIONS[correction]
	
	if (verbose) {
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(distMatrix)
}
