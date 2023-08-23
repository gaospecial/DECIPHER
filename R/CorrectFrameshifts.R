#' Corrects Frameshift Errors In Protein Coding Sequences
#' 
#' Corrects the reading frame to mitigate the impact of frameshift errors
#' caused by insertions or deletions in unaligned nucleotide sequences.
#' 
#' Accurate translation of protein coding sequences can be greatly disrupted by
#' one or two nucleotide phase shifts that occasionally occur during DNA
#' sequencing.  These frameshift errors can potentially be corrected through
#' comparison with other unshifted protein sequences.  This function uses a set
#' of reference amino acid sequences (\code{AAStringSet}) to find and correct
#' frameshift errors in a set of nucleotide sequences (\code{myXStringSet}).
#' First, three frame translation of the nucleotide sequences is performed, and
#' the nearest reference sequence is selected.  Then the optimal reading frame
#' at each position is determined based on a variation of the Guan & Uberbacher
#' (1996) method.  Putative insertions and/or deletions (indels) are returned
#' in the result, typically with close proximity to the true indel locations.
#' For a comparison of this method to others, see Wang et al. (2013).
#' 
#' If \code{type} is \code{"sequences"} or \code{"both"}, then frameshifts are
#' corrected by adding \code{N}'s and/or removing nucleotides.  Note that this
#' changes the nucleotide sequence, and the new sequence often has minor errors
#' because the exact location of the indel(s) cannot be determined.  However,
#' the original frameshifts that disrupted the entire downstream sequence are
#' reduced to local perturbations.  All of the returned nucleotide sequences
#' will have a reading frame starting from the first position.  This allows
#' direct translation, and in practice works well if there is a similar
#' reference \code{myAAStringSet} with the correct reading frame.  Hence it is
#' more important that \code{myAAStringSet} contain a wide variety of sequences
#' than it is that it contain a lot of sequences.
#' 
#' Multiple inputs control the time required for frameshift correction.  The
#' number of sequences in the reference set (\code{myAAStringSet}) will affect
#' the speed of the first search for similar sequences.  Assessing frameshifts
#' in the second step requires order \code{N*M} time, where \code{N} and
#' \code{M} are the lengths of the query (\code{myXStringSet}) and reference
#' sequences.  Two parameters control the number of assessments that are made
#' for each sequence: (1) \code{maxComparisons} determines the maximum number
#' of reference sequences to compare to each query sequence, and (2)
#' \code{acceptDist} defines the maximum distance between a query and reference
#' that is acceptable before continuing to the next query sequence.  A lower
#' value for \code{maxComparisons} or a higher value for \code{acceptDist} will
#' accelerate frameshift correction, potentially at the expense of some
#' accuracy.
#' 
#' @name CorrectFrameshifts
#' @param myXStringSet A \code{DNAStringSet} or \code{RNAStringSet} of
#' unaligned protein coding sequences in 5' to 3' orientation.
#' @param myAAStringSet An \code{AAStringSet} of reference protein sequences.
#' Ideally this would consist of a small set of diverse amino acid sequences
#' belonging to the same group of protein coding sequences as
#' \code{myXStringSet}.
#' @param type Character string indicating the type of result desired.  This
#' should be (an abbreviation of) one of \code{"indels"}, \code{"sequences"},
#' or \code{"both"}.  (See details section below.)
#' @param acceptDistance Numeric giving the maximum distance from a reference
#' sequence that is acceptable to skip any remaining comparisons.
#' @param rejectDistance Numeric giving the maximum distance from a reference
#' sequence that is allowed when correcting frameshifts.  Sequences in
#' \code{myXStringSet} that are greater than \code{rejectDistance} from the
#' nearest reference sequence will only have their length trimmed from the
#' 3'-end to a multiple of three nucleotides without any frameshift correction.
#' @param maxComparisons The number of reference comparisons to make before
#' stopping the search for a closer reference sequence.
#' @param gapOpening Numeric giving the cost for opening a gap between the
#' query and reference sequences.
#' @param gapExtension Numeric giving the cost for extending an open gap
#' between the query and reference sequences.
#' @param frameShift Numeric giving the cost for shifting between frames of the
#' query sequence.
#' @param geneticCode Named character vector in the same format as
#' \code{GENETIC_CODE} (the default), which represents the standard genetic
#' code.
#' @param substitutionMatrix Either a substitution matrix representing the
#' substitution scores for matching two amino acids or the name of the amino
#' acid substitution matrix.  The latter may be one of the following:
#' ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'',
#' ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', ``PFASUM50'' (the
#' default), or ``MIQS''.
#' @param verbose Logical indicating whether to display progress.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @return If \code{type} is \code{"indels"} then the returned object is a list
#' with the same length as \code{myXStringSet}.  Each element is a list with
#' four components: \item{"insertions"}{ Approximate positions of inserted
#' nucleotides, which could be removed to correct the reading frame, or excess
#' nucleotides at the 3'-end that make the length longer than a multiple of
#' three. } \item{"deletions"}{ Approximate positions of deleted nucleotides,
#' which could be added back to correct the reading frame. } \item{"distance"}{
#' The amino acid distance from the nearest reference sequence, between 0 and
#' 1. } \item{"index"}{ The integer index of the reference sequence that was
#' used for frame correction, or \code{0} if no reference sequence was within
#' \code{rejectDistance}. } Note that positions in \code{insertions} and
#' \code{deletions} are sometimes repeated to indicate that the same position
#' needs to be shifted successively more than once to correct the reading
#' frame.
#' 
#' If \code{type} is \code{"sequences"} then the returned object is an
#' \code{XStringSet} of the same type as the input (\code{myXStringSet}).
#' Nucleotides are added or deleted as necessary to correct for frameshifts.
#' The returned sequences all have a reading frame starting from position 1, so
#' that they can be translated directly.
#' 
#' If \code{type} is \code{"both"} then the returned object is a list with two
#' components: one for the \code{"indels"} and the other for the
#' \code{"sequences"}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignTranslation}}, \code{\link{OrientNucleotides}},
#' \code{\link{PFASUM}}
#' @references Guan, X., & Uberbacher, E. C. (1996). Alignments of DNA and
#' protein sequences containing frameshift errors. Computer Applications in the
#' Biosciences : CABIOS, \bold{12(1)}, 31-40.
#' 
#' Wang, Q., et al. (2013). Ecological Patterns of nifH Genes in Four
#' Terrestrial Climatic Zones Explored with Targeted Metagenomics Using
#' FrameBot, a New Informatics Tool. mBio, \bold{4(5)}, e00592-13-e00592-13.
#' @examples
#' 
#' fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # introduce artificial indels
#' n_ins <- 2 # insertions per sequence
#' shifted <- replaceAt(dna,
#' 	lapply(width(dna),
#' 		sample,
#' 		n_ins),
#' 	sample(DNA_BASES,
#' 		n_ins,
#' 		replace=TRUE))
#' n_dels <- 1 # deletions per sequence
#' shifted <- replaceAt(shifted,
#' 	as(lapply(width(shifted),
#' 		function(x) {
#' 			IRanges(sample(x,
#' 					n_dels),
#' 				width=1)
#' 		}), "IRangesList"))
#' 
#' # to make frameshift correction more challenging,
#' # only supply 20 reference amino acid sequences
#' s <- sample(length(dna), 20)
#' x <- CorrectFrameshifts(shifted,
#' 	translate(dna[s]),
#' 	type="both")
#' 
#' # there was a wide range of distances
#' # to the nearest reference sequence
#' quantile(unlist(lapply(x[[1]], `[`, "distance")))
#' 
#' # none of the sequences were > rejectDistance
#' # from the nearest reference sequence
#' length(which(unlist(lapply(x[[1]], `[`, "index"))==0))
#' 
#' # the number of indels was generally correct
#' table(unlist(lapply(x[[1]], function(x) {
#' 	length(x$insertions)})))/length(shifted)
#' table(unlist(lapply(x[[1]], function(x) {
#' 	length(x$deletions)})))/length(shifted)
#' 
#' # align and display the translations
#' AA <- AlignTranslation(x$sequences,
#' 	readingFrame=1,
#' 	type="AAStringSet")
#' BrowseSeqs(AA)
#' 
#' @export CorrectFrameshifts
CorrectFrameshifts <- function(myXStringSet,
	myAAStringSet,
	type="indels",
	acceptDistance=0.01,
	rejectDistance=0.60,
	maxComparisons=10,
	gapOpening=-13,
	gapExtension=-1,
	frameShift=-15,
	geneticCode=GENETIC_CODE,
	substitutionMatrix="PFASUM50",
	verbose=TRUE,
	processors=1) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (!is(myAAStringSet, "AAStringSet"))
		stop("myAAStringSet must be an AAStringSet.")
	if (length(myXStringSet)==0)
		stop("At least one sequence is required in myXStringSet.")
	if (length(myAAStringSet)==0)
		stop("At least one sequence is required in myAAStringSet.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("-", myAAStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myAAStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myAAStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myAAStringSet must be removed before correcting frameshifts.")
	org_index <- which(!duplicated(myAAStringSet))
	myAAStringSet <- myAAStringSet[org_index]
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
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	TYPES <- c("indels", "sequences", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(frameShift))
		stop("frameShift must be a numeric.")
	if (!is.numeric(rejectDistance))
		stop("rejectDistance must be a numeric.")
	if (rejectDistance > 1)
		stop("rejectDistance can be at most 1.")
	if (rejectDistance <= 0)
		stop("rejectDistance must be greater than zero.")
	if (!is.numeric(acceptDistance))
		stop("acceptDistance must be a numeric.")
	if (acceptDistance > rejectDistance)
		stop("acceptDistance can be at most rejectDistance.")
	if (acceptDistance < 0)
		stop("acceptDistance must be at least zero.")
	if (!is.numeric(maxComparisons))
		stop("maxComparisons must be a numeric.")
	if (maxComparisons < 1)
		stop("maxComparisons must be at least one.")
	if (floor(maxComparisons)!=maxComparisons)
		stop("maxComparisons must be a whole number.")
	AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
	if (is.character(substitutionMatrix)) {
		if (substitutionMatrix=="PFASUM50") {
			subMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")) {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment(), package=ifelse(substitutionMatrix=="MIQS", "DECIPHER", "Biostrings"))))
		} else {
			stop("Invalid substitutionMatrix.")
		}
	} else if (is.matrix(substitutionMatrix)) {
		if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
			any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
			stop("substitutionMatrix is incomplete.")
		subMatrix <- substitutionMatrix
	} else {
		stop("Invalid substitutionMatrix.")
	}
	subMatrix <- subMatrix[AAs, AAs] + 0 # convert to numeric
	
	# de-replicate
	w <- which(!duplicated(myXStringSet))
	l <- length(myXStringSet)
	ns <- names(myXStringSet)
	if (length(w) < l) {
		dupes <- TRUE
		derep <- numeric(l)
		derep[w] <- seq_along(w)
		derep[-w] <- match(myXStringSet[-w],
			myXStringSet[w])
		myXStringSet <- myXStringSet[w]
		l <- length(w)
	} else {
		dupes <- FALSE
	}
	
	ends <- width(myXStringSet)
	if (min(ends) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	
	frames <- AAStringSet()
	for (i in 1:3) {
		end <- ends
		offset <- end - i + 1
		end <- end - offset %% 3
		end <- ifelse(end < i - 1,
			i - 1,
			end)
		
		AA <- translate(subseq(myXStringSet,
				i,
				end),
			genetic.code=geneticCode,
			if.fuzzy.codon="solve")
		frames <- c(frames,
			AA)
	}
	
	if (length(myAAStringSet) > 1) {
		if (verbose) {
			time.1 <- Sys.time()
			cat("Finding the closest reference amino acid sequences:\n",
				sep="")
			flush.console()
			pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
		} else {
			pBar <- NULL
		}
		
		v1 <- .Call("enumerateSequenceAA",
			frames,
			7L,
			PACKAGE="DECIPHER")
		v2 <- .Call("enumerateSequenceAA",
			myAAStringSet,
			7L,
			PACKAGE="DECIPHER")
		
		v1 <- lapply(v1,
			sort)
		v2 <- lapply(v2,
			sort)
		
		d <- .Call("matchListsDual",
			v1,
			v2,
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		d <- rowsum(d,
			rep(seq_len(l),
				3),
			na.rm=TRUE)
		
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		}
	} else {
		d <- matrix(0, nrow=l)
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Assessing frameshifts in nucleotide sequences:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
	}
	
	widths <- width(myAAStringSet)
	
	index <- matrix(nrow=nrow(d), ncol=ncol(d))
	for (i in seq_len(nrow(d)))
		index[i,] <- order(d[i,], widths, decreasing=TRUE)
	if (ncol(index) > maxComparisons) {
		index <- index[, seq_len(maxComparisons), drop=FALSE]
	} else {
		maxComparisons <- ncol(index)
	}
	
	X <- .Call("findFrameshifts",
		myAAStringSet,
		ends,
		frames,
		index,
		org_index,
		maxComparisons,
		gapOpening,
		gapExtension,
		frameShift,
		acceptDistance,
		rejectDistance,
		subMatrix,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	if (type != 2) {
		X <- lapply(X,
			setNames,
			c("insertions", "deletions", "distance", "index"))
	}
	
	if (type > 1) {
		pos <- lapply(X,
			function(x) {
				ins <- rle(x[[1]])
				dels <- rle(x[[2]])
				IRanges(start=c(ins[[2]],
						dels[[2]]),
					width=c(ins[[1]],
						rep(0,
							length(dels[[2]]))))
			})
		pos <- unname(pos)
		
		val <- lapply(X,
			function(x) {
				ins <- rep("",
					length(unique(x[[1]])))
				dels <- rle(x[[2]])
				types <- c("N", "NN")
				dels <- types[dels[[1]]]
				c(ins, dels)
			})
		
		if (type==2L) { # sequences
			X <- replaceAt(myXStringSet,
				as(pos, "IRangesList"),
				val)
		} else { # both
			X <- list(indels=X,
				sequences=replaceAt(myXStringSet,
					as(pos, "IRangesList"),
					val))
		}
	}
	
	if (dupes) { # re-replicate
		if (type==1 || type==2) {
			X <- X[derep]
		} else {
			X[[1]] <- X[[1]][derep]
			X[[2]] <- X[[2]][derep]
		}
	}
	
	if (type==1 || type==2) {
		names(X) <- ns
	} else {
		names(X[[1]]) <- names(X[[2]]) <- ns
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
	}
	
	return(X)
}
