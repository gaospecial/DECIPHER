#' Design a Set of DNA Microarray Probes for Detecting Sequences
#' 
#' Chooses the set of microarray probes maximizing sensitivity and specificity
#' to each target consensus sequence.
#' 
#' The algorithm begins by determining the optimal length of probes required to
#' meet the input constraints while maximizing sensitivity to the target
#' consensus sequence at the specified hybridization formamide concentration.
#' This set of potential target sites is then scored based on the possibility
#' of cross-hybridizing to the other non-target sequences.  The set of probes
#' is returned with the minimum possibility of cross-hybridizing.
#' 
#' @name DesignArray
#' @param myDNAStringSet A \code{DNAStringSet} object of aligned consensus
#' sequences.
#' @param maxProbeLength The maximum length of probes, not including the poly-T
#' spacer.  Ideally less than 27 nucleotides.
#' @param minProbeLength The minimum length of probes, not including the poly-T
#' spacer.  Ideally more than 18 nucleotides.
#' @param maxPermutations The maximum number of probe permutations required to
#' represent a target site.  For example, if a target site has an 'N' then 4
#' probes are required because probes cannot be ambiguous.  Typically fewer
#' permutations are preferably because this requires less space on the
#' microarray and simplifies interpretation of the results.
#' @param numRecordedMismatches The maximum number of recorded potential
#' cross-hybridizations for any target site.
#' @param numProbes The target number of probes on the microarray per input
#' consensus sequence.
#' @param start Integer specifying the starting position in the alignment where
#' potential forward primer target sites begin.  Preferably a position that is
#' included in most sequences in the alignment.
#' @param end Integer specifying the ending position in the alignment where
#' potential reverse primer target sites end.  Preferably a position that is
#' included in most sequences in the alignment.
#' @param maxOverlap Maximum overlap in nucleotides between target sites on the
#' sequence.
#' @param hybridizationFormamide The formamide concentration (\%, vol/vol) used
#' in hybridization at 42 degrees Celsius.  Note that this concentration is
#' used to approximate hybridization efficiency of cross-amplifications.
#' @param minMeltingFormamide The minimum melting point formamide concentration
#' (\%, vol/vol) of the designed probes.  The melting point is defined as the
#' concentration where half of the template is bound to probe.
#' @param maxMeltingFormamide The maximum melting point formamide concentration
#' (\%, vol/vol) of the designed probes.  Must be greater than the
#' \code{minMeltingFormamide}.
#' @param minScore The minimum score of designed probes before exclusion.  A
#' greater \code{minScore} will accelerate the code because more target sites
#' will be excluded from consideration.  However, if the \code{minScore} is too
#' high it will prevent any target sites from being recorded.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @return A \code{data.frame} with the optimal set of probes matching the
#' specified constraints.  Each row lists the probe's target sequence
#' (\code{name}), \code{start} position, \code{length} in nucleotides, start
#' and end position in the sequence alignment, number of \code{permutations},
#' \code{score}, melt point in percent \code{formamide} at 42 degrees Celsius,
#' hybridization efficiency (\code{hyb_eff}), target site, and probe(s).
#' Probes are designed such that the stringency is determined by the
#' equilibrium hybridization conditions and not subsequent washing steps.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{Array2Matrix}}, \code{\link{NNLS}}
#' @references ES Wright et al. (2013) Identification of Bacterial and Archaeal
#' Communities From Source to Tap. Water Research Foundation, Denver, CO.
#' 
#' DR Noguera, et al. (2014). Mathematical tools to optimize the design of
#' oligonucleotide probes and primers. Applied Microbiology and Biotechnology.
#' doi:10.1007/s00253-014-6165-x.
#' @examples
#' 
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' names(dna) <- 1:length(dna)
#' probes <- DesignArray(dna)
#' probes[1,]
#' 
#' @export DesignArray
DesignArray <- function(myDNAStringSet,
	maxProbeLength=24,
	minProbeLength=20,
	maxPermutations=4,
	numRecordedMismatches=500,
	numProbes=10,
	start=1,
	end=NULL,
	maxOverlap=5,
	hybridizationFormamide=10,
	minMeltingFormamide=15,
	maxMeltingFormamide=20,
	minScore=-1e12,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(hybridizationFormamide))
		stop("hybridizationFormamide must be a positive numeric.")
	if (hybridizationFormamide <= 0)
		stop("hybridizationFormamide must be a positive numeric.")
	if (!is.numeric(minMeltingFormamide))
		stop("minMeltingFormamide must be a positive numeric.")
	if (minMeltingFormamide <= 0)
		stop("minMeltingFormamide must be a positive numeric.")
	if (!is.numeric(maxMeltingFormamide))
		stop("maxMeltingFormamide must be a positive numeric.")
	if (maxMeltingFormamide <= minMeltingFormamide)
		stop("maxMeltingFormamide must be a greater than minMeltingFormamide.")
	if (!is.numeric(minScore))
		stop("minScore must be a numeric.")
	if (minScore >= 100)
		stop("minScore must be less than the maximum possible score (100).")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	numF <- length(myDNAStringSet)
	if (!is.numeric(maxProbeLength))
		stop("maxProbeLength must be a positive integer.")
	if (maxProbeLength < 1)
		stop("maxProbeLength must be a positive integer.")
	if ((maxProbeLength < 19) || (maxProbeLength > 26))
		warning("maxProbeLength is not in ideal range from 19 to 26.")
	if (!is.numeric(minProbeLength))
		stop("minProbeLength must be a positive integer.")
	if (minProbeLength < 1)
		stop("minProbeLength must be a positive integer.")
	if ((minProbeLength < 19) || (minProbeLength > 26))
		warning("minProbeLength is not in ideal range from 19 to 26.")
	if (!is.numeric(maxPermutations))
		stop("maxPermutations must be a positive integer.")
	if (maxPermutations < 1)
		stop("maxPermutations must be a positive integer")
	if (!is.numeric(numRecordedMismatches))
		stop("numRecordedMismatches must be a positive integer.")
	if (numRecordedMismatches < 1)
		stop("numRecordedMismatches must be a positive integer")
	if (!is.numeric(numProbes))
		stop("numProbes must be a positive integer.")
	if (numProbes < 1)
		stop("numProbes must be a positive integer")
	if (!is.numeric(start))
		stop("start must be a positive integer.")
	if (start < 1)
		stop("start must be a positive integer")
	if (!is.numeric(maxOverlap))
		stop("maxOverlap must be a positive integer.")
	if (maxOverlap < 0)
		stop("maxOverlap must be a positive integer")
	if (numRecordedMismatches*length(myDNAStringSet) > (2^31 - 1)/2)
		stop("numRecordedMismatches is too large.")
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
	
	if (is.null(end)) {
		end <- min(width(myDNAStringSet))
	} else {
		if (!is.numeric(end))
			stop("end must be a positive integer.")
		if (end > min(width(myDNAStringSet)))
			stop("end is longer than the shortest sequence.")
	}
	if (!is.numeric(start))
		stop("start must be a positive integer.")
	if (start < 1)
		stop("start must be a positive integer")
	if (end < (start + maxProbeLength))
		stop("end must be greater than start.")
	
	if (numF < 2)
		stop("myDNAStringSet must contain more than one sequence.")
	
	# initialize a progress bar
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
	} else {
		pBar <- NULL
	}
	
	maxW <- unique(width(myDNAStringSet))
	if (length(maxW)!=1)
		stop("\nSequences are not aligned.\n")
	
	probes <- .Call("designProbes",
		myDNAStringSet,
		maxProbeLength,
		minProbeLength,
		maxPermutations,
		numRecordedMismatches,
		numProbes,
		start,
		end,
		maxOverlap,
		hybridizationFormamide,
		minMeltingFormamide,
		maxMeltingFormamide,
		minScore,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	w <- which(probes[[1]][,9] != minScore)
	if (is.null(names(myDNAStringSet)))
		names(myDNAStringSet) <- 1:length(myDNAStringSet)
	f <- function(x, myNames) {
		w1 <- which(x[(numRecordedMismatches + 1):(2*numRecordedMismatches)] != -1)
		w1 <- w1[order(x[w1], decreasing=TRUE)]
		if (length(w1) > 0) {
			return(paste(myNames[x[numRecordedMismatches + w1] + 1],
				" (",
				round(x[w1],
					digits=1),
				"%)",
				sep="",
				collapse=", "))
		} else {
			return("")
		}
	}
	mismatches <- unlist(apply(probes[[4]][w,, drop=FALSE], 1, f, names(myDNAStringSet)))
	
	p <- data.frame(name=I(names(myDNAStringSet)[probes[[1]][w,1] + 1]),
		start=I(probes[[1]][w,2]),
		length=I(probes[[1]][w,3]),
		start_aligned=I(probes[[1]][w,6]),
		end_aligned=I(probes[[1]][w,7]),
		permutations=I(probes[[1]][w,4]),
		score=I(probes[[1]][w,5]),
		formamide=I(probes[[1]][w,8]),
		hyb_eff=I(probes[[1]][w,9]),
		target_site=I(probes[[2]][w]),
		probes=I(probes[[3]][w]),
		mismatches=I(mismatches))
	
	if (verbose) {
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(p)
}
