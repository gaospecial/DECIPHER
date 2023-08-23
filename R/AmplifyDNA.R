#' Simulate Amplification of DNA by PCR
#' 
#' Predicts the amplification efficiency of theoretical PCR products
#' (amplicons) given one or more primer sequences.
#' 
#' Exponential amplification in PCR requires the annealing and elongation of
#' two primers from target sites on opposing strands of the template DNA.  If
#' the template DNA sequence (e.g., chromosome) is known then predictions of
#' theoretical amplicons can be obtained from in silico simulations of
#' amplification.  \code{AmplifyDNA} first searches for primer target sites on
#' the template DNA, and then calculates an amplification efficiency from each
#' target site using \code{\link{CalculateEfficiencyPCR}}.  Ambiguity codes
#' (\code{IUPAC_CODE_MAP}) are supported in the \code{primers}, but not in
#' \code{myDNAStringSet} to prevent trivial matches (e.g., runs of N's).
#' 
#' If \code{taqEfficiency} is \code{TRUE} (the default), the amplification
#' efficiency of each primer is defined as the product of hybridization
#' efficiency and elongation efficiency.  Amplification efficiency must be at
#' least \code{minEfficiency} for a primer to be amplified in silico.  Overall
#' amplification efficiency of the PCR product is then calculated as the
#' geometric mean of the two (i.e., forward and reverse) primers' efficiencies.
#' Finally, amplicons are generated if the two primers are within
#' \code{maxProductSize} nucleotides downstream of each other.
#' 
#' Potential PCR products are returned, either with or without including the
#' primer sequences in the amplicon.  The default (\code{includePrimers=TRUE})
#' is to incorporate the primer sequences as would normally occur during
#' amplification.  The alternative is to return the complete template sequence
#' including the target sites, which may not exactly match the primer
#' sequences.  Note that amplicons may be duplicated when different input
#' \code{primers} can amplify the same region of DNA.
#' 
#' @name AmplifyDNA
#' @param primers A \code{DNAStringSet} object or character vector with one or
#' more unaligned primer sequences in 5' to 3' orientation.
#' @param myDNAStringSet A \code{DNAStringSet} object or character vector with
#' unaligned template DNA sequences in 5' to 3' orientation.
#' @param maxProductSize Integer specifying the maximum length of PCR products
#' (amplicons) in nucleotides.
#' @param annealingTemp Numeric specifying the annealing temperature used in
#' the PCR reaction.
#' @param P Numeric giving the molar concentration of primers in the reaction.
#' @param ions Numeric giving the molar sodium equivalent ionic concentration.
#' Values may range between 0.01M and 1M.
#' @param includePrimers Logical indicating whether to include the primer
#' sequences in the theoretical PCR products.  (See details section below.)
#' @param minEfficiency Numeric giving the minimum amplification efficiency of
#' PCR products to include in the output (default \code{0.1%}).  (See details
#' section below.)
#' @param \dots Additional arguments to be passed directly to
#' \code{\link{CalculateEfficiencyPCR}}, including \code{batchSize},
#' \code{taqEfficiency}, \code{maxDistance}, \code{maxGaps}, and
#' \code{processors}.
#' @return A \code{DNAStringSet} object containing potential PCR products,
#' sorted from highest-to-lowest amplification efficiency.  The sequences are
#' named by their predicted amplification efficiency followed by the index of
#' each primer in \code{primers} and the name (or index if names are missing)
#' of the amplified sequence in \code{myDNAStringSet}.  (See examples section
#' below.)
#' @note The program OligoArrayAux
#' (\url{http://www.unafold.org/Dinamelt/software/oligoarrayaux.php}) must be
#' installed in a location accessible by the system.  For example, the
#' following code should print the installed OligoArrayAux version when
#' executed from the R console:
#' 
#' \code{system("hybrid-min -V")}
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{CalculateEfficiencyPCR}}, \code{\link{DesignPrimers}},
#' \code{\link{DesignSignatures}}, \code{\link{MeltDNA}}
#' @references ES Wright et al. (2013) "Exploiting Extension Bias in PCR to
#' Improve Primer Specificity in Ensembles of Nearly Identical DNA Templates."
#' Environmental Microbiology, doi:10.1111/1462-2920.12259.
#' @examples
#' 
#' data(yeastSEQCHR1)
#' 
#' # not run (must have OligoArrayAux installed first):
#' 
#' # match a single primer that acts as both the forward and reverse
#' primer1 <- "TGGAAGCTGAAACG"
#' \dontrun{AmplifyDNA(primer1, yeastSEQCHR1, annealingTemp=55, P=4e-7, maxProductSize=500)}
#' 
#' # perform a typical amplification with two primer sequences:
#' primer2 <- c("GGCTGTTGTTGGTGTT", "TGTCATCAGAACACCAA")
#' \dontrun{AmplifyDNA(primer2, yeastSEQCHR1, annealingTemp=55, P=4e-7, maxProductSize=500)}
#' 
#' # perform a multiplex PCR amplification with multiple primers:
#' primers <- c(primer1, primer2)
#' \dontrun{AmplifyDNA(primers, yeastSEQCHR1, annealingTemp=55, P=4e-7, maxProductSize=500)}
#' 
#' @export AmplifyDNA
AmplifyDNA <- function(primers,
	myDNAStringSet,
	maxProductSize,
	annealingTemp,
	P,
	ions=0.2,
	includePrimers=TRUE,
	minEfficiency=0.001,
	...) {
	
	# error checking
	if (is.character(primers))
		primers <- DNAStringSet(primers)
	if (!is(primers, "DNAStringSet"))
		stop("primers must be a DNAStringSet.")
	if (any(primers==""))
		stop("primers cannot contain 0 character elements.")
	if (is.character(myDNAStringSet))
		myDNAStringSet <- DNAStringSet(myDNAStringSet)
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (length(myDNAStringSet)==0)
		stop("myDNAStringSet is empty.")
	if (!is.numeric(ions))
		stop("ions must be a numeric.")
	if (ions < .01 || is.nan(ions))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	if (!is.numeric(P))
		stop("P must be a numeric.")
	if (!(P > 0))
		stop("P must be greater than zero.")
	if (!is.numeric(annealingTemp))
		stop("annealingTemp must be a numeric.")
	if (!is.numeric(maxProductSize))
		stop("maxProductSize must be a numeric.")
	if (maxProductSize <= 0)
		stop("maxProductSize must be greater than zero")
	maxProductSize <- as.integer(maxProductSize)
	if (!is.logical(includePrimers))
		stop("includePrimers must be a logical.")
	a <- vcountPattern("-", myDNAStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before amplification.")
	a <- vcountPattern("+", myDNAStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') must be removed before amplification.")
	a <- vcountPattern(".", myDNAStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before amplification.")
	
	primers <- Disambiguate(primers)
	l <- unlist(lapply(primers, length))
	l <- rep(seq_along(primers), l)
	primers <- unlist(primers)
	
	ns <- names(myDNAStringSet)
	w <- width(myDNAStringSet)
	w <- cumsum(w)
	w <- w[-length(w)]
	starts <- c(1L,
		w + seq_along(w)*maxProductSize + 1L)
	myDNAStringSet <- unlist(myDNAStringSet)
	if (length(w) > 0)
		myDNAStringSet <- replaceAt(myDNAStringSet,
			IRanges(start=w + 1L,
				width=0),
			paste(rep("-", maxProductSize),
				collapse=""))
	
	f <- IRanges()
	fp <- integer()
	for (i in 1:length(primers)) {
		temp <- matchPattern(primers[[i]],
			myDNAStringSet,
			max.mismatch=ceiling(0.25*nchar(primers[[i]])),
			with.indels=TRUE)
		temp <- as(temp, "IRanges")
		f <- c(f, temp)
		fp <- c(fp, rep(i, length(temp)))
	}
	f <- Views(myDNAStringSet, start(f), end(f))
	if (length(f)==0)
		return(DNAStringSet())
	
	r <- IRanges()
	rp <- integer()
	for (i in 1:length(primers)) {
		temp <- matchPattern(reverseComplement(primers[[i]]),
			myDNAStringSet,
			max.mismatch=ceiling(.25*nchar(primers[[i]])),
			with.indels=TRUE)
		temp <- as(temp, "IRanges")
		r <- c(r, temp)
		rp <- c(rp, rep(i, length(temp)))
	}
	r <- Views(myDNAStringSet, start(r), end(r))
	if (length(r)==0)
		return(DNAStringSet())
	
	ends <- end(r)
	o <- order(ends)
	r <- r[o]
	rp <- rp[o]
	
	targets <- extractAt(myDNAStringSet,
		at=IRanges(start=start(f),
			end=end(f)))
	fe <- CalculateEfficiencyPCR(primers[fp],
		reverseComplement(targets),
		annealingTemp, P, ions, ...)
	fw <- which(fe >= minEfficiency)
	if (length(fw)==0)
		return(DNAStringSet())
	
	targets <- extractAt(myDNAStringSet,
		at=IRanges(start=start(r),
			end=end(r)))
	re <- CalculateEfficiencyPCR(primers[rp],
		targets,
		annealingTemp, P, ions, ...)
	rw <- which(re >= minEfficiency)
	if (length(rw)==0)
		return(DNAStringSet())
	
	sf <- start(f)[fw]
	sr <- start(r)[rw]
	ef <- end(f)[fw]
	er <- end(r)[rw]
	b <- e <- fs <- rs <- integer(1e6)
	effs <- numeric(1e6)
	c <- 0L
	for (i in 1:length(sf)) {
		w <- .Call("boundedMatches",
			sr,
			as.integer(sf[i] + 1L),
			as.integer(sf[i] + maxProductSize - 1L),
			PACKAGE="DECIPHER")
		if (length(w) > 0) {
			r <- (c + 1L):(c + length(w))
			c <- c + length(w)
			if (includePrimers) {
				b[r] <- ef[i] + 1L
				e[r] <- sr[w] - 1L
			} else {
				b[r] <- sf[i]
				e[r] <- er[w]
			}
			effs[r] <- sqrt(fe[fw[i]]*re[rw[w]])
			fs[r] <- fp[fw[i]]
			rs[r] <- rp[rw[w]]
		}
	}
	
	if (c > 0) {
		length(b) <- length(e) <- length(effs) <- c
		
		o <- order(effs, decreasing=TRUE)
		
		index <- sapply(b[o],
			function(x) {
				tail(which(starts <= x), n=1)
			})
		if (!is.null(ns))
			index <- ns[index]
		
		if (includePrimers) {
			extra <- ifelse(b > e, b - e - 1L, 0) # primers overlap
			e <- ifelse(b > e, b - 1L, e) # primers overlap
			v <- Views(myDNAStringSet,
				start=b[o],
				end=e[o])
			v <- as(v, "DNAStringSet")
			v <- xscat(substring(primers[fs[o]],
					1L,
					width(primers[fs[o]]) - extra[o]),
				v,
				reverseComplement(primers[rs[o]]))
			names(v) <- paste(round(100*effs[o], 1),
				"% (",
				l[fs[o]],
				" x ",
				l[rs[o]],
				") ",
				index,
				sep="")
		} else {
			v <- Views(myDNAStringSet,
				start=b[o],
				end=e[o],
				names=paste(round(100*effs[o], 1),
					"% (",
					l[fs[o]],
					" x ",
					l[rs[o]],
					") ",
					index,
					sep=""))
			v <- as(v, "DNAStringSet")
		}
		
		w <- which(width(v) > maxProductSize)
		if (length(w) > 0)
			v <- v[-w]
	} else {
		v <- DNAStringSet()
	}
	
	return(v)
}
