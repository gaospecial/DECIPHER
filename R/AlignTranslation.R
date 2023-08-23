#' Align Sequences By Their Amino Acid Translation
#' 
#' Performs alignment of a set of DNA or RNA sequences by aligning their
#' corresponding amino acid sequences.
#' 
#' Alignment of proteins is often more accurate than alignment of their coding
#' nucleic acid sequences.  This function aligns the input nucleic acid
#' sequences via aligning their translated amino acid sequences.  First, the
#' input sequences are translated according to the specified \code{sense},
#' \code{direction}, and \code{readingFrame}.  The resulting amino acid
#' sequences are aligned using \code{AlignSeqs}, and this alignment is
#' (conceptually) reverse translated into the original sequence type,
#' \code{sense}, and \code{direction}.  Not only is alignment of protein
#' sequences generally more accurate, but aligning translations will ensure
#' that the reading frame is maintained in the nucleotide sequences.
#' 
#' If the \code{readingFrame} is \code{NA} (the default) then an attempt is
#' made to guess the reading frame of each sequence based on the number of stop
#' codons in the translated amino acids.  For each sequence, the first reading
#' frame will be chosen (either \code{1}, \code{2}, or \code{3}) without stop
#' codons, except in the final position.  If the number of stop codons is
#' inconclusive for a sequence then the reading frame will default to \code{1}.
#' The entire length of each sequence is translated in spite of any stop codons
#' identified.  Note that this method is only constructive in circumstances
#' where there is a substantially long coding sequence with at most a single
#' stop codon expected in the final position, and therefore it is preferable to
#' specify the reading frame of each sequence if it is known.
#' 
#' @name AlignTranslation
#' @param myXStringSet A \code{DNAStringSet} or \code{RNAStringSet} object of
#' unaligned sequences.
#' @param sense Single character specifying sense of the input sequences,
#' either the positive (\code{"+"}) coding strand or negative (\code{"-"})
#' non-coding strand.
#' @param direction Direction of the input sequences, either \code{"5' to 3'"}
#' or \code{"3' to 5'"}.
#' @param readingFrame Numeric vector giving a single reading frame for all of
#' the sequences, or an individual reading frame for each sequence in
#' \code{myXStringSet}.  The \code{readingFrame} can be either \code{1},
#' \code{2}, \code{3} to begin translating on the first, second, and third
#' nucleotide position, or \code{NA} (the default) to guess the reading frame.
#' (See details section below.)
#' @param type Character string indicating the type of output desired.  This
#' should be (an abbreviation of) one of \code{"DNAStringSet"},
#' \code{"RNAStringSet"}, \code{"AAStringSet"}, or \code{"both"}.  (See value
#' section below.)
#' @param geneticCode Either a character vector giving the genetic code in the
#' same format as \code{GENETIC_CODE} (the default), or a list containing one
#' genetic code for each sequence in \code{myXStringSet}.
#' @param \dots Further arguments to be passed directly to
#' \code{\link{AlignSeqs}}, including \code{gapOpening}, \code{gapExtension},
#' \code{gapPower}, \code{terminalGap}, \code{restrict}, \code{anchor},
#' \code{normPower}, \code{substitutionMatrix}, \code{structureMatrix},
#' \code{standardize}, \code{alphabet}, \code{guideTree}, \code{iterations},
#' \code{refinements}, \code{useStructures}, \code{structures}, \code{FUN}, and
#' \code{levels}.
#' @return An \code{XStringSet} of the class specified by \code{type}, or a
#' list of two components (nucleotides and amino acids) if \code{type} is
#' \code{"both"}.  Note that incomplete starting and ending codons will be
#' translated into the mask character ("+") if the result includes an
#' \code{AAStringSet}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignDB}}, \code{\link{AlignProfiles}},
#' \code{\link{AlignSeqs}}, \code{\link{AlignSynteny}},
#' \code{\link{CorrectFrameshifts}}
#' @references Wright, E. S. (2015). DECIPHER: harnessing local sequence
#' context to improve protein multiple sequence alignment. BMC Bioinformatics,
#' 16, 322. http://doi.org/10.1186/s12859-015-0749-z
#' @examples
#' 
#' # first three sequences translate to  MFITP*
#' # and the last sequence translates as MF-TP*
#' rna <- RNAStringSet(c("AUGUUCAUCACCCCCUAA", "AUGUUCAUAACUCCUUGA",
#' 	"AUGUUCAUUACACCGUAG", "AUGUUUACCCCAUAA"))
#' RNA <- AlignSeqs(rna, verbose=FALSE)
#' RNA
#' 
#' RNA <- AlignTranslation(rna, verbose=FALSE)
#' RNA
#' 
#' AA <- AlignTranslation(rna, type="AAStringSet", verbose=FALSE)
#' AA
#' 
#' BOTH <- AlignTranslation(rna, type="both", verbose=FALSE)
#' BOTH
#' 
#' # example of aligning many protein coding sequences:
#' fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' DNA <- AlignTranslation(dna) # align the translation then reverse translate
#' DNA
#' 
#' # using a mixture of standard and non-standard genetic codes
#' gC1 <- getGeneticCode(id_or_name2="1", full.search=FALSE, as.data.frame=FALSE)
#' # Mollicutes use an alternative genetic code
#' gC2 <- getGeneticCode(id_or_name2="4", full.search=FALSE, as.data.frame=FALSE)
#' w <- grep("Mycoplasma|Ureaplasma", names(dna))
#' gC <- vector("list", length(dna))
#' gC[-w] <- list(gC1)
#' gC[w] <- list(gC2)
#' AA <- AlignTranslation(dna, geneticCode=gC, type="AAStringSet")
#' BrowseSeqs(AA)
#' 
#' @export AlignTranslation
AlignTranslation <- function(myXStringSet,
	sense="+",
	direction="5' to 3'",
	readingFrame=NA,
	type=class(myXStringSet),
	geneticCode=GENETIC_CODE,
	...) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (min(width(myXStringSet)) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	if (length(myXStringSet) < 2)
		stop("At least two sequences are required.")
	if (sense != "+" && sense != "-")
		stop('sense must be either "+" or "-".')
	if (direction != "5' to 3'" && direction != "3' to 5'")
		stop('direction must be either "5\' to 3\'" or "3\' to 5\'".')
	if (!is.na(readingFrame) && !is.numeric(readingFrame))
		stop('readingFrame must be a numeric.')
	if (length(readingFrame) != 1 && length(readingFrame) != length(myXStringSet))
		stop('readingFrame must be a single numeric or the length of myXStringSet.')
	if (any(!is.na(readingFrame) & (readingFrame > 3 | readingFrame < 1)))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	if (any(!is.na(readingFrame) & floor(readingFrame) != readingFrame))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before alignment.")
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') must be removed before alignment.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before alignment.")
	if (length(type)==0)
		stop("No type specified.")
	if (length(type) > 1)
		stop("Only one type may be specified.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (type==1) {
		if (!is(myXStringSet, "DNAStringSet"))
			stop("type cannot be 'RNAStringSet' when myXStringSet is a DNAStringSet.")
	} else if (type==2) {
		if (!is(myXStringSet, "RNAStringSet"))
			stop("type cannot be 'DNAStringSet' when myXStringSet is a RNAStringSet.")
	}
	if (is.list(geneticCode)) {
		if (length(geneticCode)!=length(myXStringSet))
			stop("The list geneticCode must have one item per sequence.")
		mapping <- selfmatch(geneticCode)
		mapping <- tapply(seq_along(mapping), mapping, c)
		group <- as.integer(names(mapping))
	} else { # all identifiers use the same geneticCode
		geneticCode <- list(geneticCode)
		mapping <- list(seq_along(myXStringSet))
		group <- 1L
	}
	
	if (sense=="-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction=="3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	if (length(readingFrame)==1)
		readingFrame <- rep(readingFrame, length(myXStringSet))
	
	index <- c(0:2, 0) # circle back around to the first reading frame
	AA <- AAStringSet(rep("", length(myXStringSet)))
	for (i in seq_along(index)) {
		w <- which(((readingFrame - 1)==index[i] & # (specified reading frame AND
			i != length(index)) | # not the last possible reading frame) OR
			is.na(readingFrame)) # reading frame is unspecified
		if (length(w)==0)
			next
		
		start <- index[i] + 1
		for (j in seq_along(group)) {
			W <- w[w %in% mapping[[j]]]
			if (length(W)==0)
				next
			
			end <- width(myXStringSet)[W]
			offset <- end - start + 1
			end <- end - offset %% 3
			end <- ifelse(end < start - 1,
				start - 1,
				end)
			AA[W] <- translate(subseq(myXStringSet[W],
					start,
					end),
				genetic.code=geneticCode[[group[j]]],
				if.fuzzy.codon="solve",
				no.init.codon=TRUE)
		}
		end <- width(myXStringSet)[w]
		offset <- end - start + 1
		end <- end - offset %% 3
		end <- ifelse(end < start - 1,
			start - 1,
			end)
		
		# mask missing positions
		AA[w] <- xscat(ifelse(start > 1, "+", ""),
			AA[w],
			ifelse(end < width(myXStringSet)[w], "+", ""))
		
		a <- vcountPattern("*", AA[w])
		lastResidue <- substring(AA[w],
			width(AA)[w],
			width(AA)[w])
		readingFrame[w] <- ifelse(a==0 | # no stops OR
			(a==1 & lastResidue=="*") | # one stop at end OR
			i==length(index), # already checked all reading frames
			index[i] + 1,
			readingFrame[w])
	}
	
	names(AA) <- names(myXStringSet)
	AA <- AlignSeqs(myXStringSet=AA, ...)
	if (type > 2) {
		if (type==4) {
			results <- list(NULL, AA)
		} else {
			return(AA)
		}
	}
	
	# correct for shifts introduced by "+"
	readingFrame <- ifelse(readingFrame==1,
		1,
		ifelse(readingFrame==2,
			-1,
			0))
	
	gaps <- vmatchPattern("-", AA)
	starts <- list()
	maxRF <- max(readingFrame)
	Ls <- numeric(length(myXStringSet))
	for (i in seq_along(myXStringSet)) {
		start <- start(gaps[[i]])
		if (length(start) > 0) {
			w <- which(start + (length(start) - 1):0==width(AA)[i])
			if (length(w) > 0)
				length(start) <- length(start) - length(w)
			start <- start*3 - 3 + readingFrame[i]
			start <- sort(c(start, start + 1, start + 2))
			start <- start - 0:(length(start) - 1)
			w <- which(start==(readingFrame[i]))
			if (length(w) > 0)
				start[w] <- 1
		}
		start <- c(rep(1, maxRF - readingFrame[i]), start)
		starts[[i]] <- start
		Ls[i] <- length(start) + width(myXStringSet)[i]
	}
	
	# add trailing gaps
	maxWidth <- max(Ls)
	for (i in seq_along(myXStringSet)) {
		starts[[i]] <- c(starts[[i]],
			rep(width(myXStringSet)[i] + 1, maxWidth - Ls[i]))
	}
	
	# remove leading common gaps
	start <- min(unlist(lapply(starts, function(x) return(length(which(x==1))))))
	if (start > 0)
		starts <- lapply(starts, function(x) return(x[-(1:start)]))
	
	myXStringSet <- replaceAt(myXStringSet,
		starts,
		"-")
	if (sense=="-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction=="3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	
	if (type==4) {
		results[[1]] <- myXStringSet
		return(results)
	} else {
		return(myXStringSet)
	}
}
