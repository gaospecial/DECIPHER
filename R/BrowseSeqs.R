#' View Sequences in a Web Browser
#' 
#' Opens an html file in a web browser to show the sequences in an
#' \code{XStringSet}.
#' 
#' \code{BrowseSeqs} converts an \code{XStringSet} into html format for viewing
#' in a web browser.  The sequences are colored in accordance with the
#' \code{patterns} that are provided, or left uncolored if \code{colorPatterns}
#' is \code{FALSE} or \code{patterns} is \code{NULL}.  Character or
#' \code{XStringSet} patterns are matched as regular expressions and colored
#' according to \code{colors}.  If \code{patterns} is a list of matrices, then
#' it must contain one element per sequence.  Each matrix is interpreted as
#' providing the fraction red, blue, and green for each letter in the sequence.
#' Thus, \code{colors} is ignored when \code{patterns} is a list.  (See
#' examples section below.)
#' 
#' Patterns are not matched across column breaks, so multi-character
#' \code{patterns} should be carefully considered when \code{colWidth} is less
#' than the maximum sequence length.  Patterns are matched sequentially in the
#' order provided, so it is feasible to use nested \code{patterns} such as
#' \code{c("ACCTG", "CC")}.  In this case the ``CC'' could be colored
#' differently inside the previously colored ``ACCTG''.  Note that
#' \code{patterns} overlapping the boundaries of a previously matched pattern
#' will not be matched.  For example, ``ACCTG'' would not be matched if
#' \code{patterns=c("CC", "ACCTG")}.
#' 
#' Some web browsers cannot quickly display a large amount colored text, so it
#' is recommended to use \code{colorPatterns = FALSE} or to \code{highlight} a
#' sequence when viewing a large \code{XStringSet}.  Highlighting will only
#' show all of the characters in the highlighted sequence, and convert all
#' matching positions in the other sequences into dots without \code{color}.
#' Also, note that some web browsers display small shifts between fixed-width
#' characters that may become noticeable as color offsets between the ends of
#' long sequences.
#' 
#' @name BrowseSeqs
#' @param myXStringSet A \code{XStringSet} object of sequences.
#' @param htmlFile Character string giving the location where the html file
#' should be written.
#' @param openURL Logical indicating whether the \code{htmlFile} should be
#' opened in a web browser.
#' @param colorPatterns Logical specifying whether to color matched
#' \code{patterns}, or an integer vector providing pairs of start and stop
#' boundaries for coloring.
#' @param highlight Numeric specifying which sequence in the set to use for
#' comparison or \code{NA} to color all sequences (default).  If
#' \code{highlight} is \code{0} then positions differing from the consensus
#' sequence are highlighted.
#' @param patterns Either an \code{AAStringSet}, \code{DNAStringSet}, or
#' \code{RNAStringSet} object, a character vector containing regular
#' expressions, a list of numeric matrices, or \code{NULL}.  (See details
#' section below.)
#' @param colors Character vector providing the color for each of the matched
#' \code{patterns}.  Typically a character vector with elements of 7
#' characters: ``#'' followed by the red, blue, green values in hexadecimal
#' (after rescaling to 0 ... 255).  Ignored when \code{patterns} is a list of
#' matrices.
#' @param colWidth Integer giving the maximum number of nucleotides wide the
#' display can be before starting a new page.  Must be a multiple of \code{20}
#' (e.g., \code{100}), or \code{Inf} (the default) to display all the sequences
#' in one set of rows.
#' @param \dots Additional arguments to adjust the appearance of the consensus
#' sequence at the base of the display.  Passed directly to
#' \code{ConsensusSequence} for an \code{AAStringSet}, \code{DNAStringSet}, or
#' \code{RNAStringSet}, or to \code{consensusString} for a \code{BStringSet}.
#' @return Creates an html file containing sequence data and (if \code{openURL}
#' is \code{TRUE}) opens it in a web browser for viewing.  The layout has the
#' sequence name on the left, position legend on the top, cumulative number of
#' nucleotides on the right, and consensus sequence on the bottom.
#' 
#' Returns \code{htmlFile} if the html file was written successfully.
#' @note Some web browsers do not display colored characters with equal widths.
#' If positions do not align across sequences then try opening the
#' \code{htmlFile} with a different web browser.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{BrowseDB}}, \code{\link{ConsensusSequence}}
#' @references ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological
#' Sequence Data in R". The R Journal, \bold{8(1)}, 352-359. Kunzmann P., et
#' al. (2020) "Substitution matrix based color schemes for sequence alignment
#' visualization". BMC Bioinformatics, \bold{21(1):209}.
#' @examples
#' 
#' # load the example DNA sequences
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' dna <- SearchDB(db) # non-coding ribosomal RNA gene sequences
#' 
#' # example of using the defaults with DNA sequences
#' BrowseSeqs(dna) # view the XStringSet
#' 
#' # color only "ACTG" and "CSC" patterns (where S is C or G)
#' BrowseSeqs(dna, patterns=DNAStringSet(c("ACTG", "CSC")))
#' 
#' # highlight (i.e., only fully-color) the first sequence
#' BrowseSeqs(dna, highlight=1) # other sequences are dots where matching
#' 
#' # highlight the consensus sequence at the bottom
#' BrowseSeqs(dna, highlight=0) # other sequences are dots where matching
#' 
#' # split the wide view into multiple vertical pages (for printing)
#' BrowseSeqs(dna, colWidth=100, highlight=1)
#' 
#' # specify an alternative color scheme for -, A, C, G, T
#' BrowseSeqs(dna, colors=c("#1E90FF", "#32CD32", "#9400D3", "black", "#EE3300"))
#' 
#' # only color the positions within certain positional ranges (100-200 & 250-500)
#' BrowseSeqs(dna, colorPatterns=c(100, 200, 250, 500))
#' 
#' # example of calling attention to letters by coloring gaps black
#' BrowseSeqs(dna, patterns="-", colors="black")
#' 
#' \dontrun{
#' # color according to base-pairing by supplying the fraction RGB for every position
#' dbn <- PredictDBN(dna, type="structures") # calculate the secondary structures
#' # dbn now contains the scores for whether a base is paired (left/right) or unpaired
#' dbn[[1]][, 1] # the scores for the first position in the first sequence
#' dbn[[2]][, 10] # the scores for the tenth position in the second sequence
#' # these positional scores can be used as shades of red, green, and blue:
#' BrowseSeqs(dna, patterns=dbn) # red = unpaired, green = left-pairing, blue = right
#' # positions in black are not part of the consensus secondary structure
#' }
#' 
#' # color all restriction sites
#' data(RESTRICTION_ENZYMES) # load dataset containing restriction enzyme sequences
#' sites <- RESTRICTION_ENZYMES
#' sites <- gsub("[^A-Z]", "", sites) # remove non-letters
#' sites <- DNAStringSet(sites) # convert the character vector to a DNAStringSet
#' rc_sites <- reverseComplement(DNAStringSet(sites))
#' w <- which(sites != rc_sites) # find non-palindromic restriction sites
#' sites <- c(sites, rc_sites[w]) # append their reverse complement
#' sites <- sites[order(nchar(sites))] # match shorter sites first
#' BrowseSeqs(dna, patterns=sites)
#' 
#' # color bases by quality score
#' fastq <- system.file("extdata", "s_1_sequence.txt", package="Biostrings")
#' reads <- readQualityScaledDNAStringSet(fastq, quality.scoring="solexa")
#' colors <- colorRampPalette(c("red", "yellow", "green"))(42)
#' colors <- col2rgb(colors)/255
#' quals <- as(quality(reads), "IntegerList")
#' quals <- lapply(quals, function(x) colors[, x])
#' BrowseSeqs(DNAStringSet(reads), patterns=quals) # green = high quality, red = low quality
#' 
#' # load the example protein coding sequences
#' fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # example of using the defaults with amino acid sequences
#' aa <- unique(translate(dna)) # the unique amino acid sequences
#' BrowseSeqs(aa)
#' 
#' # example of highlighting the consensus amino acid sequence
#' AA <- AlignSeqs(aa)
#' BrowseSeqs(AA, highlight=0)
#' 
#' # example of highlighting positions that differ from the majority consensus
#' BrowseSeqs(AA, highlight=0, threshold=0.5)
#' 
#' # specify an alternative color scheme for amino acids (from Kunzmann et al.)
#' colors <- c(`-`="#000000", `A`="#BDB1E8", `R`="#EFA2C5", `N`="#F6602F",
#'     `D`="#FD5559", `C`="#12C7FE", `Q`="#DDACB4", `E`="#FEA097", `G`="#F46802",
#'     `H`="#FCA708", `I`="#369BD9", `L`="#2E95EC", `K`="#CF7690", `M`="#4B8EFE",
#'     `F`="#76997D", `P`="#FD2AE3", `S`="#A08A9A", `T`="#9A84D5", `W`="#74C80D",
#'     `Y`="#9BB896", `V`="#89B9F9")
#' BrowseSeqs(AA, colors=colors, patterns=names(colors))
#' 
#' # example of coloring in a reduced amino acid alphabet
#' alpha <- AA_REDUCED[[15]]
#' alpha # clustering of amino acids based on similarity
#' BrowseSeqs(AA, patterns=c("-", paste("[", alpha, "]", sep="")))
#' 
#' # color amino acids according to their predicted secondary structure
#' hec <- PredictHEC(AA, type="probabilities") # calculate the secondary structures
#' # hec now contains the probability that a base is in an alpha-helix or beta-sheet
#' hec[[3]][, 18] # the 18th position in sequence 3 is likely part of a beta-sheet (E)
#' # the positional probabilities can be used as shades of red, green, and blue:
#' BrowseSeqs(AA, patterns=hec) # red = alpha-helix, green = beta-sheet, blue = coil
#' 
#' # color codons according to their corresponding amino acid
#' DNA <- AlignTranslation(dna) # align the translation then reverse translate
#' colors <- rainbow(21, v=0.8, start=0.9, end=0.7) # codon colors
#' m <- match(GENETIC_CODE, unique(GENETIC_CODE)) # corresponding amino acid
#' codonBounds <- matrix(c(seq(1, width(DNA)[1], 3), # start of codons
#' 	seq(3, width(DNA)[1], 3)), # end of codons
#' 	nrow=2,
#' 	byrow=TRUE)
#' BrowseSeqs(DNA,
#' 	colorPatterns=codonBounds,
#' 	patterns=c("---", names(GENETIC_CODE)), # codons to color
#' 	colors=c("black", substring(colors[m], 1, 7)))
#' 
#' @export BrowseSeqs
BrowseSeqs <- function(myXStringSet,
	htmlFile=paste(tempdir(),"/myXStringSet.html",sep=""),
	openURL=interactive(),
	colorPatterns=TRUE,
	highlight=NA,
	patterns=c("-", alphabet(myXStringSet, baseOnly=TRUE)),
	colors=substring(rainbow(length(patterns), v=0.8, start=0.9, end=0.7), 1, 7),
	colWidth=Inf,
	...) {
	
	# error checking
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
	if (length(myXStringSet)==0)
		stop("No sequence information to display.")
	if (is(patterns, "DNAStringSet")) {
		if (!is(myXStringSet, "DNAStringSet"))
			stop("patterns must be a DNAStringSet.")
		type <- 1L
	} else if (is(patterns, "RNAStringSet")) {
		if (!is(myXStringSet, "RNAStringSet"))
			stop("patterns must be a RNAStringSet.")
		type <- 2L
	} else if (is(patterns, "AAStringSet")) {
		if (!is(myXStringSet, "AAStringSet"))
			stop("patterns must be a AAStringSet.")
		type <- 3L
	} else if (is(patterns, "list")) {
		type <- -1L
	} else {
		type <- 0L
	}
	if (type > 0L) {
		patterns <- as.character(patterns)
		if (type==1L) { # DNAStringSet
			patterns <- gsub("M",
				"[ACM]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("R",
				"[AGR]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("W",
				"[ATW]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("S",
				"[CGS]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Y",
				"[CTY]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("K",
				"[GTK]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("V",
				"[ACGMRSV]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("H",
				"[ACTMWYH]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("D",
				"[AGTRWKD]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("B",
				"[CGTSYKB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("N",
				"[ACGTMRWSYKVHDBN]",
				patterns,
				fixed=TRUE)
		} else if (type==2L) { # RNAStringSet
			patterns <- gsub("M",
				"[ACM]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("R",
				"[AGR]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("W",
				"[AUW]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("S",
				"[CGS]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Y",
				"[CUY]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("K",
				"[GUK]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("V",
				"[ACGMRSV]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("H",
				"[ACUMWYH]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("D",
				"[AGURWKD]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("B",
				"[CGUSYKB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("N",
				"[ACGUMRWSYKVHDBN]",
				patterns,
				fixed=TRUE)
		} else { # AAStringSet
			patterns <- gsub("B",
				"[NDB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Z",
				"[QEZ]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("J",
				"[ILJ]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("X",
				"[ARNDCQEGHILKMFPSTWYVUOBJZX]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("*",
				"\\*",
				patterns,
				fixed=TRUE)
		}
	} else if (type==0) { # character vector
		if (!is.null(patterns) && !is.character(patterns))
			stop("patterns must be a character vector.")
		if (any(grepl("=|\"|<|>|[1-9]|[a-z]", patterns)))
			stop("patterns cannot contain numbers, lower case characters, or the characters (=, <, >, \").")
	}
	if (any(patterns==""))
		stop("patterns cannot be empty (i.e., '').")
	if (type < 0) {
		if (length(myXStringSet) != length(patterns))
			stop("patterns is not the same length as myXStringSet.")
	} else {
		w <- which(patterns %in% c("?", "*", "+", "."))
		if (length(w) > 0)
			patterns[w] <- paste("\\", patterns[w], sep="")
		if (length(colors) != length(patterns) || !is.character(colors))
			stop("colors must be a character vector of the same length as patterns.")
	}
	# check that the file exist
	if (is.character(htmlFile)) {
		htmlfile <- file(htmlFile, "w")
		on.exit(close(htmlfile))
	} else if (!inherits(htmlFile, "connection")) {
		stop("htmlFile must be a character string or connection.")
	}
	if (!is.logical(openURL) || is.na(openURL))
		stop("openURL must be TRUE or FALSE.")
	if (!is.logical(colorPatterns) && !is.numeric(colorPatterns))
		stop("colorPatterns must be a logical or numeric.")
	if (is.numeric(colorPatterns)) {
		if ((length(colorPatterns) %% 2) == 1 || length(colorPatterns) == 0)
			stop("colorPatterns must specify all start and endpoints.")
		if (any((colorPatterns[seq(2, length(colorPatterns), 2)] - colorPatterns[seq(1, length(colorPatterns), 2)]) < 0))
			stop("colorPatterns specifies a negative range.")
		if (any(colorPatterns <= 0))
			stop("colorPatterns must be a positive numeric.")
		if (length(colorPatterns) > 2)
			if (any((colorPatterns[seq(3, length(colorPatterns), 2)] - colorPatterns[seq(2, length(colorPatterns) - 2, 2)]) <= 0))
				stop("Ranges specified in colorPatterns must be non-overlapping.")
		if (max(colorPatterns) > max(width(myXStringSet)))
			stop("Ranges specified in colorPatterns are out-of-bounds.")
	}
	if (is.numeric(colorPatterns) & !is.infinite(colWidth))
		stop("colWidth must be Inf if colorPatterns is numeric.")
	if (is.null(names(myXStringSet))) {
		names(myXStringSet) <- 1:length(myXStringSet)
	} else {
		names(myXStringSet) <- iconv(names(myXStringSet),
			from="",
			to="UTF-8",
			sub=" ")
		names(myXStringSet) <- gsub("\t", " ", names(myXStringSet), fixed=TRUE)
		names(myXStringSet) <- gsub("\n", " ", names(myXStringSet), fixed=TRUE)
	}
	if (!is.na(highlight)) {
		if (highlight < 0 || highlight > length(myXStringSet))
			stop("highlight must be 0 or the index of a sequence in myXStringSet.")
		if (highlight != floor(highlight))
			stop("highlight be be a whole number.")
	}
	# add a consensus sequence to myXStringSet
	if (is(myXStringSet, "DNAStringSet") || is(myXStringSet, "RNAStringSet") || is(myXStringSet, "AAStringSet")) {
		myXStringSet <- c(myXStringSet,
			ConsensusSequence(myXStringSet,
				...))
	} else {
		myXStringSet <- c(myXStringSet,
			as(consensusString(myXStringSet,
					...),
				class(myXStringSet)))
	}
	l <- length(myXStringSet)
	names(myXStringSet)[l] <- "Consensus"
	
	# convert the XStringSet to character strings
	html <- as.character(myXStringSet)
	
	if (type < 0) { # record positions of letter colors
		s <- strsplit(html, "", fixed=TRUE)
		n <- nchar(html)
		v <- vector("list", length(html) - 1) # ignore the consensus
		for (i in seq_len(length(html) - 1)) {
			if (!is.matrix(patterns[[i]]))
				stop("All elements of the list patterns contain a single matrix.")
			if (nrow(patterns[[i]]) != 3)
				stop("All elements of the list patterns must be a matrix with 3 rows.")
			if (any(patterns[[i]] < 0) || any(patterns[[i]] > 1))
				stop("All elements of patterns[[", i, "]] are not between 0 and 1.")
			v[[i]] <- character(n[i])
			if (any(patterns[[i]] > 1) || any(patterns[[i]] < 0))
				stop("All values of patterns[[", i, "]] must be between 0 and 1.")
			w <- which(s[[i]] %in% c(LETTERS, letters, "*"))
			if (length(w) != ncol(patterns[[i]]))
				stop("The number of columns in patterns[[", i, "]] is different than the number of letters in myXStringSet[", i, "]")
			if (length(w) > 0)
				v[[i]][w] <- apply(patterns[[i]],
					2,
					function(x) {
						rgb(x[1], x[2], x[3])
					})
			if (is.numeric(colorPatterns)) {
				start <- 0L
				cPs <- c(colorPatterns, length(v[[i]]) + 1)
				for (j in seq_len(length(cPs))) {
					if (j %% 2) {
						end <- cPs[j]
						if (end > length(v[[i]]))
							end <- length(v[[i]]) + 1L
						if (start < end - 1)
							v[[i]][(start + 1L):(end - 1L)] <- ""
					} else {
						start <- cPs[j]
						if (start > length(v[[i]]) + 1)
							break
					}
				}
			}
		}
	}
	
	if (!is.na(highlight)) { # highlight a sequence
		if (highlight==0) {
			highlight <- length(html)
			index <- 1:(length(html) - 1L)
		} else {
			index <- (1:(length(html) - 1L))[-highlight]
		}
		html <- sapply(html, strsplit, split="", fixed=TRUE)
		for (i in index) {
			L <- min(length(html[[highlight]]), length(html[[i]]))
			w <- which(html[[i]][1:L]==html[[highlight]][1:L])
			if (length(w) > 0)
				html[[i]][w] <- "."
		}
		html <- sapply(html, paste, collapse="")
	}
	
	# pad shorter sequences with spaces
	maxW <- max(width(myXStringSet))
	if (maxW==0)
		stop("No sequence information to display.")
	if (colWidth > maxW)
		colWidth <- maxW
	for (i in seq_len(l)) {
		html[i] <- paste(html[i],
			paste(rep(" ", maxW - nchar(html[i])), collapse=""),
			sep="")
	}
	
	# create a legend that gives position
	if (maxW < 20) {
		if (maxW < 10) {
			counter <- maxW
		} else {
			counter <- 10
		}
	} else {
		counter <- 20
	}
	offset <- (counter - 1) - nchar(maxW)
	if (offset < 0)
		offset <- 0
	legend <- paste(paste(rep(" ", offset), collapse=""), format(seq(counter, maxW, by=counter)), collapse="")
	counter <- ceiling(counter/2)
	tickmarks <- paste(paste(rep("'", counter - 1), collapse=""), rep("|", floor(maxW/counter)), collapse="", sep="")
	tickmarks <- paste(tickmarks, paste(rep("'", maxW - counter*floor(maxW/counter)), collapse=""), sep="")
	
	# split the html into multiple pages
	starts <- seq(1, maxW, colWidth)
	stops <- starts + colWidth - 1L
	stops[length(stops)] <- maxW
	temp <- character((length(html) + 5L)*length(starts))
	count <- 1L
	for (i in 1:length(starts)) {
		temp[count:(count + 1L)] <- substring(c(legend, tickmarks), starts[i], stops[i])
		count <- count + 2L
		temp[c(count:(count + length(html) - 2L), count + length(html))] <- substring(html, starts[i], stops[i])
		count <- count + length(html) + 3L
	}
	html <- temp
	
	# add the cumulative nucleotide lengths on the right
	if (length(starts) > 1) {
		lengths <- numeric(length(starts)*l)
		myLengths <- character(length(starts)*(l + 5L))
		for (i in 1:length(starts)) {
			s <- substring(myXStringSet,
				starts[i],
				stops[i])
			lengths[((i - 1L)*l + 1L):(i*l)] <- nchar(s) - rowSums(letterFrequency(BStringSet(s),
				c("-", " ", ".", "+")))
			if (i > 1)
				lengths[((i - 1L)*l + 1L):(i*l)] <- lengths[((i - 1L)*l + 1L):(i*l)] + lengths[((i - 2L)*l + 1L):((i - 1L)*l)]
			myLengths[((i - 1L)*(l + 5L) + 3L):(i*(l + 5L) - 4L)] <- lengths[((i - 1L)*l + 1L):(i*l - 1L)]
			myLengths[(i*(l + 5L) - 2L)] <- lengths[i*l]
		}
	} else {
		lengths <- width(myXStringSet) - rowSums(letterFrequency(myXStringSet, c("-", ".", "+")))
		myLengths <- c("", "", lengths[-seq(l, length(lengths), l)], "", lengths[seq(l, length(lengths), l)], "", "")
	}
	
	if (is.numeric(colorPatterns) || colorPatterns) {
		if (type < 0) {
			s <- strsplit(html, "", fixed=TRUE)
			for (i in seq_along(v)) {
				count <- 2L
				for (j in seq_along(starts)) {
					w <- which(v[[i]][starts[j]:stops[j]] != "")
					if (length(w) > 0)
						s[[i + count]][w] <- paste("<span style=\"color:#FFF;background:",
							v[[i]][w + starts[j] - 1],
							"\">",
							s[[i + count]][w],
							"</span>",
							sep="")
					count <- count + length(myXStringSet) + 5L
				}
			}
			html <- sapply(s, paste, collapse="")
		} else {
			patterns <- paste("(",
				patterns,
				ifelse(nchar(patterns)==1, "+)", ")"),
				sep="")
			classes <- paste("<span class=\"_",
				seq_along(patterns),
				"\">\\1</span>",
				sep="")
			if (is.numeric(colorPatterns)) {
				htm <- substring(html, 0, colorPatterns[1] - 1)
				for (i in seq(1, length(colorPatterns), 2)) {
					htmi <- substring(html, colorPatterns[i], colorPatterns[i + 1])
					for (j in seq_along(colors)) {
						htmi <- gsub(patterns[j],
							classes[j],
							htmi)
					}
					end <- ifelse(i==(length(colorPatterns) - 1),
						max(nchar(html)),
						colorPatterns[i + 2] - 1)
					
					htm <- paste(htm, htmi, substring(html, colorPatterns[i + 1] + 1, end), sep="")
				}
				html <- htm
			} else {
				for (j in seq_along(colors)) {
					html <- gsub(patterns[j],
						classes[j],
						html)
				}
			}
		}
		
		html <- paste(html,
			myLengths,
			"", # post-spacing
			sep="    ")
		
		# add the legend and consensus sequence
		html <- paste("", # pre-spacing
			format(c("", "", names(myXStringSet)[1:(l - 1)], "", "Consensus", "", ""), justify="right"),
				html,
			sep="    ")
		html <- gsub("&", "&amp;", html, fixed=TRUE)
		
		styles <- character()
		for (i in seq_along(colors)) {
			styles <- paste(styles,
				"span._", i,
				" {background:", colors[i],
				"; color:#FFF;} ",
				sep="")
		}
		styles <- paste("<style type=text/css> ",
			styles,
			"</style>",
			sep="")
		html <- c('<html><head><meta http-equiv="Content-Type" content="text/html; charset=utf-8"></head>',styles,"<pre>", html, "</pre></html>")
	} else {
		html <- paste(html,
			myLengths,
			"", # post-spacing
			sep="    ")
		
		# add the legend and consensus sequence
		html <- paste("", # pre-spacing
			format(c("", "", names(myXStringSet)[1:(l - 1)], "", "Consensus", "", ""), justify="right"),	
			html,
			sep="    ")
		
		html <- c("<html>","<pre>", html, "</pre></html>")
	}
	
	# replace unicode 'middle dot' with html entity
	html <- gsub("\u00B7", "&#183;", html, fixed=TRUE)
	
	writeLines(html, htmlfile)
	
	if (openURL)
		browseURL(path.expand(htmlFile))
	
	invisible(htmlFile)
}
