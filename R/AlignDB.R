#' Align Two Sets of Aligned Sequences in a Sequence Database
#' 
#' Merges the two separate sequence alignments in a database.  The aligned
#' sequences must have separate identifiers in the same table or be located in
#' different database tables.
#' 
#' Sometimes it is useful to align two large sets of sequences, where each set
#' of sequences is already aligned but the two sets are not aligned to each
#' other.  \code{AlignDB} first builds a profile of each sequence set in
#' increments of \code{batchSize} so that the entire sequence set is not
#' required to fit in memory.  Next the two profiles are aligned using dynamic
#' programming.  Finally, the new alignment is applied to all the sequences as
#' they are incrementally added to the \code{add2tbl}.
#' 
#' Two \code{identifier}s or \code{tblName}s must be provided, indicating the
#' two sets of sequences to align.  The sequences corresponding to the first
#' \code{identifier} and \code{tblName} will be aligned to those of the second
#' \code{identifier} or \code{tblName}.  The aligned sequences are added to
#' \code{add2tbl} under a new identifier formed from the concatenation of the
#' two \code{identifier}s or \code{tblName}s.  (See examples section below.)
#' 
#' @name AlignDB
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table(s) where the sequences
#' are located.  If two \code{tblName}s are provided then the sequences in both
#' tables will be aligned.
#' @param identifier Optional character string used to narrow the search
#' results to those matching a specific identifier.  If "" then all identifiers
#' are selected.  If two \code{identifier}s are provided then the set of
#' sequences matching each \code{identifier} will be aligned.
#' @param type The type of \code{XStringSet} being processed.  This should be
#' (an abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"}, or
#' \code{"RNAStringSet"}.
#' @param add2tbl Character string specifying the table name in which to add
#' the aligned sequences.
#' @param batchSize Integer specifying the number of sequences to process at a
#' time.
#' @param perfectMatch Numeric giving the reward for aligning two matching
#' nucleotides in the alignment.  Only used when \code{type} is
#' \code{DNAStringSet} or \code{RNAStringSet}.
#' @param misMatch Numeric giving the cost for aligning two mismatched
#' nucleotides in the alignment.  Only used when \code{type} is
#' \code{DNAStringSet} or \code{RNAStringSet}.
#' @param gapOpening Numeric giving the cost for opening a gap in the
#' alignment.
#' @param gapExtension Numeric giving the cost for extending an open gap in the
#' alignment.
#' @param gapPower Numeric specifying the exponent to use in the gap cost
#' function.
#' @param terminalGap Numeric giving the cost for allowing leading and trailing
#' gaps ("-" or "." characters) in the alignment.  Either two numbers, the
#' first for leading gaps and the second for trailing gaps, or a single number
#' for both.
#' @param normPower Numeric giving the exponent that controls the degree of
#' normalization applied to scores by column occupancy.  If two numerics are
#' provided, the first controls the normalization power of terminal gaps, while
#' the second controls that of internal gaps.  A \code{normPower} of \code{0}
#' does not normalize the scores, which results in all columns of the profiles
#' being weighted equally, and is the optimal value for aligning fragmentary
#' sequences.  A \code{normPower} of \code{1} normalizes the score for aligning
#' two positions by their column occupancy (1 - fraction of gaps).  A
#' \code{normPower} greater than \code{1} more strongly discourages aligning
#' with ``gappy'' regions of the alignment.
#' @param standardize Logical determining whether scores are standardized to be
#' in units of per matching site. Standardization effectively divides the score
#' of each possible alignment by its length so that scores are relative rather
#' than absolute.
#' @param substitutionMatrix Either a substitution matrix representing the
#' substitution scores for an alignment or the name of the amino acid
#' substitution matrix to use in alignment.  The latter may be one of the
#' following: ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'',
#' ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', or
#' ``MIQS''.  The default (NULL) will use the \code{perfectMatch} and
#' \code{misMatch} penalties for DNA/RNA or \code{PFASUM50} for AA.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @param \dots Further arguments to be passed directly to \code{\link{Codec}}.
#' @return Returns the number of newly aligned sequences added to the database.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignProfiles}}, \code{\link{AlignSeqs}},
#' \code{\link{AlignTranslation}}, \code{\link{PFASUM}}
#' @references Wright, E. S. (2015). DECIPHER: harnessing local sequence
#' context to improve protein multiple sequence alignment. BMC Bioinformatics,
#' 16, 322. http://doi.org/10.1186/s12859-015-0749-z
#' 
#' Wright, E. S. (2020). RNAconTest: comparing tools for noncoding RNA multiple
#' sequence alignment based on structural consistency. RNA 2020, 26, 531-540.
#' @examples
#' 
#' gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' 
#' # Align two tables and place result into a third
#' dbConn <- dbConnect(SQLite(), ":memory:")
#' Seqs2DB(gen, "GenBank", dbConn, "Seqs1", tblName="Set1")
#' Seqs2DB(fas, "FASTA", dbConn, "Seqs2", tblName="Set2")
#' AlignDB(dbConn, tblName=c("Set1", "Set2"), add2tbl="AlignedSets")
#' l <- IdLengths(dbConn, "AlignedSets", add2tbl=TRUE)
#' BrowseDB(dbConn, tblName="AlignedSets") # all sequences have the same width
#' dbDisconnect(dbConn)
#' 
#' # Align two identifiers and place the result in the same table
#' dbConn <- dbConnect(SQLite(), ":memory:")
#' Seqs2DB(gen, "GenBank", dbConn, "Seqs1")
#' Seqs2DB(fas, "FASTA", dbConn, "Seqs2")
#' AlignDB(dbConn, identifier=c("Seqs1", "Seqs2"))
#' l <- IdLengths(dbConn, add2tbl=TRUE)
#' BrowseDB(dbConn) # note the sequences with a new identifier
#' dbDisconnect(dbConn)
#' 
#' @export AlignDB
AlignDB <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="DNAStringSet",
	add2tbl="Seqs",
	batchSize=10000,
	perfectMatch=5,
	misMatch=0,
	gapOpening=-14,
	gapExtension=-2,
	gapPower=-1,
	terminalGap=0,
	normPower=c(1, 0),
	standardize=TRUE,
	substitutionMatrix=NULL,
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!xor(length(identifier)==2, length(tblName)==2))
		stop("tblName or identifier must be length two.")
	if (length(identifier)==2 && identifier[1]==identifier[2])
		stop("The first and second identifier are the same.")
	if (length(tblName)==2 && tblName[1]==tblName[2])
		stop("The first and second tblName are the same.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(add2tbl))
		stop("add2tbl must be a character specifying a table name.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
	if (!is.numeric(normPower))
		stop("normPower must be a numeric.")
	if (any(normPower < 0))
		stop("normPower must be at least zero.")
	if (length(normPower) < 2) {
		normPower <- rep(normPower, 2)
	} else if (length(normPower) > 2) {
		stop("Length of normPower must be 1 or 2.")
	}
	normPower <- as.double(normPower)
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
	# initialize database
	driver = dbDriver("SQLite")
	if (is.character(dbFile)) {
		dbConn = dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn = dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	gapOpening <- gapOpening/2 # split into gap opening and closing
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(gapPower))
		stop("gapPower must be a numeric.")
	if (!is.numeric(terminalGap))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("Length of terminalGap must be 1 or 2.")
	if (any(is.infinite(terminalGap)))
		stop("terminalGap must be finite.")
	if (length(terminalGap)==1)
		terminalGap[2] <- terminalGap[1]
	if (!isTRUEorFALSE(standardize))
		stop("standardize must be a logical.")
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
	
	if (type==3L) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use PFASUM50
			substitutionMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
				stop("Invalid substitutionMatrix.")
		}
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment(), package=ifelse(substitutionMatrix=="MIQS", "DECIPHER", "Biostrings"))))
		}
		subMatrix <- subMatrix[AAs, AAs]
		subMatrix <- as.numeric(subMatrix)
	} else {
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				bases <- c("A", "C", "G",
					ifelse(type==2L, "U", "T"))
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
				substitutionMatrix <- as.numeric(substitutionMatrix)
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else if (type==2L && missing(perfectMatch) && missing(misMatch)) {
			substitutionMatrix <- matrix(c(10, 3, 5, 3, 3, 12, 3, 5, 5, 3, 12, 3, 3, 5, 3, 10),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
		}
	}
	
	count1 <- SearchDB(dbFile=dbConn,
		identifier=identifier[1],
		tblName=tblName[1],
		countOnly=TRUE,
		processors=processors,
		verbose=FALSE)
	if (count1==0)
		stop("No sequences found in dbFile matching the specified criteria.")
	count2 <- SearchDB(dbFile=dbConn,
		identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
		tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
		countOnly=TRUE,
		processors=processors,
		verbose=FALSE)
	if (count2==0)
		stop("No sequences found in dbFile matching the specified criteria.")
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(max=2*(count1 + count2), style=ifelse(interactive(), 3, 1))
		totSeqs <- 0L
	}
	
	for (i in 1:ceiling(count1/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[1],
			tblName=tblName[1],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			processors=processors,
			verbose=FALSE)
		
		temp <- unique(width(myXStringSet))
		if (length(temp) > 1)
			stop("Multiple width sequences found.")
		if (i==1) {
			uw1 <- temp
		} else {
			if (uw1!=temp)
				stop("Multiple width sequences found.")
		}
		
		if (type==3L) {
			temp <- .Call("consensusProfileAA",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		} else {
			temp <- .Call("consensusProfile",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		}
		
		if (i==1) {
			p.profile <- temp*length(myXStringSet)/count1
		} else {
			p.profile <- p.profile + temp*length(myXStringSet)/count1
		}
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	for (i in 1:ceiling(count2/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
			tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			processors=processors,
			verbose=FALSE)
		
		temp <- unique(width(myXStringSet))
		if (length(temp) > 1)
			stop("Multiple width sequences found.")
		if (i==1) {
			uw2 <- temp
			if (uw1*uw2 > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				uw1*uw2,
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		} else {
			if (uw2!=temp)
				stop("Multiple width sequences found.")
		}
		
		if (type==3L) {
			temp <- .Call("consensusProfileAA",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		} else {
			temp <- .Call("consensusProfile",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		}
		
		if (i==1) {
			s.profile <- temp*length(myXStringSet)/count2
		} else {
			s.profile <- s.profile + temp*length(myXStringSet)/count2
		}
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	restrict <- c(-1e10, 1e10, 1e10) # unrestricted
	if (type==3) { # AAStringSet
		inserts <- .Call("alignProfilesAA",
			p.profile,
			s.profile,
			subMatrix,
			numeric(),
			gapOpening,
			gapExtension,
			gapPower,
			normPower,
			terminalGap[1],
			terminalGap[2],
			restrict,
			standardize,
			processors,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		inserts <- .Call("alignProfiles",
			p.profile,
			s.profile,
			type,
			substitutionMatrix,
			numeric(),
			perfectMatch,
			misMatch,
			gapOpening,
			gapExtension,
			gapPower,
			normPower,
			terminalGap[1],
			terminalGap[2],
			restrict,
			standardize,
			processors,
			PACKAGE="DECIPHER")
	}
	
	if (length(identifier)==2) {
		id <- paste(identifier[1],
			identifier[2],
			sep="_")
	} else {
		id <- paste(tblName[1],
			tblName[2],
			sep="_")
	}
	
	for (i in 1:ceiling(count1/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[1],
			tblName=tblName[1],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			processors=processors,
			verbose=FALSE)
		ns <- names(myXStringSet)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[1]]),
			as.integer(inserts[[2]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(ns,
			ifelse(length(identifier)==2, identifier[1], tblName[1]),
			sep="_")
		
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			processors=processors,
			verbose=FALSE,
			...)
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	for (i in 1:ceiling(count2/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
			tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			processors=processors,
			verbose=FALSE)
		ns <- names(myXStringSet)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[3]]),
			as.integer(inserts[[4]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(ns,
			ifelse(length(identifier)==2, identifier[2], tblName[2]),
			sep="_")
		
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			processors=processors,
			verbose=FALSE,
			...)
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	if (verbose) {
		close(pBar)
		cat(strwrap(paste("\nAdded ",
					count1 + count2,
					" aligned sequences to table ",
					add2tbl,
					" with identifier '",
					id,
					"'.\n\n",
					sep="",
					collapse=""),
				width=getOption("width") - 1L),
			sep="\n")
	}
	
	invisible(count1 + count2)
}
