#' Determine the Number of Characters in Each Sequence
#' 
#' Counts the number of standard and non-standard characters in each sequence.
#' 
#' \code{IdLengths} is designed to efficiently determine the number of standard
#' and non-standard characters in every sequence within a database.  Standard
#' and non-standard characters are defined with respect to the \code{type} of
#' the sequences.  For DNA and RNA sequences there are four standard characters
#' and 11 non-standard characters (i.e., ambiguity codes).  For amino acid
#' sequences there are 20 standard and seven non-standard characters (including
#' stops).  Gap (``-''), missing (``.''), and mask (``+'') characters count
#' toward the \code{width} but not the number of standard or non-standard
#' characters.
#' 
#' @name IdLengths
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table where the sequences are
#' located.
#' @param type The type of \code{XStringSet} being processed.  This should be
#' (an abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"}, or
#' \code{"RNAStringSet"}.
#' @param add2tbl Logical or a character string specifying the table name in
#' which to add the result.
#' @param batchSize Integer specifying the number of sequences to process at a
#' time.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @return A \code{data.frame} with the number of \code{standard} characters,
#' \code{nonstandard} characters, and \code{width} of each sequence.  The
#' \code{row.names} of the \code{data.frame} correspond to the "row_names" in
#' the \code{tblName} of the \code{dbFile}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{Add2DB}}
#' @references ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological
#' Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' l <- IdLengths(db)
#' head(l)
#' 
#' @export IdLengths
IdLengths <- function(dbFile,
	tblName="Seqs",
	type="DNAStringSet",
	add2tbl=FALSE,
	batchSize=10000,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
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
	
	# initialize database
	driver <- dbDriver("SQLite")
	if (is.character(dbFile)) {
		dbConn <- dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn <- dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	if (verbose)
		time.1 <- Sys.time()
	
	if (type == 3L) { # AAStringSet
		standard <- AA_STANDARD
		nonstandard <- c("U", "O", "B", "J", "Z", "X", "*")
	} else {
		if (type == 1L) { # DNAStringSet
			standard <- DNA_BASES
		} else { # RNAStringSet
			standard <- RNA_BASES
		}
		nonstandard <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N")
	}
	
	count <- SearchDB(dbFile=dbFile,
		tblName=tblName,
		countOnly=TRUE,
		processors=processors,
		verbose=FALSE)
	
	lengths <- data.frame(standard=integer(count),
		nonstandard=integer(count),
		width=integer(count),
		names=character(count))
	
	# initialize a progress bar
	if (verbose)
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
	
	for (i in seq_len(ceiling(count/batchSize))) {
		myXStringSet <- SearchDB(dbFile=dbFile,
			tblName=tblName,
			type=TYPES[type],
			replaceChar="+",
			clause=paste("rowid > ",
				(i - 1)*batchSize,
				" and rowid <= ",
				i*batchSize,
				sep=""),
			processors=processors,
			verbose=FALSE)
		
		numF <- length(myXStringSet)
		s <- seq(from=(i - 1L)*batchSize + 1L,
			length.out=numF)
		
		alphabetTable <- alphabetFrequency(myXStringSet)
		lengths$standard[s] <- as.integer(rowSums(alphabetTable[, standard, drop=FALSE]))
		lengths$nonstandard[s] <- as.integer(rowSums(alphabetTable[, nonstandard, drop=FALSE]))
		lengths$width[s] <- width(myXStringSet)
		lengths$names[s] <- names(myXStringSet)
		
		if (verbose)
			setTxtProgressBar(pBar,
				floor(100*i/ceiling(count/batchSize)))
	}
	
	rownames(lengths) <- lengths$names
	lengths$names <- NULL
	
	if (is.character(add2tbl) || add2tbl)
		Add2DB(myData=lengths,
			dbFile=dbFile,
			tblName=ifelse(is.character(add2tbl), add2tbl, tblName),
			verbose=FALSE)
	
	if (verbose) {
		cat("\nLengths counted for ", count, " sequences.", sep="")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				":  \"standard\", \"nonstandard\", and \"width\".",
				sep="")
		cat("\n\n")
		
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(lengths)
}
