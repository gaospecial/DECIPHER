#' Export Database Sequences to a FASTA or FASTQ File
#' 
#' Exports a database containing sequences to a FASTA or FASTQ formatted file
#' of sequence records.
#' 
#' Sequences are exported into either a FASTA or FASTQ file as determined by
#' the \code{type} of sequences.  If \code{type} is an \code{XStringSet} then
#' sequences are exported to FASTA format.  Quality information for
#' \code{QualityScaledXStringSet}s are interpreted as \code{PredQuality} scores
#' before export to FASTQ format.
#' 
#' If \code{type} is \code{"BStringSet"} (the default) then sequences are
#' exported to a FASTA file exactly the same as they were when imported.  If
#' \code{type} is \code{"DNAStringSet"} then all U's are converted to T's
#' before export, and vise-versa if \code{type} is \code{"RNAStringSet"}.  All
#' remaining characters not in the \code{XStringSet}'s alphabet are converted
#' to \code{replaceChar} or removed if \code{replaceChar} is \code{""}.  Note
#' that if \code{replaceChar} is \code{NA} (the default), it will result in an
#' error when an unexpected character is found.
#' 
#' @name DB2Seqs
#' @param file Character string giving the location where the file should be
#' written.
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table in which to extract the
#' data.
#' @param identifier Optional character string used to narrow the search
#' results to those matching a specific identifier.  If "" then all identifiers
#' are selected.
#' @param type The type of \code{XStringSet} (sequences) to export to a FASTA
#' formatted file or \code{QualityScaledXStringSet} to export to a FASTQ
#' formatted file.  This should be (an unambiguous abbreviation of) one of
#' \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"},
#' \code{"BStringSet"}, \code{"QualityScaledDNAStringSet"},
#' \code{"QualityScaledRNAStringSet"}, \code{"QualityScaledAAStringSet"}, or
#' \code{"QualityScaledBStringSet"}.  (See details section below.)
#' @param limit Number of results to display.  The default (\code{-1}) does not
#' limit the number of results.
#' @param replaceChar Optional character used to replace any characters of the
#' sequence that are not present in the \code{XStringSet}'s alphabet.  Not
#' applicable if \code{type=="BStringSet"}.  The default (\code{NA}) results in
#' an error if an incompatible character exist.  (See details section below.)
#' @param nameBy Character string giving the column name(s) for identifying
#' each sequence record.  If more than one column name is provided, the
#' information in each column is concatenated, separated by \code{sep}, in the
#' order specified.
#' @param orderBy Character string giving the column name for sorting the
#' results.  Defaults to the order of entries in the database.  Optionally can
#' be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the
#' default) or descending order.
#' @param removeGaps Determines how gaps ("-" or "." characters) are removed in
#' the sequences.  This should be (an unambiguous abbreviation of) one of
#' \code{"none"}, \code{"all"} or \code{"common"}.
#' @param append Logical indicating whether to append the output to the
#' existing \code{file}.
#' @param width Integer specifying the maximum number of characters per line of
#' sequence.  Not applicable when exporting to a FASTQ formatted file.
#' @param compress Logical specifying whether to compress the output file using
#' gzip compression.
#' @param chunkSize Number of sequences to write to the \code{file} at a time.
#' Cannot be less than the total number of sequences if \code{removeGaps} is
#' \code{"common"}.
#' @param sep Character string providing the separator between fields in each
#' sequence's name, by default pairs of colons (``::'').
#' @param clause An optional character string to append to the query as part of
#' a ``where clause''.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display status.
#' @return Writes a FASTA or FASTQ formatted file containing the sequence
#' records in the database.
#' 
#' Returns the number of sequence records written to the \code{file}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @references ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological
#' Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' tf <- tempfile()
#' DB2Seqs(tf, db, limit=10)
#' file.show(tf) # press 'q' to exit
#' unlink(tf)
#' 
#' @export DB2Seqs
DB2Seqs <- function(file,
	dbFile,
	tblName="Seqs",
	identifier="",
	type="BStringSet",
	limit=-1,
	replaceChar=NA,
	nameBy="description",
	orderBy="row_names",
	removeGaps="none",
	append=FALSE,
	width=80,
	compress=FALSE,
	chunkSize=1e5,
	sep="::",
	clause="",
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(file))
		stop("file must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
		"QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(limit))
		stop("limit must be a numeric.")
	if (floor(limit)!=limit)
		stop("limit must be a whole number.")
	if (!is.character(nameBy))
		stop("nameBy must be a character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(append))
		stop("append must be a logical.")
	if (!is.logical(compress))
		stop("compress must be a logical.")
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(chunkSize))
		stop("chunkSize must be a numeric.")
	if (floor(chunkSize)!=chunkSize)
		stop("chunkSize must be a whole number.")
	if (chunkSize <= 0)
		stop("chunkSize must be greater than zero.")
	if (!is.numeric(width))
		stop("width must be a numeric.")
	if (floor(width)!=width)
		stop("width must be a whole number.")
	if (width < 1)
		stop("width must be at least 1.")
	if (width > 20001)
		stop("width can be at most 20001.")
	if (!is.character(sep))
		stop("sep must be a character string.")
	if (length(sep)!=1)
		stop("sep must be a single character string.")
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps, GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	if (is.na(replaceChar)) {
			replaceChar <- NA_character_
	} else if (type==1 || type==5) {
		if (is.na(pmatch(replaceChar, DNA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the DNA_ALPHABET or empty character.")
	} else if (type==2 || type==6) {
		if (is.na(pmatch(replaceChar, RNA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the RNA_ALPHABET or empty character.")
	} else if (type==3 || type==7) {
		if (is.na(pmatch(replaceChar, AA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the AA_ALPHABET or empty character.")
	}
	if (type > 4 && removeGaps > 1)
		stop(paste('removeGaps must be "none" when type is ', TYPES[type], '.', sep=''))
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	searchExpression <- tblName
	
	if (identifier!="")
		searchExpression <- paste(searchExpression,
			' where identifier is "',
			identifier,
			'"',
			sep="")
	if (clause!="")
		searchExpression <- paste(searchExpression,
			ifelse(identifier=="", " where ", " and "),
			clause,
			sep="")
	
	if (limit > 0) {
		searchExpression1 <- paste(searchExpression,
			'limit',
			limit)
	} else {
		searchExpression1 <- searchExpression
	}
	searchExpression1 <- paste('select count(*) from ',
		searchExpression1,
		sep="")
	rs <- dbSendQuery(dbConn, searchExpression1)
	count <- as.numeric(dbFetch(rs, n=-1, row.names=FALSE))
	dbClearResult(rs)
	
	if (count < 1)
		stop("No sequences matched the specified parameters.")
	if (limit!=-1 && count > limit)
		count <- limit
	if (count > chunkSize && removeGaps==3)
		stop("chunkSize must be at least ",
			count,
			" with this query if removeGaps is 'common'.")
	
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	
	if (verbose)
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	s <- seq(1, count, chunkSize)
	for (i in 1:length(s)) {
		# build the search expression
		searchExpression1 <- paste(searchExpression,
			'limit',
			ifelse(i==length(s), count - s[i] + 1, chunkSize),
			'offset',
			ifelse(i==1, 0, s[i] - 1))
		if (all(nameBy=="row_names")) {
			searchExpression1 <- paste("select row_names from",
				searchExpression1)
		} else {
			temp <- 'select row_names'
			for (j in 1:length(nameBy)) {
				if (nameBy[j]=="row_names")
					next
				temp <- paste(temp, ", ", nameBy[j], sep="")
			}
			searchExpression1 <- paste(temp,
				"from",
				searchExpression1)
		}
		
		rs <- dbSendQuery(dbConn, searchExpression1)
		searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
		dbClearResult(rs)
		
		searchExpression2 <- paste(ifelse(type > 4,
				'select row_names, sequence, quality from _',
				'select row_names, sequence from _'),
			tblName,
			" where row_names in (select row_names from ",
			searchExpression,
			")",
			sep="")
		rs <- dbSendQuery(dbConn, searchExpression2)
		searchResult2 <- dbFetch(rs, n=-1, row.names=FALSE)
		dbClearResult(rs)
		
		m <- match(searchResult$row_names,
			searchResult2$row_names)
		searchResult2 <- searchResult2[m,]
		
		# decompress the resulting sequences
		searchResult2$sequence <- Codec(searchResult2$sequence,
			processors=processors)
		
		if (type!=4 && type!=8) {
			# replace characters that are not in the DNA_ALPHABET
			searchResult2$sequence <- .Call("replaceChars",
				searchResult2$sequence,
				replaceChar,
				type,
				PACKAGE="DECIPHER")
		}
		
		# remove gaps if applicable
		if (removeGaps==2) {
			searchResult2$sequence <- .Call("replaceChar",
				searchResult2$sequence,
				"-",
				"",
				PACKAGE="DECIPHER")
			searchResult2$sequence <- .Call("replaceChar",
				searchResult2$sequence,
				".",
				"",
				PACKAGE="DECIPHER")
		} else if (removeGaps==3) {
			searchResult2$sequence <- .Call("commonGaps",
				searchResult2$sequence,
				PACKAGE="DECIPHER")
		}
		
		if (type > 4) {
			# decompress the resulting qualities
			searchResult2$quality <- Codec(searchResult2$quality,
				processors=processors)
		}
		
		if (type==1 || type==5) {
			myXStringSet <- DNAStringSet(searchResult2$sequence)
		} else if (type==2 || type==6) {
			myXStringSet <- RNAStringSet(searchResult2$sequence)
		} else if (type==3 || type==7) {
			myXStringSet <- AAStringSet(searchResult2$sequence)
		} else if (type==4 || type==8) {
			myXStringSet <- BStringSet(searchResult2$sequence)
		}
		
		if (length(nameBy) > 1) {
			names(myXStringSet) <- do.call(paste,
				c(searchResult[, nameBy],
					sep=sep))
		} else {
			names(myXStringSet) <- searchResult[, nameBy]
		}
		
		if (type > 4) {
			writeXStringSet(myXStringSet,
				file,
				append=ifelse(i==1, append, TRUE),
				format="FASTQ",
				compress=compress,
				qualities=PhredQuality(searchResult2$quality))
		} else {
			writeXStringSet(myXStringSet,
				file,
				append=ifelse(i==1, append, TRUE),
				format="FASTA",
				compress=compress,
				width=width)
		}
		
		if (verbose)
			setTxtProgressBar(pBar, ifelse(i==length(s), 1, (s[i + 1] - 1)/count))
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n\nWrote ",
			count,
			ifelse(count > 1, " sequences.", " sequence."),
			"\n",
			sep="")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(count)
}
