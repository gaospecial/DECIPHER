#' Obtain Specific Sequences from a Database
#' 
#' Returns the set of sequences meeting the search criteria.
#' 
#' If \code{type} is \code{"DNAStringSet"} then all U's are converted to T's
#' before creating the \code{DNAStringSet}, and vise-versa if \code{type} is
#' \code{"RNAStringSet"}.  All remaining characters not in the
#' \code{XStringSet}'s alphabet are converted to \code{replaceChar} or removed
#' if \code{replaceChar} is \code{""}.  Note that if \code{replaceChar} is
#' \code{NA} (the default), it will result in an error when an unexpected
#' character is found.  Quality information is interpreted as
#' \code{PredQuality} scores.
#' 
#' @name SearchDB
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table where the sequences are
#' located.
#' @param identifier Optional character string used to narrow the search
#' results to those matching a specific identifier.  If "" (the default) then
#' all identifiers are selected.
#' @param type The type of \code{XStringSet} (sequences) to return.  This
#' should be (an unambiguous abbreviation of) one of \code{"XStringSet"},
#' \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"},
#' \code{"BStringSet"}, \code{"QualityScaledXStringSet"},
#' \code{"QualityScaledDNAStringSet"}, \code{"QualityScaledRNAStringSet"},
#' \code{"QualityScaledAAStringSet"}, or \code{"QualityScaledBStringSet"}.  If
#' \code{type} is \code{"XStringSet"} or \code{"QualityScaledXStringSet"} then
#' an attempt is made to guess the type of sequences based on their
#' composition.
#' @param limit Number of results to display.  The default (\code{-1}) does not
#' limit the number of results.
#' @param replaceChar Optional character used to replace any characters of the
#' sequence that are not present in the \code{XStringSet}'s alphabet.  Not
#' applicable if \code{type=="BStringSet"}.  The default (\code{NA}) results in
#' an error if an incompatible character exist.  (See details section below.)
#' @param nameBy Character string giving the column name for naming the
#' \code{XStringSet}.
#' @param orderBy Character string giving the column name for sorting the
#' results.  Defaults to the order of entries in the database.  Optionally can
#' be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the
#' default) or descending order.
#' @param countOnly Logical specifying whether to return only the number of
#' sequences.
#' @param removeGaps Determines how gaps ("-" or "." characters) are removed in
#' the sequences.  This should be (an unambiguous abbreviation of) one of
#' \code{"none"}, \code{"all"} or \code{"common"}.
#' @param clause An optional character string to append to the query as part of
#' a ``where clause''.
#' @param quality The type of quality object to return if \code{type} is a
#' \code{QualityScaledXStringSet}.  This should be (an unambiguous abbreviation
#' of) one of \code{"Phred"}, \code{"Solexa"}, or \code{"Illumina"}.  Note that
#' recent versions of Illumina software provide \code{"Phred"} formatted
#' quality scores.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display queries as they are
#' sent to the database.
#' @return An \code{XStringSet} or \code{QualityScaledXStringSet} with the
#' sequences that meet the specified criteria.  The \code{names} of the object
#' correspond to the value in the \code{nameBy} column of the database.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{Seqs2DB}}, \code{\link{DB2Seqs}}
#' @references ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological
#' Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' # get all sequences in the default table:
#' dna <- SearchDB(db)
#' # select a random sequence:
#' dna <- SearchDB(db, orderBy="random()", limit=1)
#' # remove gaps from "Mycobacterium" sequences:
#' dna <- SearchDB(db, identifier="Mycobacterium", removeGaps="all")
#' # provide a more complex query:
#' dna <- SearchDB(db, nameBy="description", orderBy="bases", removeGaps="common",
#'                 clause="nonbases is 0")
#' 
#' @export SearchDB
SearchDB <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="XStringSet",
	limit=-1,
	replaceChar=NA,
	nameBy="row_names",
	orderBy="row_names",
	countOnly=FALSE,
	removeGaps="none",
	quality="Phred",
	clause="",
	processors=1,
	verbose=TRUE) {
	
	# error checking
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
		"QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet", "XStringSet", "QualityScaledXStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (length(identifier) != 1)
		stop("identifier must be a single character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(countOnly))
		stop("countOnly must be a logical.")
	QUALS <- c("Phred", "Solexa", "Illumina")
	quality <- pmatch(quality[1], QUALS)
	if (is.na(quality))
		stop("Invalid quality.")
	if (quality == -1)
		stop("Ambiguous quality.")
	if (quality==1) {
		quality <- PhredQuality
	} else if (quality==2) {
		quality <- SolexaQuality
	} else if (quality==3) {
		quality <- IlluminaQuality
	}
	if (!is.character(clause))
		stop("clause must be a character string.")
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
	if (type > 4 && type != 9 && removeGaps > 1)
		stop(paste('removeGaps must be "none" when type is ', TYPES[type], '.', sep=''))
	if (is.numeric(limit)) {
		if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
	} else {
		if (!grepl("[0-9],[0-9]", limit, perl=TRUE)) {
			limit <- as.numeric(limit)
			if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
		}
	}
	if (countOnly && (is.character(limit) || limit > 0))
		stop("limit cannot be specified when countOnly is TRUE.")
	
	if (verbose)
		time.1 <- Sys.time()
	
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
	
	# build the search expression
	if (countOnly) {
		searchExpression <- paste('select count(*) from ',
			tblName,
			sep="")
	} else if (nameBy=="row_names" && orderBy=="row_names") {
		searchExpression <- paste('select row_names, sequence',
			ifelse(type > 4 && type != 9, ', quality', ''),
			' from _',
			tblName,
			" where row_names in (select row_names from ",
			tblName,
			sep="")
	} else {
		searchExpression <- paste('select ',
			ifelse(nameBy=="row_names",
				paste(tblName, ".row_names", sep=""),
				nameBy),
			', _',
			tblName,
			'.sequence',
			ifelse(type > 4 && type != 9,
				paste(', _',
					tblName,
					'.quality',
					sep=""),
				''),
			' from ',
			tblName,
			' join _',
			tblName,
			' on ',
			tblName,
			'.row_names',
			' = _',
			tblName,
			'.row_names where _',
			tblName,
			'.row_names in (select row_names from ',
			tblName,
			sep="")
	}
	
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
	
	if (!countOnly)
		searchExpression <- paste(searchExpression, ")", sep="")
	
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	if (limit > 0)
		searchExpression <- paste(searchExpression,
			'limit',
			limit)
	
	if (verbose)
		cat("Search Expression:",
			strwrap(searchExpression,
					width=getOption("width") - 1L),
				sep="\n")
	
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
	dbClearResult(rs)
	
	if (countOnly) {
		count <- as.integer(searchResult)
	} else {
		# decompress the resulting sequences
		searchResult$sequence <- Codec(searchResult$sequence,
			processors=processors)
		
		if (type==9 || type==10) {
			# guess the input type of XStringSet
			freqs <- .Call("composition",
				searchResult$sequence,
				PACKAGE="DECIPHER")
			
			if (all(freqs[1:2] < 0.6)) { # not DNA/RNA
				if (freqs[3] > 0.9) { # AA
					if (type==9) {
						type <- 3
					} else {
						type <- 7
					}
				} else {
					if (type==9) {
						type <- 4
					} else {
						type <- 8
					}
				}
			} else if (freqs[1] > freqs[2]) { # DNA
				if (type==9) {
					type <- 1
				} else {
					type <- 5
				}
			} else { # RNA
				if (type==9) {
					type <- 2
				} else {
					type <- 6
				}
			}
		}
		
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
		
		if (type!=4 && type!=8) {
			# replace characters that are not in the alphabet
			searchResult$sequence <- .Call("replaceChars",
				searchResult$sequence,
				replaceChar,
				type,
				PACKAGE="DECIPHER")
		}
		
		# remove gaps if applicable
		if (removeGaps==2) {
			searchResult$sequence <- .Call("replaceChar",
				searchResult$sequence,
				"-",
				"",
				PACKAGE="DECIPHER")
			searchResult$sequence <- .Call("replaceChar",
				searchResult$sequence,
				".",
				"",
				PACKAGE="DECIPHER")
		} else if (removeGaps==3) {
			searchResult$sequence <- .Call("commonGaps",
				searchResult$sequence,
				PACKAGE="DECIPHER")
		}
		
		# build an XStringSet based on the database sequences
		if (type==1) {
			myXStringSet <- DNAStringSet(searchResult$sequence)
		} else if (type==2) {
			myXStringSet <- RNAStringSet(searchResult$sequence)
		} else if (type==3) {
			myXStringSet <- AAStringSet(searchResult$sequence)
		} else if (type==4) {
			myXStringSet <- BStringSet(searchResult$sequence)
		} else {
			w <- which(lengths(searchResult$quality)==0)
			if (length(w)==length(searchResult$quality)) {
				stop("All sequences are missing quality scores.")
			} else if (length(w) > 0) {
				stop("Some sequences are missing quality scores.")
			} else {
				searchResult$quality <- Codec(searchResult$quality)
			}
			
			if (type==5) {
				myXStringSet <- QualityScaledDNAStringSet(DNAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else if (type==6) {
				myXStringSet <- QualityScaledRNAStringSet(RNAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else if (type==7) {
				myXStringSet <- QualityScaledAAStringSet(AAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else { # type==8
				myXStringSet <- QualityScaledBStringSet(BStringSet(searchResult$sequence),
					quality(searchResult$quality))
			}
		}
		names(myXStringSet) <- searchResult[, 1]
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		if (countOnly) {
			cat("\nCount = ",
				count,
				" matches.\n",
				sep="")	
		} else {
			cat("\n",
				TYPES[type],
				" of length: ",
				length(myXStringSet),
				"\n",
				sep="")
		}
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	if (countOnly) {
		return(count)
	} else {
		return(myXStringSet)
	}
}
