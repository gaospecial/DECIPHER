#' Create Consensus Sequences by Groups
#' 
#' Forms a consensus sequence representing the sequences in each group.
#' 
#' Creates a consensus sequence for each of the distinct groups defined in
#' \code{colName}.  The resulting \code{XStringSet} contains as many consensus
#' sequences as there are distinct groups in \code{colName}.  For example, it
#' is possible to create a set of consensus sequences with one consensus
#' sequence for each \code{"id"} in the \code{tblName}.
#' 
#' @name IdConsensus
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table in which to form
#' consensus.
#' @param identifier Optional character string used to narrow the search
#' results to those matching a specific identifier.  If "" then all identifiers
#' are selected.
#' @param type The type of \code{XStringSet} (sequences) to use in forming
#' consensus.  This should be (an abbreviation of) one of
#' \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, or
#' \code{"BStringSet"}.
#' @param colName Column containing the group name of each sequence.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @param \dots Additional arguments to be passed directly to
#' \code{\link{ConsensusSequence}} for an \code{AAStringSet},
#' \code{DNAStringSet}, or \code{RNAStringSet}, or to \code{consensusString}
#' for a \code{BStringSet}.
#' @return An \code{XStringSet} object containing the consensus sequence for
#' each group.  The \code{names} of the \code{XStringSet} contain the number of
#' sequences and name of each group.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{Seqs2DB}}
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' con <- IdConsensus(db, colName="identifier", noConsensusChar="N")
#' BrowseSeqs(con)
#' 
#' @export IdConsensus
IdConsensus <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="DNAStringSet",
	colName="identifier",
	processors=1,
	verbose=TRUE,
	...) {
	
	# initialize variables
	if (verbose)
		time.1 <- Sys.time()
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(colName))
		stop("colName must be a character string.")
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
	
	f <- dbListFields(dbConn, tblName)
	w <- which(f==colName)
	if (length(w)==0)
		stop(paste("The colName '", colName, "' does not exist.", sep=""))
	
	searchExpression <- paste('select distinct ',
		colName,
		' from ',
		tblName,
		sep="")
	if (identifier != "")
		searchExpression <- paste(searchExpression,
			" where identifier is '",
			identifier,
			"'",
			sep="")
	rs <- dbSendQuery(dbConn, searchExpression)
	groups <- dbFetch(rs, n=-1, row.names=FALSE)[,eval(colName)]
	dbClearResult(rs)
	
	# remove any null groups in database
	w <- which(is.na(groups))
	if (length(w) > 0)
		groups <- groups[-w]
	
	if (length(groups) < 1)
		stop("No groups in which to form consensus.")
	
	# initialize a progress bar
	if (verbose)
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
	
	if (type==1) {
		consensus <- DNAStringSet()
	} else if (type==2) {
		consensus <- RNAStringSet()
	} else if (type==3) {
		consensus <- AAStringSet()
	} else { # type==4
		consensus <- BStringSet()
	}
	
	seqCount <- numeric(length(groups))
	j <- 0L
	for (i in groups) {
		j <- j + 1L
		x_subset <- SearchDB(dbFile,
			tblName=tblName,
			type=TYPES[type],
			verbose=FALSE,
			identifier=identifier,
			processors=processors,
			clause=paste(colName,
				"= '",
				gsub("'", "''", i, fixed=TRUE),
				"'",
				sep=""))
		
		if (length(consensus)==0) {
			if (type!=4) {
				consensus <- ConsensusSequence(myXStringSet=x_subset,
					...)
			} else {
				consensus <- consensusString(x=x_subset,
					...)
			}
		} else {
			if (type!=4) {
				consensus <- c(consensus,
					 ConsensusSequence(myXStringSet=x_subset,
						...))
			} else {
				consensus <- c(consensus,
					 consensusString(x=x_subset,
						...))
			}
		}
		
		seqCount[j] <- length(x_subset)
		if (verbose)
			setTxtProgressBar(pBar, 100*j/length(groups))
	}
	
	names(consensus) <- groups
	
	if (verbose) {
		close(pBar)
		
		cat("\nFound consensus for ",
			length(groups),
			" groups.",
			sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(consensus)
}
