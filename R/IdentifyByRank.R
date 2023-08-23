#' Identify By Taxonomic Rank
#' 
#' Identifies sequences by a specific level of their taxonomic rank.
#' 
#' \code{IdentifyByRank} simply identifies a sequence by a specific level of
#' its taxonomic rank.  Requires that \code{rank} information be present in the
#' \code{tblName}, such as that created by default when importing sequences
#' from a GenBank formatted file.
#' 
#' The input parameter \code{level} should be an integer giving the ``level''
#' of the taxonomic rank to choose as the identifier.  Negative \code{level}s
#' are interpreted as being that many levels from the last level in each rank.
#' The \code{level} zero selects the base level (see below).
#' 
#' If the specified level of rank does not exist then the closest rank is
#' chosen.  Therefore, setting \code{level} to \code{Inf} will always select
#' the last taxonomic level (i.e., genus).
#' 
#' For example, a representative ``rank'' imported from a GenBank file is:\cr
#' Saccharomyces cerevisiae\cr Eukaryota; Fungi; Ascomycota; Saccharomycotina;
#' Saccharomycetes;\cr Saccharomycetales; Saccharomycetaceae; Saccharomyces.
#' 
#' Setting \code{level} to \code{0} would result in an \code{identifier} of
#' ``Saccharomyces cerevisiae'', because it is on the first line.  A
#' \code{level} of \code{2} would return ``Fungi'', and \code{-2} (second to
#' last) would return ``Saccharomycetaceae''.  A \code{level} of \code{Inf}
#' would find the nearest level to the end, ``Saccharomyces''.
#' 
#' @name IdentifyByRank
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table where the rank
#' information is located.
#' @param level Level of the taxonomic rank.  (See details section below.)
#' @param add2tbl Logical or a character string specifying the table name in
#' which to add the result.
#' @param verbose Logical indicating whether to print database queries and
#' other information.
#' @return A \code{data.frame} with the \code{rank} and corresponding
#' identifier as \code{identifier}.  Note that quotes are stripped from
#' identifiers to prevent problems that they may cause.  The \code{origin}
#' gives the \code{rank} preceding the \code{identifier}.  If \code{add2tbl} is
#' not \code{FALSE} then the ``identifier'' column is updated in \code{dbFile}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{FormGroups}}
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' ids <- IdentifyByRank(db, level=Inf)
#' head(ids)
#' 
#' @export IdentifyByRank
IdentifyByRank <- function(dbFile,
	tblName="Seqs",
	level=0,
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking:
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(level) || floor(level) != level)
		stop("level must be an integer.")
	
	if (verbose)
		time.1 <- Sys.time()
	
	# initialize databases
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
	
	if (is.na(match("rank",
		dbListFields(dbConn,
			tblName))))
		stop("No rank column in table.")
	
	searchExpression <- paste("select distinct rank from", tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	x <- dbFetch(rs, n=-1, row.names=FALSE)
	dbClearResult(rs)
	
	.change <- function(id) {
		id <- .Call("replaceChar", id, '"', "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, "'", "", PACKAGE="DECIPHER")
		id <- gsub("^\\s+|\\s+$", "", id) # trim flanking white space
		id <- gsub("\\.+$", "", id)
		return(id)
	}
	
	z <- x
	if (level==0) {
		x <- strsplit(x$rank, "\n", fixed=TRUE)
		z$origin <- unlist(lapply(x,
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		z$identifier <- .change(unlist(lapply(x, `[`, 1L)))
	} else {
		x$rank <- unlist(lapply(strsplit(x$rank,
				"\n",
				fixed=TRUE),
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		for (j in seq_along(x$rank)) {
			a <- strsplit(x$rank[j], ";")[[1]]
			l <- length(a)
			if (level < 0) {
				temp_level <- l + level + 1L
			} else {
				temp_level <- level
			}
			
			if (temp_level > l) {
				id <- as.character(a[l])
			} else if (temp_level < 1) {
				id <- as.character(a[1])
			} else {
				id <- as.character(a[temp_level])
			}
			z$origin[j] <- unlist(strsplit(as.character(x$rank[j]),
				id,
				fixed=TRUE))[1]
			
			z$identifier[j] <- .change(id)
		}
	}
	
	if (is.character(add2tbl) || add2tbl) {
		dbWriteTable(dbConn, "taxa", z)
		
		if (verbose)
			cat("\nUpdating column: \"identifier\"...")
		
		searchExpression <- paste("update ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			" set identifier = (select identifier from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank)",
			sep="")
		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)
		
		searchExpression <- "drop table taxa"
		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)
	}
	
	if (verbose) {
		cat("\nFormed",
			length(unique(z$identifier)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				": \"identifier\".",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	names(z)[1] <- "rank"
	return(z)
}
