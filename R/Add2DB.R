# below function from RSQLite package version 0.11.4
.sqliteDataType <- function(obj, ...) {
	rs.class <- data.class(obj)
	rs.mode <- storage.mode(obj)
	switch(rs.class,
		numeric = if (rs.mode=="integer") "INTEGER" else "REAL",
		character = "TEXT",
		logical = "INTEGER",
		factor = "TEXT",
		ordered = "TEXT",
		## list maps to BLOB. Although not checked, the list must
		## either be empty or contain only raw vectors or NULLs.
		list = "BLOB",
		## attempt to store obj according to its storage mode if it has
		## an unrecognized class.
		switch(rs.mode,
			integer = "INTEGER",
			double = "REAL",
			## you'll get this if class is AsIs for a list column
			## within a data.frame
			list = if (rs.class == "AsIs") "BLOB" else "TEXT",
			"TEXT"))
}



















#' Add Data to a Database
#' 
#' Adds a \code{data.frame} to a database table by its \code{row.names}.
#' 
#' Data contained in \code{myData} will be added to the \code{tblName} by its
#' respective \code{row.names}.
#' 
#' @name Add2DB
#' @param myData Data frame containing information to be added to the
#' \code{dbFile}.
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table in which to add the
#' data.
#' @param clause An optional character string to append to the query as part of
#' a ``where clause''.
#' @param verbose Logical indicating whether to display each query as it is
#' sent to the database.
#' @return Returns \code{TRUE} if the data was added successfully, or
#' \code{FALSE} otherwise.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{Seqs2DB}}, \code{\link{SearchDB}},
#' \code{\link{BrowseDB}}
#' @references ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological
#' Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
#' @examples
#' 
#' \dontrun{
#' # Create a sequence database
#' gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
#' dbConn <- dbConnect(SQLite(), ":memory:")
#' Seqs2DB(gen, "GenBank", dbConn, "Bacteria")
#' 
#' # Identify the sequence lengths
#' l <- IdLengths(dbConn)
#' 
#' # Add lengths to the database
#' Add2DB(l, dbConn)
#' 
#' # View the added lengths
#' BrowseDB(dbConn)
#' 
#' # Change the value of existing columns
#' ids <- data.frame(identifier=rep("Bacteroidetes", 18), stringsAsFactors=FALSE)
#' rownames(ids) <- 10:27
#' Add2DB(ids, dbConn)
#' BrowseDB(dbConn)
#' 
#' # Add data to a subset of rows using a clause
#' ids[[1]][] <- "Changed"
#' nrow(ids) # 18 rows
#' Add2DB(ids, dbConn, clause="accession like 'EU808318%'")
#' BrowseDB(dbConn) # only 1 row effected
#' 
#' dbDisconnect(dbConn)
#' }
#' @export Add2DB
Add2DB <- function(myData,
	dbFile,
	tblName="Seqs",
	clause="",
	verbose=TRUE) {
	
	time.1 <- Sys.time()
	
	# error checking
	if (!is.data.frame(myData))
		stop("myData must be a data frame object.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (substr(tblName, 1, 1) == "_")
		stop("Invalid tblName.")
	if (tblName == "taxa")
		stop("taxa is a reserved tblName.")
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (ncol(myData)==0)
		stop("myData contains no columns.")
	if (any(grepl(".", names(myData), fixed=TRUE)))
		stop("Column names cannot contain periods.")
	
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
	
	result  <- dbListTables(dbConn)
	w <- which(result==tblName)
	if (length(w)==0) { # need to make the table
		dbWriteTable(dbConn,
			tblName,
			myData,
			row.names=TRUE,
			overwrite=TRUE,
			append=FALSE)
		if (verbose)
			cat("Created table '", tblName, "'.\n", sep="")
	} else {
		result <- dbListFields(dbConn, tblName)
		colIDs <- names(myData)
		
		if (is.na(match("row_names", colIDs)))
			myData$row_names=row.names(myData)
		
		x <- sqliteQuickColumn(dbConn, tblName, "row_names")
		m <- match(myData$row_names, x)
		if (any(is.na(m)))
			stop("row.names of myData are missing from '", tblName, "'.")
		
		# loop through each of the columns to add
		for (i in 1:length(colIDs)) {
			colName <- colIDs[i]
			
			if (is.na(match(colName, result))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table ",
					tblName,
					" add column ",
					colName,
					" ",
					.sqliteDataType(myData[colIDs[i], 1]),
					sep="")
				if (verbose)
					cat("Expression:\n",
						paste(strwrap(expression1,
								width=getOption("width") - 1L),
							collapse="\n"),
						"\n\n",
						sep="")
				rs <- dbSendStatement(dbConn, expression1)
				dbClearResult(rs)
			}
			
			# next update the column with new data
			expression2 <- paste("update ",
				tblName,
				" set ",
				colName,
				" = :",
				colName,
				" where row_names = :row_names",
				sep="")
			if (clause!="")
				expression2 <- paste(expression2,
					" and ",
					clause,
					sep="")
			if (verbose)
				cat("Expression:\n",
					paste(strwrap(expression2,
							width=getOption("width") - 1L),
						collapse="\n"),
					"\n\n",
					sep="")
			rs <- dbSendQuery(dbConn, expression2)
			dbBind(rs, myData[, c("row_names", colIDs[i])])
			dbClearResult(rs)
		}
		if (verbose) {
			cat("Added to table ",
				tblName,
				":  \"",
				paste(names(myData)[-match("row_names",
					names(myData))],
				collapse="\" and \""),
				"\".\n\n",
				sep="")
		}
	}
	
	if (verbose) { # print the elapsed time to update table
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(TRUE)
}
