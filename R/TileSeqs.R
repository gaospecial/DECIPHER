#' Form a Set of Tiles for Each Group of Sequences.
#' 
#' Creates a set of tiles that represent each group of sequences in the
#' database for downstream applications.
#' 
#' \code{TileSeqs} will create a set of overlapping tiles representing each
#' target site in an alignment of sequences.  The most common tile permutations
#' are added until the desired minimum group coverage is obtained.  The
#' \code{dbFile} is assumed to contain \code{DNAStringSet} sequences (any U's
#' are converted to T's).
#' 
#' Target sites with one more more tiles not meeting a set of requirements are
#' marked with \code{misprime} equals \code{TRUE}.  Requirements include
#' minimum group coverage, minimum length, and maximum length.  Additionally,
#' tiles are required not to contain more than four runs of a single base or
#' four di-nucleotide repeats.
#' 
#' @name TileSeqs
#' @param dbFile A SQLite connection object or a character string specifying
#' the path to the database file.
#' @param tblName Character string specifying the table of sequences to use for
#' forming tiles.
#' @param identifier Optional character string used to narrow the search
#' results to those matching a specific identifier.  If "" then all identifiers
#' are selected.
#' @param minLength Integer providing the minimum number of nucleotides in each
#' tile.  Typically the same or slightly less than \code{maxLength}.
#' @param maxLength Integer providing the maximum number of nucleotides in each
#' tile.  Tiles are designed primarily for this length, which should ideally be
#' slightly greater than the maximum length of oligos used in downstream
#' functions.
#' @param maxTilePermutations Integer specifying the maximum number of tiles in
#' each target site.
#' @param minCoverage Numeric providing the fraction of coverage that is
#' desired for each target site in the group.  For example, a
#' \code{minCoverage} of 0.9 request that additional tiles are added until 90\%
#' of the group is represented by the tile permutations.
#' @param add2tbl Logical or a character string specifying the table name in
#' which to add the result.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @param \dots Additional arguments to be passed directly to \code{SearchDB}.
#' @return A \code{data.frame} with a row for each tile, and multiple columns
#' of information.  The \code{row_names} column gives the row number.  The
#' \code{start}, \code{end}, \code{start_aligned}, and \code{end_aligned}
#' columns provide positioning of the tile in a consensus sequence formed from
#' the group.  The column \code{misprime} is a logical specifying whether the
#' tile meets the specified constraints.  The columns \code{width} and
#' \code{id} indicate the tile's length and group of origin, respectively.
#' 
#' The \code{coverage} field gives the fraction of sequences containing the
#' tile in the group that encompass the tile's start and end positions in the
#' alignment, whereas \code{groupCoverage} contains the fraction of all
#' sequences in the group containing a tile at their respective target site.
#' For example, if only a single sequence out of 10 has information (no gap) in
#' the first alignment position, then \code{coverage} would be 100\% (1.0),
#' while \code{groupCoverage} would be 10\% (0.1).
#' 
#' The final column, \code{target_site}, provides the sequence of the tile.
#' @note If \code{add2tbl} is \code{TRUE} then the tiles will be added to the
#' database table that currently contains the sequences used for tiling.  The
#' added tiles may cause interference when querying a table of sequences.
#' Therefore, it is recommended to add the tiles to their own table, for
#' example, by using \code{add2tbl="Tiles"}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{DesignPrimers}}
#' @examples
#' 
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' tiles <- TileSeqs(db, identifier="Pseudomonadales")
#' 
#' @export TileSeqs
TileSeqs <- function(dbFile,
	tblName="Seqs",
	identifier="",
	minLength=26,
	maxLength=27,
	maxTilePermutations=10,
	minCoverage=0.9,
	add2tbl=FALSE,
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
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
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(maxTilePermutations))
		stop("maxTilePermutations must be a numeric.")
	if (floor(maxTilePermutations)!=maxTilePermutations)
		stop("maxTilePermutations must be a whole number.")
	if (maxTilePermutations < 1)
		stop("maxTilePermutations must be at least 1.")
	if (!is.numeric(minCoverage))
		stop("minCoverage must be a numeric.")
	if (minCoverage > 1 || minCoverage < 0)
		stop("minCoverage must be between zero and one.")
	if (!is.numeric(minLength))
		stop("minLength must be a numeric.")
	if (floor(minLength)!=minLength)
		stop("minLength must be a whole number.")
	if (minLength < 1)
		stop("minLength must be at least 1.")
	if (!is.numeric(maxLength))
		stop("maxLength must be a numeric.")
	if (floor(maxLength)!=maxLength)
		stop("maxLength must be a whole number.")
	if (maxLength < 1)
		stop("maxLength must be at least 1.")
	if (minLength > maxLength)
		stop("minLength must be less than or equal to maxLength.")
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
	
	searchExpression <- paste("select distinct identifier from",
		tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
	ids <- searchResult$identifier
	dbClearResult(rs)
	
	if (identifier[1]=="") {
		identifier <- ids
	} else {
		w <- which(!(identifier %in% ids))
		if (length(w) > 0)
			stop("identifier not in tiles: ",
				paste(identifier[w], collapse=", "))
	}
	
	if (is.character(add2tbl) || add2tbl) {
		result <- dbListTables(dbConn)
		w <- which(result==ifelse(is.character(add2tbl), add2tbl, tblName))
		if (length(w)==1) { # add to existing table
			searchExpression <- paste("select max(row_names) from ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				sep="")
			row_start <- as.integer(dbGetQuery(dbConn, searchExpression))
		} else { # create new table
			row_start <- 0
		}
	} else { # don't add to table
		row_start <- 0
	}
	
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	}
	count <- 0
	for (k in seq_along(identifier)) {
		target <- SearchDB(dbConn,
			tblName=tblName,
			type="DNAStringSet",
			identifier=identifier[k],
			processors=processors,
			verbose=FALSE,
			...)
		numF <- length(target)
		
		if (numF==0) {
			warning("No sequences: ", identifier[k])
			next
		}
		
		uw <- unique(width(target))
		if (length(uw) > 1) {
			warning("Skipped due to multiple width sequences: ", identifier[k])
			next
		}
		
		a <- alphabetFrequency(target)
		if (length(target) > 1 && all(rowSums(a[,1:15]) < maxLength)) {
			warning("Skipped because sequences shorter than maxLength: ", identifier[k])
			next
		} else if (length(target)==1 && sum(a[,1:15]) < maxLength) {
			warning("Skipped because sequence shorter than maxLength: ", identifier[k])
			next
		}
		
		consensus <- ConsensusSequence(target)
		
		pos <- which(strsplit(toString(consensus), "", fixed=TRUE)[[1]] != "-")
		l <- length(pos) - maxLength + 1
		
		tiles <- data.frame(row_names=(row_start + count + 1):(row_start + count + l*maxTilePermutations),
			start=I(integer(l*maxTilePermutations)),
			end=I(integer(l*maxTilePermutations)),
			start_aligned=I(integer(l*maxTilePermutations)),
			end_aligned=I(integer(l*maxTilePermutations)),
			misprime=I(logical(l*maxTilePermutations)),
			width=I(rep(uw, l*maxTilePermutations)),
			id=I(rep(identifier[k], l*maxTilePermutations)),
			coverage=I(numeric(l*maxTilePermutations)),
			groupCoverage=I(numeric(l*maxTilePermutations)),
			target_site=I(character(l*maxTilePermutations)))
		
		start <- integer(l*maxTilePermutations)
		end <- integer(l*maxTilePermutations)
		start_aligned <- integer(l*maxTilePermutations)
		end_aligned <- integer(l*maxTilePermutations)
		misprimes <- logical(l*maxTilePermutations)
		coverage <- numeric(l*maxTilePermutations)
		groupCoverage <- numeric(l*maxTilePermutations)
		target_sites <- character(l*maxTilePermutations)
		
		tGaps <- TerminalChar(target)
		
		count <- 0
		for (i in seq_len(l)) {
			# find all target_sites within terminal gaps
			w <- which(pos[i] > tGaps[,1] &
				pos[i + maxLength - 1] <= tGaps[,1] + tGaps[,3])
			if (length(w) == 0)
				next
			target_site <- subseq(target[w],
				start=pos[i],
				end=pos[i + maxLength - 1])
			target_site <- as.character(target_site)
			target_site <- .Call("replaceChar",
				target_site,
				"-",
				"",
				PACKAGE="DECIPHER")
			target_site <- .Call("replaceChar",
				target_site,
				".",
				"",
				PACKAGE="DECIPHER")
			if (all(target_site=="")) # only gaps
				next
			
			# tablulate the target_sites
			t <- table(target_site)
			w <- which(names(t)=="")
			if (length(w) > 0)
				t <- t[-w]
			t <- t[order(t, decreasing=TRUE)]
			
			# choose the top target sites
			thresh <- minCoverage*sum(t)
			j <- ifelse(length(t) > maxTilePermutations, maxTilePermutations, length(t))
			w <- which(cumsum(t[1:j]) > thresh)
			if (length(w) != 0) {
				index <- 1:min(w[1],j)
			} else {
				index <- 1:j
			}
			coverages <- as.integer(t[index])/sum(t)
			groupCoverages <- as.integer(t[index])/numF
			target_site <- names(t[index])
			
			misprime <- FALSE
			n <- nchar(target_site)
			w <- which(n < minLength | n > maxLength)
			if (length(w) > 0) {
				target_site <- target_site[-w]
				coverages <- coverages[-w]
				groupCoverages <- groupCoverages[-w]
				misprime <- TRUE
			}
			w <- which(grepl("[^A|C|T|G]", target_site))
			if (length(w) > 0) {
				target_site <- target_site[-w]
				coverages <- coverages[-w]
				groupCoverages <- groupCoverages[-w]
				misprime <- TRUE
			}
			if (length(target_site)==0)
				next
			
			index <- 1:length(target_site)
			coverage[count + index] <- coverages
			groupCoverage[count + index] <- groupCoverages
			start_aligned[count + index] <- pos[i]
			end_aligned[count + index] <- pos[i + maxLength - 1]
			start[count + index] <- i
			end[count + index] <- i + maxLength - 1
			target_sites[count + index] <- target_site
			
			for (j in seq_along(target_site)) {
				count <- count + 1
				ts <- strsplit(target_site[j], "", fixed=FALSE)[[1]]
				repeats <- 0
				runs <- 0
				
				if (!misprime) {
					for (p in 2:length(ts)) {
						if (p > 3) {
							if ((ts[p]==ts[p - 2]) && (ts[p - 1]==ts[p - 3])) {
								repeats <- repeats + 1
								if (repeats > 5) { # more than 4 di-nucleotides
									misprime <- TRUE
									break
								}
							} else {
								repeats <- 0
							}
						}
						
						if (ts[p]==ts[p - 1]) {
							runs <- runs + 1
							if (runs > 3) { # more than 4 of the same base
								misprime <- TRUE
								break
							}
						} else {
								runs <- 0
						}
					}
				}
			}
			misprimes[(count - j + 1):count] <- misprime
		}
		
		tiles$start <- start
		tiles$end <- end
		tiles$start_aligned <- start_aligned
		tiles$end_aligned <- end_aligned
		tiles$misprime <- misprimes
		tiles$coverage <- coverage
		tiles$groupCoverage <- groupCoverage
		tiles$target_site <- target_sites
		
		w <- which(tiles$target_site=="")
		if (length(w) > 0)
			tiles <- tiles[-w,]
		if (exists("tiles_all")) {
			tiles_all <- rbind(tiles_all, tiles)
		} else {
			tiles_all <- tiles
		}
		
		if (is.character(add2tbl) || add2tbl) {
			if (count==0 && row_start==0) {
				ft <- list(row_names="INTEGER PRIMARY KEY ASC",
					start="INTEGER",
					end="INTEGER",
					start_aligned="INTEGER",
					end_aligned="INTEGER",
					misprime="LOGICAL",
					width="INTEGER",
					id="TEXT",
					coverage="REAL",
					groupCoverage="REAL",
					target_site="TEXT")
				append <- FALSE
			} else {
				ft <- NULL
				append <- TRUE
			}
			dbWriteTable(dbConn,
				ifelse(is.character(add2tbl),
					add2tbl,
					tblName),
				tiles,
				row.names=FALSE,
				overwrite=FALSE,
				append=append,
				field.types=ft)
		}
		count <- dim(tiles_all)[1]
		
		if (verbose)
			setTxtProgressBar(pBar,
				floor(100*k/length(identifier)))
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		close(pBar)
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(tiles_all)
}
