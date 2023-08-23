#' Write Genes to a File
#' 
#' Writes predicted genes to a file in GenBank (gbk) or general feature format
#' (gff).
#' 
#' \code{WriteGenes} will write a \code{"Genes"} object to a GenBank (if
#' \code{format} is \code{"gbk"}) or general feature format (if \code{format}
#' is \code{"gff"}) file.
#' 
#' @name WriteGenes
#' @param x An object of class \code{Genes}.
#' @param file A connection or a character string naming the file path where
#' the tree should be exported.  If "" (the default), the tree is printed to
#' the standard output connection, the console unless redirected by sink.
#' @param format Character specifying \code{"gbk"} or \code{"gff"} output
#' format.
#' @param append Logical indicating whether to append to an existing
#' \code{file}.  Only applicable if \code{file} is a character string.  If
#' \code{FALSE} (the default), then the file is overwritten.
#' @return \code{NULL}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{ExtractGenes}}, \code{\link{FindGenes}},
#' \code{\link{Genes-class}}
#' @examples
#' 
#' # import a test genome
#' fas <- system.file("extdata",
#' 	"Chlamydia_trachomatis_NC_000117.fas.gz",
#' 	package="DECIPHER")
#' genome <- readDNAStringSet(fas)
#' 
#' x <- FindGenes(genome)
#' WriteGenes(x[1:10,], format="gbk")
#' WriteGenes(x[1:10,], format="gff")
#' 
#' @export WriteGenes
WriteGenes <- function(x,
	file="",
	format="gbk",
	append=FALSE) {
	
	# error checking
	if (!is(x, "Genes"))
		stop("x must be an object of class 'Genes'.")
	FORMATS <- c("gbk", "gff")
	if (length(format)==0)
		stop("No format specified.")
	format <- pmatch(format, FORMATS)
	if (is.na(format))
		stop("Invalid format.")
	if (format==-1)
		stop("Ambiguous format.")
	if (!is.logical(append))
		stop("append must be a logical.")
	
	if (is.character(file)) {
		if (file == "") {
			file <- stdout()
		} else if (substring(file, 1L, 1L) == "|") {
			file <- pipe(substring(file, 2L), "w")
			on.exit(close(file))
		} else {
			file <- file(file, "w")
			on.exit(close(file))
		}
	}
	
	ns <- names(attr(x, "widths"))
	if (is.null(ns)) {
		ns <- sn <- seq_along(attr(x, "widths"))
	} else {
		sn <- strsplit(ns , " ", fixed=TRUE)
		sn <- sapply(sn, head, n=1)
	}
	
	w <- which(x[, "Gene"]==1)
	if (length(w)==0L)
		stop("No genes specified by x.")
	w <- w[order(x[w, "Index"], x[w, "Begin"])]
	x <- x[w,]
	
	if (format==1L) { # gbk
		if (!append) # overwrite the file
			cat("", file=file)
		t <- tapply(seq_len(nrow(x)),
			x[, "Index"],
			c)
		for (i in seq_along(t)) {
			cat("DEFINITION  ",
				ns[as.numeric(names(t)[i])],
				"\nFEATURES             Location/Qualifiers\n",
				sep="",
				file=file,
				append=TRUE)
			for (j in seq_along(t[[i]])) {
				if (x[t[[i]][j], "Strand"]==1L) {
					cat(ifelse(1/x[t[[i]][j], "Gene"] < 0,
							"     misc_RNA        complement(",
							"     CDS             complement("),
						x[t[[i]][j], "Begin"],
						"..",
						x[t[[i]][j], "End"],
						")",
						sep="",
						file=file,
						append=TRUE)
				} else {
					cat(ifelse(1/x[t[[i]][j], "Gene"] < 0,
							"     misc_RNA        ",
							"     CDS             "),
						x[t[[i]][j], "Begin"],
						"..",
						x[t[[i]][j], "End"],
						sep="",
						file=file,
						append=TRUE)
				}
				cat("\n                     /note=\"ID=",
					names(t)[i],
					"_",
					j,
					";\"\n",
					sep="",
					file=file,
					append=TRUE)
			}
		}
	} else { # gff
		cat("##gff-version 3\n",
			file=file,
			append=append)
		t <- tapply(seq_len(nrow(x)),
			x[, "Index"],
			c)
		for (i in seq_along(t)) {
			j <- as.numeric(names(t)[i])
			cat("##sequence-region ",
				sn[j],
				" 1 ",
				attr(x, "widths")[j],
				"\n",
				sep="",
				file=file,
				append=TRUE)
			tab <- cbind(sn[j],
				paste("DECIPHER_v",
					packageVersion("DECIPHER"),
					sep=""),
				ifelse(1/x[t[[i]], "Gene"] < 0,
					"transcript",
					"CDS"),
				x[t[[i]], "Begin"],
				x[t[[i]], "End"],
				round(100*x[t[[i]], "FractionReps"]),
				ifelse(x[t[[i]], "Strand"]==1L,
					"-",
					"+"),
				ifelse(1/x[t[[i]], "Gene"] < 0,
					".",
					"0"),
				paste("ID=",
					x[t[[i]], "Index"],
					"_",
					seq_along(t[[i]]),
					sep=""))
			write.table(tab,
				file=file,
				sep="\t",
				quote=FALSE,
				row.names=FALSE,
				col.names=FALSE,
				append=TRUE)
		}
	}
	
	invisible(NULL)
}
