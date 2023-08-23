#' Read a Dendrogram from a Newick Formatted File
#' 
#' Reads a dendrogram object from a file in Newick (also known as New
#' Hampshire) parenthetic format.
#' 
#' \code{ReadDendrogram} will create a dendrogram object from a Newick
#' formatted tree.  Note that all edge lengths must be specified, but labels
#' are optional.  Leaves will be numbered by their labels in alphabetical
#' order.
#' 
#' @name ReadDendrogram
#' @param file a connection object or a character string.
#' @param convertBlanks Logical specifying whether to convert underscores in
#' unquoted leaf labels to spaces.
#' @param internalLabels Logical indicating whether to keep internal node
#' labels as ``edgetext'' preceding the node in the \code{dendrogram}.
#' @param keepRoot Logical specifying whether to keep the root node (if one is
#' present) as a dendrogram leaf.
#' @return An object of class \code{dendrogram}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{TreeLine}}, \code{\link{WriteDendrogram}}
#' @examples
#' 
#' tf <- tempfile()
#' dists <- matrix(c(0, 10, 20, 10, 0, 5, 20, 5, 0),
#'     nrow=3,
#'     dimnames=list(c("dog", "elephant", "horse")))
#' dend1 <- TreeLine(myDistMatrix=dists, method="NJ", type="dendrogram")
#' WriteDendrogram(dend1, file=tf)
#' 
#' dend2 <- ReadDendrogram(tf)
#' layout(matrix(1:2))
#' plot(dend1, main="Dendrogram Written")
#' plot(dend2, main="Dendrogram Read")
#' 
#' unlink(tf)
#' 
#' @export ReadDendrogram
ReadDendrogram <- function(file,
	convertBlanks=TRUE,
	internalLabels=TRUE,
	keepRoot=TRUE) {
	
	# error checking
	if (!is.logical(convertBlanks))
		stop("convertBlanks must be a logical.")
	if (!is.logical(internalLabels))
		stop("internalLabels must be a logical.")
	if (!is.logical(keepRoot))
		stop("keepRoot must be a logical.")
	
	r <- readLines(file, warn=FALSE)
	w <- which(nchar(r) > 0)
	if (length(w) > 1) {
		r <- paste(r[w], collapse="")
	} else if (length(w)==0) {
		stop("file is empty.")
	} else {
		r <- r[w]
	}
	
	r <- strsplit(r,
		'(?=[\\[\\](),:;])(?=([^"]*"[^"]*")*[^"]*$)',
		perl=TRUE)[[1]]
	r <- gsub("^\\s+|\\s+$", "", r)
	w <- which(r=="")
	if (length(w) > 0)
		r <- r[-w]
	
	getLab <- function(LAB) {
		# convert underscores to spaces in unquoted labels
		lab <- gsub("^\"(.*)\"$", "\\1", LAB)
		if (convertBlanks && nchar(lab)==nchar(LAB))
			lab <- gsub("_", " ", lab, fixed=TRUE)
		if (nchar(lab)!=nchar(LAB))
			lab <- gsub("''", "'", lab, fixed=TRUE)
		return(lab)
	}
	
	warned <- FALSE
	i <- k <- 1L
	.newick2list <- function() {
		x <- vector("list")
		j <- 0L
		while (i <= length(r) && r[i] != ";") {
			if (r[i]=="[") { # comment
				count <- 1L
				i <<- i + 1L
				while (count > 0) {
					if (i > length(r))
						stop("Improperly formatted comment.")
					if (r[i]=="]") {
						count <- count - 1L
					} else if (r[i]=="[") {
						count <- count + 1L
					}
					i <<- i + 1L
				}
			} else if (r[i]==")") {
				i <<- i + 1L
				if (i <= length(r) && r[i]==":") {
					# internal node
					i <<- i + 1L
					if (i <= length(r))
						attr(x, "height") <- as.numeric(r[i])
					i <<- i + 1L
				} else if (i < length(r) && r[i + 1L]==":") {
					# labeled internal node
					if (internalLabels)
						attr(x, "edgetext") <- getLab(r[i])
					i <<- i + 2L
					if (i <= length(r))
						attr(x, "height") <- as.numeric(r[i])
					i <<- i + 1L
				} else if (i < length(r) && r[i + 1L]==";") {
					i <<- i + 1L
				} else if (i <= length(r) && r[i] != ";") {
					stop("Unsupported file formatting.")
				}
				break
			} else if (r[i]=="(") {
				j <- j + 1L
				i <<- i + 1L
				x[[j]] <- .newick2list()
			} else if (r[i]==",") {
				i <<- i + 1L
			} else if ((i + 2L) <= length(r) && r[i + 1L] == ":") {
				j <- j + 1L
				x[[j]] <- k
				k <<- k + 1L
				attr(x[[j]], "leaf") <- TRUE
				
				attr(x[[j]], "label") <- getLab(r[i])
				attr(x[[j]], "height") <- as.numeric(r[i + 2L])
				attr(x[[j]], "members") <- 1L
				i <<- i + 3L
			} else if ((i + 1L) <= length(r) && r[i] == ":") {
				j <- j + 1L
				x[[j]] <- k
				k <<- k + 1L
				attr(x[[j]], "leaf") <- TRUE
				
				attr(x[[j]], "label") <- ""
				attr(x[[j]], "height") <- as.numeric(r[i + 1L])
				attr(x[[j]], "members") <- 1L
				i <<- i + 2L
			} else if (r[i]==" ") {
				i <<- i + 1L
			} else {
				stop("Unsupported file formatting.")
			}
		}
		
		return(x)
	}
	
	.list2dendrogram <- function(x) {
		memb <- 0L
		tot <- 0
		for (i in seq_along(x)) {
			if (is.null(attr(x[[i]], "leaf"))) {
				x[[i]] <- .list2dendrogram(x[[i]])
				tot <- tot + memb + attr(x[[i]], "midpoint") + 1
				memb <- memb + attr(x[[i]], "members")
			} else {
				memb <- memb + 1L
				tot <- tot + memb
			}
		}
		attr(x, "members") <- memb
		attr(x, "midpoint") <- tot/length(x) - 1
		
		return(x)
	}
	
	.sumHeights <- function(x) {
		maxH <- 0
		for (i in seq_along(x)) {
			if (is.null(attr(x[[i]], "leaf"))) {
				x[[i]] <- .sumHeights(x[[i]])
				newH <- attr(x[[i]], "cum")
			} else {
				newH <- attr(x[[i]], "height")
			}
			if (newH > maxH)
				maxH <- newH
		}
		
		if (!is.null(attr(x, "height"))) {
			attr(x, "cum") <- maxH + attr(x, "height")
		} else {
			attr(x, "height") <- maxH
		}
		
		return(x)
	}
	
	.setHeights <- function(x, current) {
		for (i in seq_along(x)) {
			new <- current - attr(x[[i]], "height")
			if (is.null(attr(x[[i]], "leaf"))) {
				x[[i]] <- .setHeights(x[[i]],
					current=new)
			} else {
				attr(x[[i]], "height") <- current - attr(x[[i]], "height")
			}
			attr(x[[i]], "height") <- new
		}
		
		return(x)
	}
	
	x <- .newick2list()[[1]]
	if (!is.null(attr(x, "height"))) {
		if (keepRoot) {
			height <- attr(x, "height")
		} else {
			height <- NULL
		}
		attr(x, "height") <- NULL
	} else {
		height <- NULL
	}
	x <- .list2dendrogram(x)
	x <- .sumHeights(x)
	x <- .setHeights(x, current=attr(x, "height"))
	if (!is.null(height)) {
		x <- list(x)
		attr(x, "members") <- attr(x[[1]], "members")
		attr(x, "midpoint") <- attr(x[[1]], "midpoint")
		attr(x, "height") <- attr(x[[1]], "height") + height
		attr(x, "label") <- NULL
	}
	class(x) <- "dendrogram"
	x <- dendrapply(x, function(x) {
			attr(x, "cum") <- NULL
			x
		})
	
	# reorder numbers by label
	labs <- rapply(x,
		function(x)
			attr(x, "label"))
	if (any(duplicated(labs))) {
		warning("Leaf numbering is unordered because of duplicated leaf labels.")
	} else {
		o <- order(order(labs))
		x <- rapply(x,
			function(x) {
				x[] <- o[x]
				x
			},
			how="replace")
	}
	
	return(x)
}
