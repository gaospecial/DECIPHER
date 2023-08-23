#' Compute cophenetic distances on dendrogram objects
#' 
#' Calculates the matrix of cophenetic distances represented by a dendrogram
#' object.
#' 
#' The cophenetic distance between two observations is defined as the branch
#' length separating them on a dendrogram.  This function differs from the
#' \code{cophenetic} function in that it does not assume the tree is
#' ultrametric and outputs the branch length separating pairs of observations
#' rather than the height of their merger. A dendrogram that better preserves a
#' distance matrix will show higher correlation between the distance matrix and
#' it cophenetic distances.
#' 
#' @name Cophenetic
#' @param x A dendrogram object.
#' @return An object of class 'dist'.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{TreeLine}}
#' @examples
#' 
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' d1 <- DistanceMatrix(dna, type="dist")
#' dend <- TreeLine(myDistMatrix=d1, method="NJ")
#' d2 <- Cophenetic(dend)
#' cor(d1, d2)
#' 
#' @export Cophenetic
Cophenetic <- function(x) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x must be an object of class 'dendrogram'.")
	
	n <- attr(x, "members")
	d <- numeric(n*(n - 1)/2)
	u <- unlist(x)
	o <- order(u)
	u <- u[o]
	labs <- rapply(x,
		function(x)
			attr(x, "label"))
	labs <- labs[o]
	x <- rapply(x,
		function(y) {
			y[] <- match(y[1L], u)
			y
		},
		how="replace")
	
	.dist <- function(dend) {
		# initialize a stack of maximum length (n)
		stack <- vector("list", n)
		visit <- logical(n) # node already visited
		parent <- integer(n) # index of parent node
		index <- integer(n) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- dend
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit
				
				for (k in seq_along(stack[[pos]])) {
					h <- attr(stack[[pos]], "height") - attr(stack[[pos]][[k]], "height")
					I <- unlist(stack[[pos]][[k]])
					d <- .Call("cophenetic", # in-place change of d (requires previous temporary copy)
						I,
						n,
						d,
						h,
						PACKAGE="DECIPHER")
				}
				
				# replace self in parent
				if (parent[pos] > 0)
					stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
				pos <- pos - 1L # pop off of stack
			} else { # descending tree
				visit[pos] <- TRUE
				p <- pos
				for (i in seq_along(stack[[p]])) {
					if (!is.leaf(stack[[p]][[i]])) {
						# push subtree onto stack
						pos <- pos + 1L
						stack[[pos]] <- stack[[p]][[i]]
						parent[[pos]] <- p
						index[[pos]] <- i
					}
				}
			}
		}
		return(stack[[1L]])
	}
	.dist(x)
	class(d) <- "dist"
	attr(d, "Size") <- n
	attr(d, "Diag") <- TRUE
	attr(d, "Upper") <- TRUE
	attr(d, "Labels") <- labs
	
	return(d)
}
