#' Map Changes in Ancestral Character States
#' 
#' Maps character changes on a phylogenetic tree containing reconstructed
#' ancestral states.
#' 
#' Ancestral state reconstruction affords the ability to identify character
#' changes that occurred along edges of a rooted phylogenetic tree.  Character
#' changes are reported according to their index in \code{refPositions}.  If
#' \code{ignoreIndels} is \code{FALSE}, adjacent insertions and deletions are
#' merged into single changes occurring at their first position.  The table of
#' changes can be used to identify parallel, convergent, and divergent
#' mutations.
#' 
#' @name MapCharacters
#' @param x An object of class \code{dendrogram} with \code{"state"} attributes
#' for each node.
#' @param refPositions Numeric vector of reference positions in the original
#' sequence alignment.  Only changes at \code{refPositions} are reported, and
#' state changes are labeled according to their position in
#' \code{refPositions}.
#' @param labelEdges Logical determining whether to label edges with the number
#' of changes along each edge.
#' @param type Character string indicating the type of output desired.  This
#' should be (an abbreviation of) one of \code{"dendrogram"}, \code{"table"},
#' or \code{"both"}.  (See value section below.)
#' @param chars Character vector specifying the characters to consider in state
#' changes at each site.  The default (\code{LETTERS}) is to consider any upper
#' case letter. Alternatively, \code{chars} could be \code{AA_STANDARD},
#' \code{DNA_BASES}, or \code{RNA_BASES}.
#' @param ignoreIndels Logical specifying whether to report insertions and
#' deletions (indels).  If \code{TRUE} (the default), only substitutions of one
#' state with another are reported.
#' @return If \code{type} is \code{"dendrogram"} (the default) then the
#' original \code{dendrogram} \code{x} is returned with the addition of
#' \code{"change"} attributes on every edge except the root.  If \code{type} is
#' \code{"table"} then a sorted \code{table} of character changes is returned
#' with the most frequent parallel changes at the beginning.  If \code{type} is
#' \code{"both"} then a \code{list} of length 2 is provided containing both the
#' \code{dendrogram} and \code{table}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{TreeLine}}
#' @examples
#' 
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # align the sequences
#' rna <- RNAStringSet(RemoveGaps(dna))
#' rna <- AlignSeqs(rna)
#' rna # input alignment
#' 
#' d <- DistanceMatrix(rna, type="dist", correction="JC")
#' tree <- TreeLine(myDistMatrix=d,
#'                  method="NJ",
#'                  type="dendrogram",
#'                  myXStringSet=rna,
#'                  reconstruct=TRUE)
#' 
#' out <- MapCharacters(tree,
#'                      labelEdges=TRUE,
#'                      type="both",
#'                      chars=RNA_BASES)
#' 
#' # plot the tree with defaults
#' tree <- out[[1]]
#' plot(tree, horiz=TRUE) # edges show number of changes
#' 
#' # color edges by number of changes
#' maxC <- 200 # changes at maximum of color spectrum
#' colors <- colorRampPalette(c("black", "darkgreen", "green"))(maxC)
#' colorEdges <- function(x) {
#'    num <- attr(x, "edgetext") + 1
#'    if (length(num)==0)
#'        return(x)
#'    if (num > maxC)
#'        num <- maxC
#'    attr(x, "edgePar") <- list(col=colors[num])
#'    attr(x, "edgetext") <- NULL
#'    return(x)
#' }
#' colorfulTree <- dendrapply(tree, colorEdges)
#' plot(colorfulTree, horiz=TRUE, leaflab="none")
#' 
#' # look at parallel changes (X->Y)
#' parallel <- out[[2]]
#' head(parallel) # parallel changes
#' 
#' # look at convergent changes (*->Y)
#' convergent <- gsub(".*?([0-9]+.*)", "\\1", names(parallel))
#' convergent <- tapply(parallel, convergent, sum)
#' convergent <- sort(convergent, decreasing=TRUE)
#' head(convergent)
#' 
#' # look at divergent changes (X->*)
#' divergent <- gsub("(.*[0-9]+).*", "\\1", names(parallel))
#' divergent <- tapply(parallel, divergent, sum)
#' divergent <- sort(divergent, decreasing=TRUE)
#' head(divergent)
#' 
#' # plot number of changes by position
#' changes <- gsub(".*?([0-9]+).*", "\\1", names(parallel))
#' changes <- tapply(parallel, changes, sum)
#' plot(as.numeric(names(changes)),
#'      changes,
#'      xlab="Position",
#'      ylab="Total independent changes")
#' 
#' # count cases of potential compensatory mutations
#' compensatory <- dendrapply(tree,
#'     function(x) {
#'         change <- attr(x, "change")
#'         pos <- as.numeric(gsub(".*?([0-9]+).*", "\\1", change))
#'         e <- expand.grid(seq_along(pos), seq_along(pos))
#'         e <- e[pos[e[, 1]] < pos[e[, 2]],]
#'         list(paste(change[e[, 1]], change[e[, 2]], sep=" & "))
#'     })
#' compensatory <- unlist(compensatory)
#' u <- unique(compensatory)
#' m <- match(compensatory, u)
#' m <- tabulate(m, length(u))
#' compensatory <- sort(setNames(m, u), decreasing=TRUE)
#' head(compensatory) # ranked list of concurrent mutations
#' 
#' @export MapCharacters
MapCharacters <- function(x,
	refPositions=seq_len(nchar(attr(x, "state")[1])),
	labelEdges=FALSE,
	type="dendrogram",
	chars=LETTERS,
	ignoreIndels=TRUE) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x is not a dendrogram.")
	if (!is.logical(labelEdges))
		stop("labelEdges must be a logical.")
	if (!is.numeric(refPositions))
		stop("refPositions must be a numeric.")
	if (any(refPositions > nchar(attr(x, "state")[1])) ||
		any(refPositions < 1))
		stop("refPositions out of bounds.")
	if (is.unsorted(refPositions))
		stop("refPositions must be in ascending order.")
	if (length(chars) <= 1)
		stop("At least two chars must be specified.")
	if (any(is.na(chars)))
		stop("Chars cannot contain NA values.")
	if (any(nchar(chars) != 1))
		stop("All chars must be a single character.")
	gaps <- "-"
	if (any(chars %in% gaps))
		stop("Gap ('-') characters cannot be in chars.")
	if (!is.logical(ignoreIndels))
		stop("ignoreIndels must be a logical.")
	TYPES <- c("dendrogram", "table", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	
	.map <- function(x) {
		# initialize a stack of maximum length (l)
		stack <- vector("list", l)
		visit <- logical(l) # node already visited
		parent <- integer(l) # index of parent node
		index <- integer(l) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- x
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit
				
				for (i in seq_along(stack[[pos]])) {
					s1 <- attr(stack[[pos]], "state") # root ancestral state
					if (is.null(s1))
						stop("Missing 'state' attribute in x.")
					s1 <- s1[1]
					if (!is.character(s1))
						stop("Incorrect 'state' attribute in x.")
					s1 <- strsplit(s1, "", fixed=TRUE)[[1]]
					
					s2 <- attr(stack[[pos]][[i]], "state")
					if (is.null(s2))
						stop("Missing 'state' attribute in x.")
					s2 <- s2[1]
					if (!is.character(s2))
						stop("Incorrect 'state' attribute in x.")
					s2 <- strsplit(s2, "", fixed=TRUE)[[1]]
					if (length(s1) != length(s2))
						stop("Inconsistent 'state' attributes in x.")
					
					if (ignoreIndels) {
						nongaps <- which(!(s1 %in% gaps) &
							!(s2 %in% gaps))
					} else {
						nongaps <- which(!(s1 %in% gaps) |
							!(s2 %in% gaps))
					}
					w <- which(s1[nongaps] != s2[nongaps] &
						s1[nongaps] %in% chars &
						s2[nongaps] %in% chars)
					
					if (length(w) > 0) {
						ini <- s1[nongaps[w]]
						fin <- s2[nongaps[w]]
						
						ins <- which(s1[nongaps[w]]=="-")
						if (length(ins) > 0) {
							counts <- rep(1L, length(ins))
							remove <- logical(length(ins))
							start <- 1L
							for (j in seq_along(ins)[-1]) {
								if (w[ins[j]]==w[ins[j - 1]] + 1L) {
									counts[start] <- counts[start] + 1L
									remove[j] <- TRUE
								} else {
									start <- j
								}
							}
							ini[ins] <- "ins("
							fin[ins] <- paste(",",
								sapply(mapply(seq,
									ins,
									ins + counts - 1L,
									SIMPLIFY=FALSE),
									function(s) {
										paste(s2[nongaps[w[s]]],
										collapse="")
									}),
								")",
								sep="")
							remove <- which(remove)
							if (length(remove) > 0) {
								w <- w[-ins[remove]]
								ini <- ini[-ins[remove]]
								fin <- fin[-ins[remove]]
							}
						}
						
						del <- which(s2[nongaps[w]]=="-")
						if (length(del) > 0) {
							counts <- rep(1L, length(del))
							remove <- logical(length(del))
							start <- 1L
							for (j in seq_along(del)[-1]) {
								if (w[del[j]]==w[del[j - 1]] + 1L) {
									counts[start] <- counts[start] + 1L
									remove[j] <- TRUE
								} else {
									start <- j
								}
							}
							ini[del] <- "del("
							fin[del] <- paste(",",
								sapply(mapply(seq,
									del,
									del + counts - 1L,
									SIMPLIFY=FALSE),
									function(s) {
										paste(s1[nongaps[w[s]]],
										collapse="")
									}),
								")",
								sep="")
							remove <- which(remove)
							if (length(remove) > 0) {
								w <- w[-del[remove]]
								ini <- ini[-del[remove]]
								fin <- fin[-del[remove]]
							}
						}
						
						w <- match(nongaps[w], refPositions)
						keep <- which(!is.na(w))
						if (length(keep) > 0) {
							attr(stack[[pos]][[i]], "change") <- paste(ini[keep],
								w[keep],
								fin[keep],
								sep="")
							if (type > 1L) {
								res[[num]] <<- attr(stack[[pos]][[i]], "change")
								num <<- num + 1L
							}
						} else {
							attr(stack[[pos]][[i]], "change") <- character()
						}
						if (labelEdges)
							attr(stack[[pos]][[i]], "edgetext") <- length(keep)
					} else {
						attr(stack[[pos]][[i]], "change") <- character()
						if (labelEdges)
							attr(stack[[pos]][[i]], "edgetext") <- 0L
					}
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
	
	l <- attr(x, "members")
	if (type > 1L) {
		res <- vector("list",
			2*l - 2)
		num <- 1L
	}
	out <- .map(x)
	attr(out, "change") <- character() # empty set at root node
	if(type > 1L) {
		res <- table(unlist(res))
		res <- sort(res, decreasing=TRUE)
		if (type==2L) { # table
			out <- res
		} else if (type==3L) { # both
			out <- list(out, res)
		}
	}
	
	invisible(out)
}
