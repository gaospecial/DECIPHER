/****************************************************************************
 *                       Remove Gaps in XStringSets                         *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// strcpy
#include <string.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("removeCommonGaps", sequences, type, processors, PACKAGE="DECIPHER")
SEXP removeCommonGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads)
{
	int i, j, k, l, w, x_length, *width, sum, start, delta;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int m = asInteger(mask);
	int nthreads = asInteger(nThreads);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	Chars_holder x_i;
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	l = x_i.length;
	
	// initialize a vector of columns that are 100% gaps
	int *remove = (int *) R_alloc(l, sizeof(int));
	sum = 0; // number of columns to remove
	if (t==3) { // AAStringSet
		for (i = 0; i < l; i++) {
			if (!(x_i.ptr[i] ^ 0x2D) || !(x_i.ptr[i] ^ 0x2E) || (m && !(x_i.ptr[i] ^ 0x2B))) { // position is a gap
				remove[sum] = i;
				sum++;
			}
		}
		
		// find columns to remove
		for (i = 1; i < x_length; i++) { // each sequence
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			if (x_i.length!=l)
				error("Sequences are not equal length.");
			
			for (j = 0; j < sum; j++) {
				if (x_i.ptr[remove[j]] ^ 0x2D && x_i.ptr[remove[j]] ^ 0x2E && (!m || x_i.ptr[remove[j]] ^ 0x2B)) { // non-gap in this position
					// stop the position from being removed
					sum--;
					// shift all positions over by one
					for (k = j; k < sum; k++)
						remove[k] = remove[k + 1];
					j--;
				}
			}
		}
	} else { // DNAStringSet or RNAStringSet
		for (i = 0; i < l; i++) {
			if (x_i.ptr[i] & 0x10 || x_i.ptr[i] & 0x40 || (m && x_i.ptr[i] & 0x20)) { // position is a gap
				remove[sum] = i;
				sum++;
			}
		}
		
		// find columns to remove
		for (i = 1; i < x_length; i++) { // each sequence
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			if (x_i.length!=l)
				error("Sequences are not equal length.");
			
			for (j = 0; j < sum; j++) {
				if (!(x_i.ptr[remove[j]] & 0x10 || x_i.ptr[remove[j]] & 0x40 || (m && x_i.ptr[remove[j]] & 0x20))) { // non-gap in this position
					// stop the position from being removed
					sum--;
					// shift all positions over by one
					for (k = j; k < sum; k++)
						remove[k] = remove[k + 1];
					j--;
				}
			}
		}
	}
	
	// determine the widths of the aligned (equal width) XStringSet
	PROTECT(ans_width = NEW_INTEGER(x_length));
	w = l - sum;
	for (i = 0, width = INTEGER(ans_width); i < x_length; i++, width++)
		*width = w;
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t==1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t==2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t==3
		strcpy(ans_classname, "AAStringSet");
	}
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	Chars_holder ans_elt_holder;
	
	#pragma omp parallel for private(i,j,ans_elt_holder,x_i,start,delta) schedule(guided) num_threads(nthreads)
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		start = 0; // position in ans_elt_holder.ptr
		for (j = 0; j < sum; j++) {
			if (j > 0) {
				delta = remove[j] - remove[j - 1] - 1;
				if (delta > 0) {
					memcpy((char *) ans_elt_holder.ptr + start, x_i.ptr + remove[j - 1] + 1, delta * sizeof(char));
					start += delta;
				}
			} else {
				memcpy((char *) ans_elt_holder.ptr, x_i.ptr, remove[j] * sizeof(char));
				start = remove[j];
			}
		}
		if (start < w) {
			delta = w - start;
			memcpy((char *) ans_elt_holder.ptr + start, x_i.ptr + l - delta, delta * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}

//ans_start <- .Call("removeGaps", sequences, type, processors, PACKAGE="DECIPHER")
SEXP removeGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads)
{
	int i, j, x_length, p, sum;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int m = asInteger(mask);
	int nthreads = asInteger(nThreads);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	Chars_holder x_i, ans_elt_holder;
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// count the sequence lengths
	PROTECT(ans_width = NEW_INTEGER(x_length));
	int *width = INTEGER(ans_width);
	#pragma omp parallel for private(i,j,x_i) schedule(guided) num_threads(nthreads)
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		width[i] = x_i.length;
		if (t==3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E) || (m && !(x_i.ptr[j] ^ 0x2B))) { // position is a gap
					width[i]--;
				}
			}
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40 || (m && x_i.ptr[j] & 0x20)) { // position is a gap
					width[i]--;
				}
			}
		}
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t==1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t==2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t==3
		strcpy(ans_classname, "AAStringSet");
	}
	
	// initialize a new XStringSet
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	
	#pragma omp parallel for private(i,j,ans_elt_holder,x_i,p,sum) schedule(guided) num_threads(nthreads)
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		sum = 0; // number of positions to copy
		p = 0; // position in ans_elt_holder
		if (t==3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E) || (m && !(x_i.ptr[j] ^ 0x2B))) { // position is a gap
					if (sum > 0) {
						memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
						p += sum;
						sum = 0;
					}
				} else {
					sum++;
				}
			}
			if (sum > 0)
				memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40 || (m && x_i.ptr[j] & 0x20)) { // position is a gap
					if (sum > 0) {
						memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
						p += sum;
						sum = 0;
					}
				} else {
					sum++;
				}
			}
			if (sum > 0)
				memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}
