/****************************************************************************
 *                      Utilities for Gene Prediction                       *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

int getBase(const char p)
{
	switch (p) {
		case 1: // A
			return 0;
		case 2: // C
			return 1;
		case 4: // G
			return 2;
		case 8: // T
			return 3;
		default: // other
			return 100000;
	}
}

int getBaseRC(const char p)
{
	switch (p) {
		case 1: // A
			return 3; // T
		case 2: // C
			return 2; // G
		case 4: // G
			return 1; // C
		case 8: // T
			return 0; // A
		default: // other
			return 100000;
	}
}

char getBaseLetter(const char p)
{
	switch (p) {
		case 1: // A
			return 'A';
		case 2: // C
			return 'C';
		case 4: // G
			return 'G';
		case 8: // T
			return 'T';
		default: // other
			return '\0';
	}
}

char getBaseLetterRC(const char p)
{
	switch (p) {
		case 1: // A
			return 'T'; // T
		case 2: // C
			return 'G'; // G
		case 4: // G
			return 'C'; // C
		case 8: // T
			return 'A'; // A
		default: // other
			return '\0';
	}
}

int nextCount(int count, int tot, int *orfs, int minL)
{
	int c = count;
	while ((c < tot &&
		(orfs[c + 3*tot] - orfs[c + 2*tot] + 1) < minL) ||
		((c + 1) < tot &&
		((orfs[c + 3*tot]==orfs[c + 3*tot + 1] &&
		orfs[c + tot]==0) ||
		(orfs[c + 2*tot]==orfs[c + 2*tot + 1] &&
		orfs[c + tot]==1)) &&
		orfs[c]==orfs[c + 1]))
		c++;
	if (c==tot)
		c--;
	
	return c;
}

SEXP getORFs(SEXP x, SEXP start_codons, SEXP stop_codons, SEXP min_gene_length, SEXP allow_edges)
{
	int i, j, k, s, rf, val, x_length, lastStop;
	int lstarts = length(start_codons);
	int lstops = length(stop_codons);
	int *starts = INTEGER(start_codons);
	int *stops = INTEGER(stop_codons);
	int minL = asInteger(min_gene_length);
	int allow = asInteger(allow_edges);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	int count = 0;
	int size = 10000;
	int *ORFstarts = Calloc(size, int);
	int *ORFstops = Calloc(size, int);
	int *ORFstrands = Calloc(size, int);
	int *ORFindices = Calloc(size, int);
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (s = 0; s <= 1; s++) { // strand
			for (rf = 1; rf <= 3; rf++) { // reading frame
				if (allow) {
					lastStop = x_i.length - rf - 3;
				} else {
					lastStop = -1;
				}
				for (j = x_i.length - rf; j >= 2;) { // position
					if (s) { // negative strand
						val = getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 4*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 16*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
					} else { // positive strand
						val = getBase(x_i.ptr[j--]);
						val += 4*getBase(x_i.ptr[j--]);
						val += 16*getBase(x_i.ptr[j--]);
					}
					
					for (k = 0; k < lstops; k++) {
						if (val==stops[k]) {
							lastStop = j;
							break;
						}
					}
					if (lastStop==j) // current position is a stop
						continue;
					
					for (k = 0; k < lstarts; k++) {
						if ((val==starts[k] || // current position is a start
							(allow && j <= 1)) &&
							(lastStop - j + 3) >= minL) {
							if (count >= size) {
								size += 10000;
								ORFstarts = Realloc(ORFstarts, size, int);
								ORFstops = Realloc(ORFstops, size, int);
								ORFstrands = Realloc(ORFstrands, size, int);
								ORFindices = Realloc(ORFindices, size, int);
							}
							if (s) { // negative strand
								ORFstarts[count] = x_i.length - lastStop - 3;
								ORFstops[count] = x_i.length - j - 1;
							} else { // positive strand
								ORFstarts[count] = j + 2;
								ORFstops[count] = lastStop + 4;
							}
							ORFstrands[count] = s;
							ORFindices[count] = i + 1;
							count++;
							break;
						}
					}
				}
			}
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, count, 4));
	int *rans = INTEGER(ans);
	for (i = 0; i < count; i++) {
		rans[i] = ORFindices[i];
		rans[i + count] = ORFstrands[i];
		rans[i + 2*count] = ORFstarts[i];
		rans[i + 3*count] = ORFstops[i];
	}
	
	Free(ORFstarts);
	Free(ORFstops);
	Free(ORFstrands);
	Free(ORFindices);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP codonModel(SEXP x, SEXP orftable, SEXP stop_codons, SEXP min_orf_length)
{
	int i, j, k, rf, s, val, x_length;
	int lstops = length(stop_codons);
	int *stops = INTEGER(stop_codons);
	int minL = asInteger(min_orf_length);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int inside = 0;
	
	int *freq = Calloc(64, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// start at first ORF that is at least minL in length
	int count = nextCount(0, tot, orfs, minL);
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (s = 0; s <= 1; s++) { // strand
			for (rf = 1; rf <= 3; rf++) { // reading frame
				for (j = x_i.length - rf; j >= 2;) { // position
					if (s) { // negative strand
						val = getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 4*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 16*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
					} else { // positive strand
						val = getBase(x_i.ptr[j--]);
						val += 4*getBase(x_i.ptr[j--]);
						val += 16*getBase(x_i.ptr[j--]);
					}
					
					if (inside) { // look for beg
						if (s) { // negative strand
							if ((x_i.length - j - 1)==orfs[count + 3*tot]) {
								inside = 0;
								count = nextCount(++count, tot, orfs, minL);
							} else if (val < 64) {
								freq[val]++;
							}
						} else { // positive strand
							if ((j + 2)==orfs[count + 2*tot]) {
								inside = 0;
								count = nextCount(++count, tot, orfs, minL);
							} else if (val < 64) {
								freq[val]++;
							}
						}
					} else { // look for end
						if (s) { // negative strand
							if ((x_i.length - j - 3)==orfs[count + 2*tot]) {
								inside = 1;
							} else if (val < 64) {
								bg[val]++;
							}
						} else { // positive strand
							if ((j + 4)==orfs[count + 3*tot]) {
								inside = 1;
							} else if (val < 64) {
								bg[val]++;
							}
						}
					}
				}
			}
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	int notstop;
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			sumfreq += freq[i];
			sumbg += bg[i];
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			rans[i] = log(((double)freq[i]/(double)sumfreq)/((double)bg[i]/(double)sumbg));
		} else {
			rans[i] = 0;
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreCodonModel(SEXP x, SEXP orftable, SEXP codon_scores)
{
	int i, j, s, fin, val, lastVal, dicodon;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *codons = REAL(codon_scores);
	if (length(codon_scores)==64) {
		dicodon = 0;
	} else if (length(codon_scores)==4096) {
		dicodon = 1;
	} else {
		error("codon_scores is the wrong length.");
	}
	int curr_i = 0;
	double score;
	//int first_i, max_i;
	//double max_score;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot;) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start and stop positions
		s = orfs[i + tot];
		if (s) { // negative strand
			fin = orfs[i + 3*tot] - 3; // finish at start codon - 1 codon
			j = orfs[i + 2*tot] + 2; // start at stop codon + 1 codon
		} else { // positive strand
			fin = orfs[i + 2*tot] + 1; // finish at start codon + 1 codon
			j = orfs[i + 3*tot] - 4; // start at stop codon - 1 codon
		}
		
		// sum score from j to fin
		score = 0;
		if (dicodon)
			lastVal = 100000;
		//first_i = i;
		//max_score = -1e52;
		while (1) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (dicodon) {
				if (val < 64 && lastVal < 64)
					score += codons[lastVal*64 + val];
				lastVal = val;
			} else {
				if (val < 64)
					score += codons[val];
			}
			
			if (j==fin) {
				// record the score
				rans[i] = score;
				//if (score > max_score) {
				//	max_score = score;
				//	max_i = i;
				//}
				i++;
				
				// check whether still within bounds
				if (s) { // negative strand
					if (i != tot &&
						orfs[i + 2*tot]==orfs[i + 2*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 3*tot] - 3;
					} else {
						break;
					}
				} else { // positive strand
					if (i != tot &&
						orfs[i + 3*tot]==orfs[i + 3*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 2*tot] + 1;
					} else {
						break;
					}
				}
			}
		}
		
		// subtract difference from top overlapping ORF
		//for (j = max_i - 1; j >= first_i; j--)
		//	rans[j] -= rans[max_i] - rans[j];
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP dicodonModel(SEXP x, SEXP orftable, SEXP stop_codons)
{
	int i, j, k, rf, s, val, lastVal, x_length;
	int lstops = length(stop_codons);
	int *stops = INTEGER(stop_codons);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int inside = 0;
	
	int *freq = Calloc(4096, int);
	int *bg = Calloc(4096, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	int count = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (s = 0; s <= 1; s++) { // strand
			for (rf = 1; rf <= 3; rf++) { // reading frame
				lastVal = 100000;
				for (j = x_i.length - rf; j >= 2;) { // position
					if (s) { // negative strand
						val = getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 4*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 16*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
					} else { // positive strand
						val = getBase(x_i.ptr[j--]);
						val += 4*getBase(x_i.ptr[j--]);
						val += 16*getBase(x_i.ptr[j--]);
					}
					
					if (inside) { // look for beg
						if (s) { // negative strand
							if ((x_i.length - j - 1)==(orfs[count + 3*tot])) {
								inside = 0;
								lastVal = 100000;
								count++;
							} else if (val < 64 && lastVal < 64) {
								freq[lastVal*64 + val]++;
							}
						} else { // positive strand
							if ((j + 2)==(orfs[count + 2*tot])) {
								inside = 0;
								lastVal = 100000;
								count++;
							} else if (val < 64 && lastVal < 64) {
								freq[lastVal*64 + val]++;
							}
						}
					} else { // look for end
						if (s) { // negative strand
							if ((x_i.length - j - 3)==orfs[count + 2*tot]) {
								inside = 1;
								lastVal = 100000;
							} else if (val < 64 && lastVal < 64) {
								bg[lastVal*64 + val]++;
							}
						} else { // positive strand
							if ((j + 4)==orfs[count + 3*tot]) {
								inside = 1;
								lastVal = 100000;
							} else if (val < 64 && lastVal < 64) {
								bg[lastVal*64 + val]++;
							}
						}
					}
					
					lastVal = val;
				}
			}
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	int notstop;
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			for (j = 0; j < 64; j++) {
				notstop = 1;
				for (k = 0; k < lstops; k++) {
					if (j==stops[k]) {
						notstop = 0;
						break;
					}
				}
				if (notstop) {
					if (freq[i*64 + j]==0) // add pseudocount
						freq[i*64 + j] = 1;
					if (bg[i*64 + j]==0) // add pseudocount
						bg[i*64 + j] = 1;
					sumfreq += freq[i*64 + j];
					sumbg += bg[i*64 + j];
				}
			}
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, 64, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			for (j = 0; j < 64; j++) {
				notstop = 1;
				for (k = 0; k < lstops; k++) {
					if (j==stops[k]) {
						notstop = 0;
						break;
					}
				}
				if (notstop) {
					rans[i*64 + j] = log(((double)freq[i*64 + j]/(double)sumfreq)/((double)bg[i*64 + j]/(double)sumbg))/2;
				} else {
					rans[i*64 + j] = 0;
				}
			}
		} else {
			for (j = 0; j < 64; j++)
				rans[i*64 + j] = 0;
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP unicodonModel(SEXP x, SEXP orftable, SEXP stop_codons)
{
	int i, j, k, rf, s, val, x_length;
	int lstops = length(stop_codons);
	int *stops = INTEGER(stop_codons);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int inside = 0;
	
	int *freq = Calloc(64, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	int count = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (s = 0; s <= 1; s++) { // strand
			for (rf = 1; rf <= 3; rf++) { // reading frame
				for (j = x_i.length - rf; j >= 2;) { // position
					if (s) { // negative strand
						val = getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 4*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 16*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
					} else { // positive strand
						val = getBase(x_i.ptr[j--]);
						val += 4*getBase(x_i.ptr[j--]);
						val += 16*getBase(x_i.ptr[j--]);
					}
					
					if (inside) { // look for beg
						if (s) { // negative strand
							if ((x_i.length - j - 1)==(orfs[count + 3*tot])) {
								inside = 0;
								count++;
							} else if (val < 64) {
								freq[val]++;
							}
						} else { // positive strand
							if ((j + 2)==(orfs[count + 2*tot])) {
								inside = 0;
								count++;
							} else if (val < 64) {
								freq[val]++;
							}
						}
					} else { // look for end
						if (s) { // negative strand
							if ((x_i.length - j - 3)==orfs[count + 2*tot]) {
								inside = 1;
							} else if (val < 64) {
								bg[val]++;
							}
						} else { // positive strand
							if ((j + 4)==orfs[count + 3*tot]) {
								inside = 1;
							} else if (val < 64) {
								bg[val]++;
							}
						}
					}
				}
			}
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	int notstop;
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			sumfreq += freq[i];
			sumbg += bg[i];
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		notstop = 1;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				notstop = 0;
				break;
			}
		}
		if (notstop) {
			rans[i] = log(((double)freq[i]/(double)sumfreq)/((double)bg[i]/(double)sumbg));
		} else {
			rans[i] = 0;
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP startCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP start_codons)
{
	int i, j, k, s, val;
	int lstarts = length(start_codons);
	int *starts = INTEGER(start_codons);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int count = 0;
	int curr_i = 0;
	
	int *freq = Calloc(64, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot] - 3; // start codon
		} else { // positive strand
			j = orfs[i + 2*tot] + 1; // start codon
		}
		
		if (s) { // negative strand
			val = getBaseRC(x_i.ptr[j++]);
			val += 4*getBaseRC(x_i.ptr[j++]);
			val += 16*getBaseRC(x_i.ptr[j]);
		} else { // positive strand
			val = getBase(x_i.ptr[j--]);
			val += 4*getBase(x_i.ptr[j--]);
			val += 16*getBase(x_i.ptr[j]);
		}
		
		if (count < l &&
			(i + 1)==index[count]) {
			count++;
			if (val < 64)
				freq[val]++;
		} else if (val < 64) {
			bg[val]++;
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	int isstart;
	for (i = 0; i < 64; i++) {
		isstart = 0;
		for (k = 0; k < lstarts; k++) {
			if (i==starts[k]) {
				isstart = 1;
				break;
			}
		}
		if (isstart) {
			if (freq[i]==0)
				freq[i] = 1; // add pseudocount
			if (bg[i]==0)
				bg[i] = 1; // add pseudocount
			sumfreq += freq[i];
			sumbg += bg[i];
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		isstart = 0;
		for (k = 0; k < lstarts; k++) {
			if (i==starts[k]) {
				isstart = 1;
				break;
			}
		}
		if (isstart) {
			rans[i] = log(((double)freq[i]/(double)sumfreq)/((double)bg[i]/(double)sumbg));
		} else {
			rans[i] = 0;
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreStartCodonModel(SEXP x, SEXP orftable, SEXP start_scores)
{
	int i, j, s, val;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *starts = REAL(start_scores);
	int curr_i = 0;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot] - 3; // start codon
		} else { // positive strand
			j = orfs[i + 2*tot] + 1; // start codon
		}
		
		if (s) { // negative strand
			val = getBaseRC(x_i.ptr[j++]);
			val += 4*getBaseRC(x_i.ptr[j++]);
			val += 16*getBaseRC(x_i.ptr[j]);
		} else { // positive strand
			val = getBase(x_i.ptr[j--]);
			val += 4*getBase(x_i.ptr[j--]);
			val += 16*getBase(x_i.ptr[j]);
		}
		
		if (val < 64) {
			rans[i] = starts[val];
		} else {
			rans[i] = 0;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP initialCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP initial_codons)
{
	int i, j, k, s, val, fin;
	int ini = asInteger(initial_codons);
	int ini_nucs = 3*ini;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int curr_i = 0;
	
	int *freq = Calloc(64*ini, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < l; i++) {
		// get the current sequence
		if (orfs[index[i] - 1] != curr_i) {
			curr_i = orfs[index[i] - 1];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start position
		s = orfs[index[i] - 1 + tot];
		if (s) { // negative strand
			j = orfs[index[i] - 1 + 3*tot] - 3 - ini_nucs; // start codon - 1 codon - ini_nucs
			fin = j;
		} else { // positive strand
			j = orfs[index[i] - 1 + 2*tot] + 1 + ini_nucs; // start codon + 1 codon + ini_nucs
			fin = j;
		}
		
		for (k = ini - 1; k >= 0; k--) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64)
				freq[64*k + val]++;
		}
		
		if (s) { // negative strand
			j = orfs[index[i] - 1 + 2*tot] + 2; // start at stop codon + 1 codon
		} else { // positive strand
			j = orfs[index[i] - 1 + 3*tot] - 4; // start at stop codon - 1 codon
		}
		
		while (j != fin) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64)
				bg[val]++;
		}
	}
	
	// normalize the frequencies
	int sumfreq[ini];
	for (k = 0; k < ini; k++) {
		sumfreq[k] = 0;
		for (i = 0; i < 64; i++) {
			sumfreq[k] += freq[64*k + i];
		}
	}
	int sumbg = 0;
	for (i = 0; i < 64; i++)
		sumbg += bg[i];
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, 64, ini));
	double *rans = REAL(ans);
	
	for (k = 0; k < ini; k++) {
		for (i = 0; i < 64; i++) {
			if (freq[64*k + i]==0 || bg[i]==0) {
				rans[64*k + i] = 0;
			} else {
				rans[64*k + i] = log(((double)freq[64*k + i]/(double)sumfreq[k])/((double)bg[i]/(double)sumbg));
			}
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreInitialCodonModel(SEXP x, SEXP orftable, SEXP ini_scores)
{
	int i, j, k, s, val;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int ini = length(ini_scores)/64;
	int ini_nucs = 3*ini;
	double *scores = REAL(ini_scores);
	int curr_i = 0;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot] - 3 - ini_nucs; // start codon - 1 codon
		} else { // positive strand
			j = orfs[i + 2*tot] + 1 + ini_nucs; // start codon + 1 codon
		}
		
		rans[i] = 0;
		
		for (k = ini - 1; k >= 0; k--) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64) {
				rans[i] += scores[64*k + val];
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP terminationCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP terminal_codons)
{
	int i, j, k, s, val, fin;
	int ter = asInteger(terminal_codons);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int curr_i = 0;
	
	int *freq = Calloc(64*ter, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < l; i++) {
		// get the current sequence
		if (orfs[index[i] - 1] != curr_i) {
			curr_i = orfs[index[i] - 1];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the end position
		s = orfs[index[i] - 1 + tot];
		if (s) { // negative strand
			j = orfs[index[i] - 1 + 2*tot] + 2; // stop codon + 1 codon
			fin = orfs[index[i] - 1 + 3*tot] - 3; // start codon - 1 codon
		} else { // positive strand
			j = orfs[index[i] - 1 + 3*tot] - 4; // stop codon - 1 codon
			fin = orfs[index[i] - 1 + 2*tot] + 1; // start codon + 1 codon
		}
		
		for (k = 0; k < ter; k++) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64)
				freq[64*k + val]++;
		}
		
		while (j != fin) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64)
				bg[val]++;
		}
	}
	
	// normalize the frequencies
	int sumfreq[ter];
	for (k = 0; k < ter; k++) {
		sumfreq[k] = 0;
		for (i = 0; i < 64; i++) {
			sumfreq[k] += freq[64*k + i];
		}
	}
	int sumbg = 0;
	for (i = 0; i < 64; i++)
		sumbg += bg[i];
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, 64, ter));
	double *rans = REAL(ans);
	
	for (k = 0; k < ter; k++) {
		for (i = 0; i < 64; i++) {
			if (freq[64*k + i]==0 || bg[i]==0) {
				rans[64*k + i] = 0;
			} else {
				rans[64*k + i] = log(((double)freq[64*k + i]/(double)sumfreq[k])/((double)bg[i]/(double)sumbg));
			}
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreTerminationCodonModel(SEXP x, SEXP orftable, SEXP ter_scores)
{
	int i, j, k, s, val;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int ter = length(ter_scores)/64;
	double *scores = REAL(ter_scores);
	int curr_i = 0;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 2*tot] + 2; // stop codon + 1 codon
		} else { // positive strand
			j = orfs[i + 3*tot] - 4; // stop codon - 1 codon
		}
		
		rans[i] = 0;
		
		for (k = 0; k < ter; k++) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64) {
				rans[i] += scores[64*k + val];
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP getRegion(SEXP x, SEXP orftable, SEXP width, SEXP offset)
{
	int i, j, k, up, s;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int w = asInteger(width);
	if (w < 0) { // downstream
		w *= -1;
		up = 0;
	} else { // upstream
		up = 1;
	}
	int curr_i = 0;
	int off = asInteger(offset);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, tot));
	char seq[w + 1]; // each sequence
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the strand
		s = orfs[i + tot];
		
		if (up) { // upstream
			if (s) { // negative strand
				j = orfs[i + 3*tot] + w - 1 - off;
			} else { // positive strand
				j = orfs[i + 2*tot] - w - 1 + off;
			}
			
			if ((s==1 && j < x_i.length) ||
				(s==0 && j >= 0)) {
				for (k = 0; k < w; k++) {
					if (s) { // negative strand
						seq[k] = getBaseLetterRC(x_i.ptr[j--]);
					} else { // positive strand
						seq[k] = getBaseLetter(x_i.ptr[j++]);
					}
				}
			} else {
				k = 0;
			}
		} else { // downstream
			if (s) { // negative strand
				j = orfs[i + 3*tot] - 4 - off;
			} else { // positive strand
				j = orfs[i + 2*tot] + 2 + off;
			}
			
			if ((s==1 && (j + w - 1) < x_i.length) ||
				(s==0 && (j - w + 1) >= 0)) {
				for (k = 0; k < w; k++) {
					if (s) { // negative strand
						seq[k] = getBaseLetterRC(x_i.ptr[j--]);
					} else { // positive strand
						seq[k] = getBaseLetter(x_i.ptr[j++]);
					}
				}
			} else {
				k = 0;
			}
		}
		
		seq[k] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(seq));
	}
	
	UNPROTECT(1);
	
	return seqs;
}

SEXP autocorrelationModel(SEXP x, SEXP orftable, SEXP indices, SEXP aatable)
{
	int i, j, fin, s, val, lastVal, dist;
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int tot = length(orftable)/4; // number of rows
	int *AAs = INTEGER(aatable);
	int l = length(indices);
	int curr_i = 0;
	
	int *freq = Calloc(4096, int);
	int *prev = Calloc(20, int);
	int *pos = Calloc(20, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < l; i++) {
		// get the current sequence
		if (orfs[index[i] - 1] != curr_i) {
			curr_i = orfs[index[i] - 1];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		for (j = 0; j < 20; j++) {
			prev[j] = 100000;
			pos[j] = -100000;
		}
		
		// get the start position
		s = orfs[index[i] - 1 + tot];
		if (s) { // negative strand
			fin = orfs[index[i] - 1 + 3*tot] - 3; // finish at start codon - 1 codon
			j = orfs[index[i] - 1 + 2*tot] + 2; // start at stop codon + 1 codon
		} else { // positive strand
			fin = orfs[index[i] - 1 + 2*tot] + 1; // finish at start codon + 1 codon
			j = orfs[index[i] - 1 + 3*tot] - 4; // start at stop codon - 1 codon
		}
		
		while (1) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64) {
				lastVal = prev[AAs[val]];
				dist = j - pos[AAs[val]];
				if (dist < 0)
					dist *= -1;
				if (val < 64 &&
					lastVal < 64 &&
					dist > 1 && // not already incorporated into dicodon scores
					dist < 20) // last observation was reasonably closeby
					freq[lastVal*64 + val]++;
				prev[AAs[val]] = val;
				pos[AAs[val]] = j;
			}
			
			if (j==fin)
				break;
		}
	}
	
	Free(prev);
	Free(pos);
	
	// calculate the row sums
	int *rowSums = Calloc(64, int);
	int *colSums = Calloc(64, int);
	int *AAsums = Calloc(20, int);
	for (i = 0; i < 64; i++) {
		for (j = 0; j < 64; j++) {
			if (freq[64*j + i] > 0) {
				rowSums[i] += freq[64*j + i];
				colSums[j] += freq[64*j + i];
				AAsums[AAs[j]] += freq[64*j + i];
			}
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, 64, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		for (j = 0; j < 64; j++) {
			if (freq[64*j + i]==0 || colSums[j]==0 || rowSums[i]==0) {
				rans[64*j + i] = 0;
			} else {
				rans[64*j + i] = log(((double)freq[64*j + i]/(double)rowSums[i])/((double)colSums[j]/(double)AAsums[AAs[j]]));
			}
		}
	}
	
	Free(freq);
	Free(rowSums);
	Free(colSums);
	Free(AAsums);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreAutocorrelationModel(SEXP x, SEXP orftable, SEXP codon_scores, SEXP aatable)
{
	int i, j, s, fin, val, lastVal, dist;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *codons = REAL(codon_scores);
	int *AAs = INTEGER(aatable);
	int curr_i = 0;
	double score;
	
	int *prev = Calloc(20, int);
	int *pos = Calloc(20, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot;) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		for (j = 0; j < 20; j++) {
			prev[j] = 100000;
			pos[j] = -100000;
		}
		
		// get the start and stop positions
		s = orfs[i + tot];
		if (s) { // negative strand
			fin = orfs[i + 3*tot] - 3; // finish at start codon - 1 codon
			j = orfs[i + 2*tot] + 2; // start at stop codon + 1 codon
		} else { // positive strand
			fin = orfs[i + 2*tot] + 1; // finish at start codon + 1 codon
			j = orfs[i + 3*tot] - 4; // start at stop codon - 1 codon
		}
		
		// sum score from j to fin
		score = 0;
		while (1) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64) {
				lastVal = prev[AAs[val]];
				dist = j - pos[AAs[val]];
				if (dist < 0)
					dist *= -1;
				if (val < 64 &&
					lastVal < 64 &&
					dist > 1 && // not already incorporated into dicodon scores
					dist < 20) // last observation was reasonably closeby
					score += codons[lastVal*64 + val];
				prev[AAs[val]] = val;
				pos[AAs[val]] = j;
			}
			
			if (j==fin) {
				// record the score
				rans[i] = score;
				i++;
				
				// check whether still within bounds
				if (s) { // negative strand
					if (i != tot &&
						orfs[i + 2*tot]==orfs[i + 2*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 3*tot] - 3;
					} else {
						break;
					}
				} else { // positive strand
					if (i != tot &&
						orfs[i + 3*tot]==orfs[i + 3*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 2*tot] + 1;
					} else {
						break;
					}
				}
			}
		}
	}
	
	Free(prev);
	Free(pos);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP nucleotideBiasModel(SEXP x, SEXP orftable, SEXP indices, SEXP positions)
{
	int i, j, p, gene, s, val;
	int pos = asInteger(positions);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int count = 0;
	int curr_i = 0;
	
	int *freq = Calloc(4*pos, int);
	int *bg = Calloc(4*pos, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		if (count < l &&
			(i + 1)==index[count]) {
			gene = 1;
		} else {
			gene = 0;
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot];
			if (j + pos > x_i.length) {
				if (gene)
					count++;
				continue;
			}
		} else { // positive strand
			j = orfs[i + 2*tot] - 2;
			if (j - pos < -1) {
				if (gene)
					count++;
				continue;
			}
		}
		
		for (p = 0; p < pos; p++) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
			}
			if (val < 4) {
				if (gene) {
					freq[p*4 + val]++;
				} else {
					bg[p*4 + val]++;
				}
			}
		}
		
		if (gene)
			count++;
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, 4, pos));
	double *rans = REAL(ans);
	
	// normalize the frequencies
	int sumfreq, sumbg;
	for (p = 0; p < pos; p++) {
		sumfreq = 0;
		sumbg = 0;
		for (i = 0; i < 4; i++) {
			if (freq[p*4 + i]==0)
				freq[p*4 + i] = 1; // add pseudocount
			if (bg[p*4 + i]==0)
				bg[p*4 + i] = 1; // add pseudocount
			sumfreq += freq[p*4 + i];
			sumbg += bg[p*4 + i];
		}
		for (i = 0; i < 4; i++)
			rans[p*4 + i] = log(((double)freq[p*4 + i]/(double)sumfreq)/((double)bg[p*4 + i]/(double)sumbg));
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreNucleotideBiasModel(SEXP x, SEXP orftable, SEXP nuc_scores)
{
	int i, j, p, s, val;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *nucs = REAL(nuc_scores);
	int pos = length(nuc_scores)/4;
	int curr_i = 0;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		rans[i] = 0;
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot];
			if (j + pos > x_i.length)
				continue;
		} else { // positive strand
			j = orfs[i + 2*tot] - 2;
			if (j - pos < -1)
				continue;
		}
		
		for (p = 0; p < pos; p++) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
			}
			if (val < 4)
				rans[i] += nucs[p*4 + val];
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP upstreamMotifModel(SEXP x, SEXP orftable, SEXP indices, SEXP begin, SEXP distance, SEXP size)
{
	int i, j, k, p, gene, s, val;
	int beg = asInteger(begin);
	int d = asInteger(distance);
	int kmer = asInteger(size); // k-mer size
	int n = pow(4, kmer);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int count = 0;
	int curr_i = 0;
	
	int *freq = Calloc(n, int);
	int *bg = Calloc(n, int);
	
	int mult[kmer];
	mult[0] = 1;
	for (i = 1; i < kmer; i++)
		mult[i] = 4*mult[i - 1];
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		if (count < l &&
			(i + 1)==index[count]) {
			gene = 1;
		} else {
			gene = 0;
		}
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot] + beg - 1;
			if (j + d > x_i.length) {
				if (gene)
					count++;
				continue;
			}
		} else { // positive strand
			j = orfs[i + 2*tot] - beg - 1;
			if (j - d < -1) {
				if (gene)
					count++;
				continue;
			}
		}
		
		for (p = 0; p < d - kmer + 1; p++) {
			val = 0;
			for (k = 1; k <= kmer; k++) {
				if (s) { // negative strand
					val += mult[k - 1]*getBaseRC(x_i.ptr[j + k - 1]);
				} else { // positive strand
					val += mult[k - 1]*getBase(x_i.ptr[j - k + 1]);
				}
			}
			if (s) { // negative strand
				j++;
			} else { // positive strand
				j--;
			}
			
			if (val < n) {
				if (gene) {
					freq[val]++;
				} else {
					bg[val]++;
				}
			}
		}
		
		if (gene)
			count++;
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, n));
	double *rans = REAL(ans);
	
	// normalize the frequencies
	int sumfreq = 0, sumbg = 0;
	for (i = 0; i < n; i++) {
		if (freq[i]==0)
			freq[i] = 1; // add pseudocount
		if (bg[i]==0)
			bg[i] = 1; // add pseudocount
		sumfreq += freq[i];
		sumbg += bg[i];
	}
	for (i = 0; i < n; i++)
		rans[i] = log(((double)freq[i]/(double)sumfreq)/((double)bg[i]/(double)sumbg));
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreUpstreamMotifModel(SEXP x, SEXP orftable, SEXP motif_scores, SEXP begin, SEXP distance, SEXP size)
{
	int i, j, k, p, s, val;
	int beg = asInteger(begin);
	int d = asInteger(distance);
	int kmer = asInteger(size); // k-mer size
	int n = pow(4, kmer);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *motif = REAL(motif_scores);
	int curr_i = 0;
	
	int mult[kmer];
	mult[0] = 1;
	for (i = 1; i < kmer; i++)
		mult[i] = 4*mult[i - 1];
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		rans[i] = 0;
		
		// get the start position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 3*tot] + beg - 1;
			if (j + d > x_i.length)
				continue;
		} else { // positive strand
			j = orfs[i + 2*tot] - beg - 1;
			if (j - d < -1)
				continue;
		}
		
		for (p = 0; p < d - kmer + 1; p++) {
			val = 0;
			for (k = 1; k <= kmer; k++) {
				if (s) { // negative strand
					val += mult[k - 1]*getBaseRC(x_i.ptr[j + k - 1]);
				} else { // positive strand
					val += mult[k - 1]*getBase(x_i.ptr[j - k + 1]);
				}
			}
			if (s) { // negative strand
				j++;
			} else { // positive strand
				j--;
			}
			
			if (val < n)
				rans[i] += motif[val]/kmer;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

int setRun(int runLength, double score, int *freq)
{
	if (score > 0) {
		if (runLength >= 0) { // elongate run
			runLength++;
		} else { // switch sign
			if (runLength < -20)
				runLength = -20;
			freq[-1 - runLength]++;
			runLength = 1;
		}
	} else if (score < 0) {
		if (runLength <= 0) { // elongate run
			runLength--;
		} else { // switch sign
			if (runLength > 20)
				runLength = 20;
			freq[runLength - 1]++;
			runLength = -1;
		}
	} // skip if score == 0
	return runLength;
}

SEXP runLengthModel(SEXP x, SEXP orftable, SEXP codon_scores)
{
	int i, j, rf, s, val, lastVal, x_length;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *codons = REAL(codon_scores);
	int dicodon = (length(codon_scores)==64) ? 0 : 1;
	int inside = 0;
	int runLength;
	
	int *freq = Calloc(20, int);
	int *bg = Calloc(20, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	int count = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (s = 0; s <= 1; s++) { // strand
			for (rf = 1; rf <= 3; rf++) { // reading frame
				lastVal = 100000;
				runLength = 0;
				for (j = x_i.length - rf; j >= 2;) { // position
					if (s) { // negative strand
						val = getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 4*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
						val += 16*getBaseRC(x_i.ptr[x_i.length - j - 1]);
						j--;
					} else { // positive strand
						val = getBase(x_i.ptr[j--]);
						val += 4*getBase(x_i.ptr[j--]);
						val += 16*getBase(x_i.ptr[j--]);
					}
					
					if (inside) { // look for beg
						if (s) { // negative strand
							if ((x_i.length - j - 1)==(orfs[count + 3*tot])) {
								inside = 0;
								lastVal = 100000;
								count++;
								runLength = 0;
							} else if (val < 64) {
								if (dicodon && lastVal < 64) {
									runLength = setRun(runLength, codons[lastVal*64 + val], freq);
								} else {
									runLength = setRun(runLength, codons[val], freq);
								}
							}
						} else { // positive strand
							if ((j + 2)==(orfs[count + 2*tot])) {
								inside = 0;
								lastVal = 100000;
								count++;
								runLength = 0;
							} else if (val < 64) {
								if (dicodon && lastVal < 64) {
									runLength = setRun(runLength, codons[lastVal*64 + val], freq);
								} else {
									runLength = setRun(runLength, codons[val], freq);
								}
							}
						}
					} else { // look for end
						if (s) { // negative strand
							if ((x_i.length - j - 3)==orfs[count + 2*tot]) {
								inside = 1;
								lastVal = 100000;
								runLength = 0;
							} else if (val < 64) {
								if (dicodon && lastVal < 64) {
									runLength = setRun(runLength, codons[lastVal*64 + val], bg);
								} else {
									runLength = setRun(runLength, codons[val], bg);
								}
							}
						} else { // positive strand
							if ((j + 4)==orfs[count + 3*tot]) {
								inside = 1;
								lastVal = 100000;
								runLength = 0;
							} else if (val < 64) {
								if (dicodon && lastVal < 64) {
									runLength = setRun(runLength, codons[lastVal*64 + val], bg);
								} else {
									runLength = setRun(runLength, codons[val], bg);
								}
							}
						}
					}
					
					lastVal = val;
				}
			}
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	for (i = 0; i < 20; i++) {
		if (freq[i]==0)
			freq[i] = 1; // add pseudocounts
		if (bg[i]==0)
			bg[i] = 1; // add pseudocounts
		sumfreq += freq[i];
		sumbg += bg[i];
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 20));
	double *rans = REAL(ans);
	
	for (i = 0; i < 20; i++)
		rans[i] = log((double)freq[i]/(double)sumfreq/((double)bg[i]/(double)sumbg));
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreRunLengthModel(SEXP x, SEXP orftable, SEXP codon_scores, SEXP run_scores)
{
	int i, j, s, fin, val, lastVal;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *codons = REAL(codon_scores);
	double *runs = REAL(run_scores);
	int curr_i = 0;
	double score;
	int runLength;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot;) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the start and stop positions
		s = orfs[i + tot];
		if (s) { // negative strand
			fin = orfs[i + 3*tot] - 3; // finish at start codon - 1 codon
			j = orfs[i + 2*tot] + 2; // start at stop codon + 1 codon
		} else { // positive strand
			fin = orfs[i + 2*tot] + 1; // finish at start codon + 1 codon
			j = orfs[i + 3*tot] - 4; // start at stop codon - 1 codon
		}
		
		// sum score from j to fin
		score = 0;
		runLength = 0;
		lastVal = 100000;
		while (1) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64 && lastVal < 64) {
				if (codons[lastVal*64 + val] > 0) {
					if (runLength >= 0) { // elongate run
						runLength++;
					} else { // switch sign
						if (runLength < -20)
							runLength = -20;
						score += runs[-1 - runLength];
						runLength = 1;
					}
				} else {
					if (runLength <= 0) { // elongate run
						runLength--;
					} else { // switch sign
						if (runLength > 20)
							runLength = 20;
						score += runs[runLength - 1];
						runLength = -1;
					}
				} // assume never == 0 inside gene
			}
			lastVal = val;
			
			if (j==fin) {
				// record the score
				rans[i] = score;
				// include score for any remaining runs
				if (runLength > 20 || runLength < 20) {
					rans[i] += runs[19];
				} else if (runLength > 0) {
					rans[i] += runs[runLength - 1];
				} else if (runLength < 0) {
					rans[i] += runs[-1 - runLength];
				}
				i++;
				
				// check whether still within bounds
				if (s) { // negative strand
					if (i != tot &&
						orfs[i + 2*tot]==orfs[i + 2*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 3*tot] - 3;
					} else {
						break;
					}
				} else { // positive strand
					if (i != tot &&
						orfs[i + 3*tot]==orfs[i + 3*tot - 1] &&
						s==orfs[i + tot] &&
						curr_i==orfs[i]) {
						// continue scoring
						fin = orfs[i + 2*tot] + 1;
					} else {
						break;
					}
				}
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP stopCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP stop_codons)
{
	int i, j, k, s, val;
	int lstops = length(stop_codons);
	int *stops = INTEGER(stop_codons);
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int count = 0;
	int curr_i = 0;
	
	int *freq = Calloc(64, int);
	int *bg = Calloc(64, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the stop position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 2*tot] - 1; // stop codon
		} else { // positive strand
			j = orfs[i + 3*tot] - 1; // stop codon
		}
		
		if (s) { // negative strand
			val = getBaseRC(x_i.ptr[j++]);
			val += 4*getBaseRC(x_i.ptr[j++]);
			val += 16*getBaseRC(x_i.ptr[j]);
		} else { // positive strand
			val = getBase(x_i.ptr[j--]);
			val += 4*getBase(x_i.ptr[j--]);
			val += 16*getBase(x_i.ptr[j]);
		}
		
		if (count < l &&
			(i + 1)==index[count]) {
			count++;
			if (val < 64)
				freq[val]++;
		} else if (val < 64) {
			bg[val]++;
		}
	}
	
	// normalize the frequencies
	int sumfreq = 0;
	int sumbg = 0;
	int isstop;
	for (i = 0; i < 64; i++) {
		isstop = 0;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				isstop = 1;
				break;
			}
		}
		if (isstop) {
			if (freq[i]==0)
				freq[i] = 1; // add pseudocount
			if (bg[i]==0)
				bg[i] = 1; // add pseudocount
			sumfreq += freq[i];
			sumbg += bg[i];
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 64));
	double *rans = REAL(ans);
	
	for (i = 0; i < 64; i++) {
		isstop = 0;
		for (k = 0; k < lstops; k++) {
			if (i==stops[k]) {
				isstop = 1;
				break;
			}
		}
		if (isstop) {
			rans[i] = log(((double)freq[i]/(double)sumfreq)/((double)bg[i]/(double)sumbg));
		} else {
			rans[i] = 0;
		}
	}
	
	Free(freq);
	Free(bg);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP scoreStopCodonModel(SEXP x, SEXP orftable, SEXP stop_scores)
{
	int i, j, s, val;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	double *stops = REAL(stop_scores);
	int curr_i = 0;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, tot));
	double *rans = REAL(ans);
	
	// loop through each candidate ORF
	for (i = 0; i < tot; i++) {
		// get the current sequence
		if (orfs[i] != curr_i) {
			curr_i = orfs[i];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the stop position
		s = orfs[i + tot];
		if (s) { // negative strand
			j = orfs[i + 2*tot] - 1; // stop codon
		} else { // positive strand
			j = orfs[i + 3*tot] - 1; // stop codon
		}
		
		if (s) { // negative strand
			val = getBaseRC(x_i.ptr[j++]);
			val += 4*getBaseRC(x_i.ptr[j++]);
			val += 16*getBaseRC(x_i.ptr[j]);
		} else { // positive strand
			val = getBase(x_i.ptr[j--]);
			val += 4*getBase(x_i.ptr[j--]);
			val += 16*getBase(x_i.ptr[j]);
		}
		
		if (val < 64) {
			rans[i] = stops[val];
		} else {
			rans[i] = 0;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// gene codon frequencies excluding start and stop codon
SEXP codonFrequencies(SEXP x, SEXP orftable, SEXP indices)
{
	int i, j, s, val, fin;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int *index = INTEGER(indices);
	int l = length(indices);
	int curr_i = 0;
	
	int *freq = Calloc(64*l, int);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	
	// loop through each candidate ORF
	for (i = 0; i < l; i++) {
		// get the current sequence
		if (orfs[index[i] - 1] != curr_i) {
			curr_i = orfs[index[i] - 1];
			x_i = get_elt_from_XStringSet_holder(&x_set, curr_i - 1);
		}
		
		// get the end position
		s = orfs[index[i] - 1 + tot];
		if (s) { // negative strand
			j = orfs[index[i] - 1 + 2*tot] + 2; // stop codon + 1 codon
			fin = orfs[index[i] - 1 + 3*tot] - 3; // start codon - 1 codon
		} else { // positive strand
			j = orfs[index[i] - 1 + 3*tot] - 4; // stop codon - 1 codon
			fin = orfs[index[i] - 1 + 2*tot] + 1; // start codon + 1 codon
		}
		
		while (j != fin) {
			if (s) { // negative strand
				val = getBaseRC(x_i.ptr[j++]);
				val += 4*getBaseRC(x_i.ptr[j++]);
				val += 16*getBaseRC(x_i.ptr[j++]);
			} else { // positive strand
				val = getBase(x_i.ptr[j--]);
				val += 4*getBase(x_i.ptr[j--]);
				val += 16*getBase(x_i.ptr[j--]);
			}
			
			if (val < 64)
				freq[l*val + i]++;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, l, 64));
	double *rans = REAL(ans);
	
	// normalize the frequencies
	for (i = 0; i < l; i++) {
		int sumfreq = 0;
		for (j = 0; j < 64; j++)
			sumfreq += freq[l*j + i];
		if (sumfreq > 0) {
			for (j = 0; j < 64; j++)
				rans[l*j + i] = (double)freq[l*j + i]/(double)sumfreq;
		} else {
			for (j = 0; j < 64; j++)
				rans[l*j + i] = 0;
		}
	}
	
	Free(freq);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP chainGenes(SEXP orftable, SEXP topScore, SEXP topLength, SEXP scoreIntergenic, SEXP maxOverlapSame, SEXP maxOverlapOpposite, SEXP maxFracOverlap, SEXP sameScores, SEXP oppoScores)
{
	int i, j;
	int tot = length(orftable)/4; // number of rows
	if (tot==0)
		return NEW_INTEGER(0);
	int *orfs = INTEGER(orftable);
	int *topL = INTEGER(topLength);
	int scoreInt = asInteger(scoreIntergenic);
	double maxOvS = asReal(maxOverlapSame);
	double maxOvO = asReal(maxOverlapOpposite);
	double maxOvF = asReal(maxFracOverlap);
	double *samS = REAL(sameScores);
	double *oppS = REAL(oppoScores);
	double *topS = REAL(topScore);
	int l = length(sameScores);
	int maxL = (l - 1)/2;
	int negMaxL = -1*maxL;
	l--;
	
//	SEXP pointers; // index of previous row in chain
//	PROTECT(pointers = allocVector(REALSXP, tot));
//	double *p = REAL(pointers);
//	for (i = 0; i < tot; i++)
//		p[i] = i + 1; // starts at 1
	
//	SEXP cumscore; // cumulative score of chain
//	PROTECT(cumscore = duplicate(topScore));
//	double *c = REAL(cumscore);
	
	int *p = Calloc(tot, int);
	double *c = Calloc(tot, double);
	for (i = 0; i < tot; i++) {
		p[i] = i;
		c[i] = topS[i];
	}
	
	int currIndex = orfs[0]; // current index in orfs
	int switchIndex = 0; // i at last index switch
	int lastPointer = -1; // pointer to top cumscore in last index
	int counter = -1; // last deactivated row
	
	int pointer; // pointer to top index in chain
	double currscore; // minimum score required to chain
	int firstCouldNotChain; // pointer to first instance of !canChain
	int sameStrand; // whether indices are on the same strand
	int delta; // difference in positions
	int canChain; // whether an index can be added to the chain
	double temp; // temporary score
	
	for (i = 0; i < tot; i++) {
		if (orfs[i] != currIndex) {
			currIndex = orfs[i];
			counter = i - 1;
			temp = 0;
			for (j = switchIndex; j <= counter; j++) {
				if (c[j] > temp) {
					lastPointer = j;
					temp = c[j];
				}
			}
			switchIndex = i;
		}
		
		pointer = lastPointer;
		if (pointer >= 0) {
			currscore = c[pointer];
		} else {
			currscore = 0;
		}
		
		j = counter + 1;
		firstCouldNotChain = -1;
		
		while (j < i) {
			sameStrand = orfs[i + tot] == orfs[j + tot];
			delta = orfs[i + 3*tot] - orfs[j + 3*tot]; // deltaend
			canChain = delta > 0 && (orfs[i + 2*tot] != orfs[j + 2*tot]);
			
			if (canChain) {
				delta = orfs[j + 3*tot] - orfs[i + 2*tot] + 1; //deltastart
				if (sameStrand) {
					canChain = delta <= maxOvS;
				} else {
					canChain = delta <= maxOvO;
				}
				if (canChain && delta > 0) {
					canChain = (double)delta/(double)topL[i] < maxOvF &&
						(double)delta/(double)topL[j] < maxOvF;
				}
			}
			
			if (canChain) { // possible chain
				if (scoreInt) {
					if (sameStrand) {
						if (delta < negMaxL) {
							temp = samS[l];
						} else {
							temp = samS[maxL - delta];
						}
					} else {
						if (delta < negMaxL) {
							temp = oppS[l];
						} else {
							temp = oppS[maxL - delta];
						}
					}
				} else {
					temp = 0;
				}
				temp += c[j];
				
				if (temp > currscore) {
					currscore = temp;
					pointer = j;
					
					if (firstCouldNotChain >= 0) {
						counter = firstCouldNotChain - 1;
					} else {
						counter = j - 1;
					}
				}
			} else if (firstCouldNotChain==-1) {
				firstCouldNotChain = j;
			}
			j++;
		}
		
		if (pointer >= 0) {
			c[i] += currscore;
			p[i] = pointer;
		}
	}
	
//	SEXP ret_list;
//	PROTECT(ret_list = allocVector(VECSXP, 2));
//	SET_VECTOR_ELT(ret_list, 0, pointers);
//	SET_VECTOR_ELT(ret_list, 1, cumscore);
	
//	return ret_list;
	
	pointer = 0;
	for (i = 1; i < tot; i++)
		if (c[i] > c[pointer])
			pointer = i;
	
	int *indices = Calloc(pointer + 1, int);
	int position = 0;
	indices[position] = pointer;
	while (pointer != p[pointer]) {
		pointer = p[pointer];
		position++;
		indices[position] = pointer;
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, position + 1));
	int *rans = INTEGER(ans);
	i = 0;
	while (position >= 0)
		rans[i++] = indices[position--] + 1;
	
	Free(p);
	Free(c);
	Free(indices);
	
//	UNPROTECT(3);
	UNPROTECT(1);
	
	return ans;
}

SEXP longestORFs(SEXP orftable)
{
	int i, lastVal;
	int tot = length(orftable)/4; // number of rows
	int *orfs = INTEGER(orftable);
	int count = 0;
	int lastI = 0;
	int *longest = Calloc(tot, int);
	
	// initialize lastVal from first row
	if (orfs[tot]) {
		lastVal = orfs[2*tot];
	} else {
		lastVal = orfs[3*tot];
	}
	
	// loop through each candidate ORF
	for (i = 1; i < tot; i++) {
		if ((orfs[i] != orfs[lastI]) ||
			(orfs[i + tot] != orfs[lastI + tot])) {
			if (orfs[i + tot]) {
				lastVal = orfs[i + 2*tot];
			} else {
				lastVal = orfs[i + 3*tot];
			}
			longest[lastI] = 1;
			count++;
		} else {
			if (orfs[i + tot]) {
				if (orfs[i + 2*tot] != lastVal) {
					longest[lastI] = 1;
					lastVal = orfs[i + 2*tot];
					count++;
				}
			} else {
				if (orfs[i + 3*tot] != lastVal) {
					longest[lastI] = 1;
					lastVal = orfs[i + 3*tot];
					count++;
				}
			}
		}
		lastI = i;
	}
	longest[lastI] = 1;
	count++;
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, count));
	int *rans = INTEGER(ans);
	
	count = 0;
	for (i = 0; i < tot; i++)
		if (longest[i])
			rans[count++] = i + 1;
	
	Free(longest);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP getIndex(SEXP start1, SEXP start2, SEXP len, SEXP score)
{
	int i, j, e;
	int count = 0;
	int l = asInteger(len);
	int n = length(start1);
	int *s1 = INTEGER(start1);
	int *s2 = INTEGER(start2);
	int *s = INTEGER(score);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	for (j = 0; j < l; j++)
		rans[j] = 0;
	
	for (j = 0; j < n; j++) {
		if (s2[j] > count) {
			if (s1[j] > l) {
				break;
			} else if (s1[j] > count) {
				if (s2[j] > l) {
					e = l;
				} else {
					e = s2[j];
				}
				for (i = s1[j] - 1; i < e; i++)
					if (s[j] > rans[i])
						rans[i] = s[j];
				count = s2[j];
			} else if (count < s2[j]) {
				if (s2[j] > l) {
					e = l;
				} else {
					e = s2[j];
				}
				for (i = count; i < e; i++)
					if (s[j] > rans[i])
						rans[i] = s[j];
				count = s2[j];
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP inBounds(SEXP vec1, SEXP vec3, SEXP lo1, SEXP hi1, SEXP hi3)
{
	int i = 0, j, count = 0;
	int l = length(vec1);
	int *v1 = INTEGER(vec1);
	double *v3 = REAL(vec3);
	
	int l1 = asInteger(lo1);
	int h1 = asInteger(hi1);
	double h3 = asReal(hi3);
	
	while (i < l) {
		if (v1[i] >= l1) {
			j = i;
			while (j < l && v1[j] <= h1) {
				if (v3[j] <= h3)
					count++;
				j++;
			}
			break;
		}
		i++;
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, count));
	int *rans = INTEGER(ans);
	
	j = 0;
	while (j < count) {
		if (v1[i] >= l1 &&
			v1[i] <= h1 &&
			v3[i] <= h3)
			rans[j++] = i + 1;
		i++;
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP getBounds(SEXP widths, SEXP start, SEXP end, SEXP minL, SEXP maxL, SEXP lenScores, SEXP kmer, SEXP Ksize, SEXP negOk, SEXP minS, SEXP partS)
{
	int j, e;
	double score, temp;
	int i = 1;
	int s = 0;
	int index = 0;
	int strand = 0;
	
	int *w = INTEGER(widths);
	double *starts = REAL(start);
	double *ends = REAL(end);
	double *kmers = REAL(kmer);
	int minLength = asInteger(minL);
	int maxLength = asInteger(maxL);
	double *lS = REAL(lenScores);
	int k = asInteger(Ksize);
	int neg = asInteger(negOk);
	double minScore = asReal(minS);
	double partScore = asReal(partS);
	
	int midLength = 0;
	for (i = 0; i < maxLength - minLength + 1; i++) {
			if (lS[i] > lS[midLength])
				midLength = i;
	}
	midLength += minLength;
	
	int lX = length(widths);
	int l = length(start);
	
	int *ws = Calloc(lX, int);
	int *offsets = Calloc(lX, int);
	int *rev_width = Calloc(lX, int);
	for (j = 0; j < lX; j++) {
		rev_width[lX - j - 1] = w[j];
		if (j > 0) {
			ws[j] = ws[j - 1] + w[j];
			offsets[j] = ws[j - 1];
		} else {
			ws[j] = w[j];
		}
	}
	
	int size = 9e2;
	double *r = Calloc(size, double);
	int count = 0;
	
	while (i < l) {
		if (i >= ws[index]) {
			index++;
			s = -1;
			if (index >= lX) {
				index = 0;
				strand = 1;
				for (j = 0; j < lX; j++) {
					if (j > 0) {
						offsets[j] = ws[j - 1];
					} else {
						offsets[j] = ws[lX - 1];
					}
					ws[j] = offsets[j] + rev_width[j];
				}
			}
		} else if (starts[i - 1] < starts[i]) { // ascended
			s = i;
		} else if (s >= 0 && starts[i - 1] > starts[i]) { // descended
			for (j = s; j < i; j++)
				if (kmers[s] > kmers[j])
					s = j; // choose min
			
			e = s + minLength - 1;
			if (e >= ws[index]) { // all out of bounds
				i = ws[index];
			} else if (starts[s] >= 0) {
				score = ends[e];
				score += kmers[e - k + 2] - kmers[s];
				score += lS[e - s - minLength + 1];
				for (j = e + 1; j < ws[index] && j < s + maxLength; j++) {
					if (ends[j] >= 0 &&
						starts[s] + ends[j] >= partScore) {
						temp = ends[j];
						if (j - s + 1 < midLength) { // add k-mer score
							temp += kmers[j - k + 2] - kmers[s];
						} else { // add average k-mer score for midLength
							temp += (kmers[j - k + 2] - kmers[s])*midLength/(j - s + 1);
						}
						temp += lS[j - s - minLength + 1];
						if (score < temp) {
							e = j;
							score = temp;
						}
					}
				}
				score += starts[s];
				if (ends[e] >= 0 &&
					starts[s] + ends[e] >= partScore &&
					score >= minScore) {
					if (count + 9 >= size) {
						size *= 2;
						r = Realloc(r, size, double);
					}
					
					if (strand) {
						r[count++] = lX - index;
					} else {
						r[count++] = index + 1;
					}
					r[count++] = strand;
					
					// coordinates in XString
					r[count++] = s + 1;
					r[count++] = e + 1;
					
					r[count++] = score;
					r[count++] = neg;
					
					// coordinates to XStringSet
					if (strand) {
						r[count++] = rev_width[index] - e + offsets[index];
						r[count++] = rev_width[index] - s + offsets[index];
					} else {
						r[count++] = s - offsets[index] + 1;
						r[count++] = e - offsets[index] + 1;
					}
					
					r[count++] = starts[s] + ends[e];
				}
				
				s = -1;
			} else {
				s = -1;
			}
		}
		
		i++;
	}
	
	Free(rev_width);
	Free(offsets);
	Free(ws);
	
	int nrow = count/9;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, nrow, 9));
	double *rans = REAL(ans);
	
	i = 0;
	j = 0;
	while (i < count) {
		rans[j] = r[i++];
		rans[nrow + j] = r[i++];
		rans[2*nrow + j] = r[i++];
		rans[3*nrow + j] = r[i++];
		rans[4*nrow + j] = r[i++];
		rans[5*nrow + j] = r[i++];
		rans[6*nrow + j] = r[i++];
		rans[7*nrow + j] = r[i++];
		rans[8*nrow + j] = r[i++];
		j++;
	}
	
	Free(r);
	UNPROTECT(1);
	
	return ans;
}

SEXP addIfElse(SEXP vec, SEXP index, SEXP scores)
{
	if (MAYBE_SHARED(vec))
		error(".Call function 'addIfElse' called in incorrect context.");
	
	double *v = REAL(vec);
	double *s = REAL(scores);
	int *ind = INTEGER(index);
	int l = length(vec);
	
	for (int i = 0; i < l; i++)
		v[i] += s[ind[i]];
	
	return vec;
}

SEXP kmerScores(SEXP oligos, SEXP ints, SEXP windowSize, SEXP kSize)
{
	int i = 0, j = 0, k = 0, count = 0;
	double *o = REAL(oligos);
	int *mer = INTEGER(ints);
	int wS = asInteger(windowSize);
	double kS = asReal(kSize); // coerce to double
	
	int hS = wS/2; // half the windowSize
	int l = length(ints);
	int *bg = Calloc(length(oligos), int); // rolling distribution of k-mers
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, l + 1));
	double *rans = REAL(ans);
	rans[0] = 0;
	
	// calculate:
	// log(fg/sum(fg)/(1/n)) - log(bg/sum(bg)/(1/n))/kS
	// log(oligos*n) - log(bg*n/count)/kS
	// log(oligos*count/bg)/kS
	
	while (count < wS && i < l) {
		if (mer[i] != NA_INTEGER) {
			count++;
			bg[mer[i] - 1]++;
			
			while (count >= wS) {
				// record score at midpoint
				while (k <= i - hS) {
					if (mer[k] != NA_INTEGER) {
						if (bg[mer[k] - 1] > 0) {
							rans[k + 1] = log(o[mer[k] - 1]*(double)count/(double)bg[mer[k] - 1])/kS;
						} else { // use pseudocount
							rans[k + 1] = log(o[mer[k] - 1]*(double)count)/kS;
						}
					} else {
						rans[k + 1] = 0;
					}
					k++;
				}
				
				// remove from positions before windowSize
				if (mer[j] != NA_INTEGER) {
					count--;
					bg[mer[j] - 1]--;
				}
				j++;
			}
		}
		i++;
	}
	while (k < l) {
		if (mer[k] != NA_INTEGER) {
			if (bg[mer[k] - 1] > 0) {
				rans[k + 1] = log(o[mer[k] - 1]*(double)count/(double)bg[mer[k] - 1])/kS;
			} else { // use pseudocount
				rans[k + 1] = log(o[mer[k] - 1]*(double)count)/kS;
			}
		} else {
			rans[k + 1] = 0;
		}
		k++;
	}
	
	// perform cumulative sum
	for (k = 2; k <= l; k++)
		rans[k] += rans[k - 1];
	
	Free(bg);
	UNPROTECT(1);
	
	return ans;
}

SEXP getHits(SEXP starts, SEXP ends, SEXP left1, SEXP left2, SEXP right1, SEXP right2, SEXP deltaG)
{
	int j, c1 = 0, c2;
	int *pS = INTEGER(starts);
	int *pE = INTEGER(ends);
	int *s1 = INTEGER(left1);
	int *s2 = INTEGER(left2);
	int *e1 = INTEGER(right1);
	int *e2 = INTEGER(right2);
	double *dG = REAL(deltaG);
	int l1 = length(starts);
	int l2 = length(left1);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l2));
	int *rans = INTEGER(ans);
	for (j = 0; j < l2; j++) {
		rans[j] = 0;
		
		while (c1 < l1 && pS[c1] < s1[j])
			c1++;
		
		c2 = c1;
		while (c2 < l1 && pS[c2] <= s2[j]) {
			if (pE[c2] >= e1[j] && pE[c2] <= e2[j]) {
				if (rans[j] > 0) {
					if (dG[c2] < dG[rans[j] - 1])
						rans[j] = c2 + 1;
				} else {
					rans[j] = c2 + 1;
				}
			}
			c2++;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
