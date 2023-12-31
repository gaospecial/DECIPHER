/****************************************************************************
 *                       Converts Sequence To Numbers                       *
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

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

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

// DECIPHER header file
#include "DECIPHER.h"

static void alphabetFrequency(const Chars_holder *P, int *bits, int position)
{
	const char *p;
	p = (P->ptr + position);
	
	switch (*p) {
		case 1: // A
			*(bits) = 0;
			break;
		case 2: // C
			*(bits) = 1;
			break;
		case 4: // G
			*(bits) = 2;
			break;
		case 8: // T
			*(bits) = 3;
			break;
		default: // other
			*(bits) = -1;
			break;
	}
}

// changes repeat regions to NAs
SEXP maskRepeats(SEXP e, int n, int l1, int l2, int l3)
{
	if (MAYBE_SHARED(e))
		error(".Call function 'maskRepeats' called in incorrect context.");
	
	int i, p, j, k;
	int *x = INTEGER(e); // enumerated sequence
	int l = length(e);
//	int n = asInteger(size);
//	int l1 = asInteger(minL); // min period
//	int l2 = asInteger(maxL); // max period
//	int l3 = asInteger(totL); // min length of repeat
	
	i = 0; // current position
	while (i < (l - l2)) {
		if (x[i]!=NA_INTEGER) {
			for (p = l1; p <= l2; p++) { // periodicity
				if (x[i]==x[i + p]) { // repeat
					j = i + 1;
					
					while (j < (l - p)) {
						if (x[j]!=x[j + p])
							break;
						j++;
					}
					
					if ((j - i + n) > p && // continuous repeat
						(j + p - i + n) > l3) {
						for (k = i; k <= (j + p - 1); k++)
							x[k] = NA_INTEGER;
						i = k - 1;
						break;
					}
				}
			}
		}
		i++;
	}
	
	return R_NilValue;
}

//ans_start <- .Call("enumerateSequence", myDNAStringSet, wordSize, FALSE, PACKAGE="DECIPHER")
SEXP enumerateSequence(SEXP x, SEXP wordSize, SEXP mask)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, wS, maskReps, sum, ambiguous, *rans;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	wS = asInteger(wordSize); // [1 to 15]
	maskReps = asInteger(mask);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, x_length));
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = 4;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*4;
	}
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		SEXP ans;
		if ((x_i.length - wS + 1) < 1) {
			PROTECT(ans = allocVector(INTSXP, 0));
		} else {
			PROTECT(ans = allocVector(INTSXP, x_i.length - wS + 1));
			rans = INTEGER(ans);
			int bases[wS];
			for (j = 0; j < (wS - 1); j++) {
				alphabetFrequency(&x_i, &bases[j], j); // fill initial numbers
			}
			for (j = wS - 1; j < x_i.length; j++) {
				alphabetFrequency(&x_i, &bases[wS - 1], j);
				sum = bases[0];
				ambiguous = 0;
				if (bases[0] < 0)
					ambiguous = 1;
				for (k = 1; k < wS; k++) {
					sum += bases[k]*pwv[k];
					if (bases[k] < 0)
						ambiguous = 1;
					bases[k - 1] = bases[k]; // shift numbers left
				}
				if (ambiguous) {
					*(rans + j - wS + 1) = NA_INTEGER;
				} else {
					*(rans + j - wS + 1) = sum;
				}
			}
		}
		
		if (maskReps)
			maskRepeats(ans, wS, 7, 12, 30);
		
		SET_VECTOR_ELT(ret_list, i, ans);
		UNPROTECT(1);
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ret_list;
}

static void alphabetFrequencyAA(const Chars_holder *P, int *bits, int position)
{
	const char *p;
	p = (P->ptr + position);
	
	switch (*p) {
		case 65: // A
			*(bits) = 0;
			break;
		case 82: // R
			*(bits) = 1;
			break;
		case 78: // N
			*(bits) = 2;
			break;
		case 68: // D
			*(bits) = 3;
			break;
		case 67: // C
			*(bits) = 4;
			break;
		case 81: // Q
			*(bits) = 5;
			break;
		case 69: // E
			*(bits) = 6;
			break;
		case 71: // G
			*(bits) = 7;
			break;
		case 72: // H
			*(bits) = 8;
			break;
		case 73: // I
			*(bits) = 9;
			break;
		case 76: // L
			*(bits) = 10;
			break;
		case 75: // K
			*(bits) = 11;
			break;
		case 77: // M
			*(bits) = 12;
			break;
		case 70: // F
			*(bits) = 13;
			break;
		case 80: // P
			*(bits) = 14;
			break;
		case 83: // S
			*(bits) = 15;
			break;
		case 84: // T
			*(bits) = 16;
			break;
		case 87: // W
			*(bits) = 17;
			break;
		case 89: // Y
			*(bits) = 18;
			break;
		case 86: // V
			*(bits) = 19;
			break;
		default: // other
			*(bits) = -1;
			break;
	}
}

//ans_start <- .Call("enumerateSequenceAA", myDNAStringSet, wordSize, PACKAGE="DECIPHER")
SEXP enumerateSequenceAA(SEXP x, SEXP wordSize)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, wS, sum, ambiguous, *rans;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	wS = asInteger(wordSize); // [1 to 7]
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, x_length));
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = 20;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*20;
	}
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		SEXP ans;
		if ((x_i.length - wS + 1) < 1) {
			PROTECT(ans = allocVector(INTSXP, 0));
		} else {
			PROTECT(ans = allocVector(INTSXP, x_i.length - wS + 1));
			rans = INTEGER(ans);
			int bases[wS];
			for (j = 0; j < (wS - 1); j++) {
				alphabetFrequencyAA(&x_i, &bases[j], j); // fill initial numbers
			}
			for (j = wS - 1; j < x_i.length; j++) {
				alphabetFrequencyAA(&x_i, &bases[wS - 1], j);
				sum = bases[0];
				ambiguous = 0;
				if (bases[0] < 0)
					ambiguous = 1;
				for (k = 1; k < wS; k++) {
					sum += bases[k]*pwv[k];
					if (bases[k] < 0)
						ambiguous = 1;
					bases[k - 1] = bases[k]; // shift numbers left
				}
				if (ambiguous) {
					*(rans + j - wS + 1) = NA_INTEGER;
				} else {
					*(rans + j - wS + 1) = sum;
				}
			}
		}
		
		SET_VECTOR_ELT(ret_list, i, ans);
		UNPROTECT(1);
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ret_list;
}

int pop(unsigned int x)
{
	x = x - ((x >> 1) & 0x55555555);
	x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
	x = x + (x >> 4) & 0xF0F0F0F;
	return((x*0x1010101) >> 24);
}

//ans_start <- .Call("enumerateGappedSequence", myDNAStringSet, wordSize, ordering, PACKAGE="DECIPHER")
SEXP enumerateGappedSequence(SEXP x, SEXP wordSize, SEXP ordering)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int i, j, k, wS, sum, ambiguous, *rans, *p;
	int *o = INTEGER(ordering);
	int l = length(ordering);
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	wS = asInteger(wordSize); // [1 to 15]
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, l*2));
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = 4;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*4;
	}
	
	for (i = 0; i < l; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, *(o + i) - 1);
		SEXP ans, pos;
		if ((x_i.length - wS + 1) < 1) {
			PROTECT(ans = allocVector(INTSXP, 0));
			PROTECT(pos = allocVector(INTSXP, 0));
		} else {
			int *POS = Calloc(x_i.length - wS + 1, int); // initialized to zero
			int *ANS = Calloc(x_i.length - wS + 1, int); // initialized to zero
			int bases[wS];
			int count = 0;
			int repeat[66] = {-1}; // array of previous occurrences
			int values[66];
			int bitCount, temp;
			for (j = 0; j < x_i.length; j++) {
				if (count < wS - 1) {
					if (x_i.ptr[j] != 16 && x_i.ptr[j] != 64) { // not gap
						alphabetFrequency(&x_i, &bases[count], j); // fill initial numbers
						count++;
					}
					continue;
				}
				if (x_i.ptr[j] == 16 || x_i.ptr[j] == 64) // gap
					continue;
				alphabetFrequency(&x_i, &bases[wS - 1], j);
				count++;
				*(POS + count - wS) = j + 1;
				sum = bases[0];
				ambiguous = 0;
				if (bases[0] < 0)
					ambiguous = 1;
				for (k = 1; k < wS; k++) {
					sum += bases[k]*pwv[k];
					if (bases[k] < 0)
						ambiguous = 1;
					bases[k - 1] = bases[k]; // shift numbers left
				}
				if (ambiguous) {
					*(ANS + count - wS) = -1;
				} else {
					*(ANS + count - wS) = sum;
				}
			}
			
			if ((count - wS + 1) < 1) {
				PROTECT(ans = allocVector(INTSXP, 0));
				PROTECT(pos = allocVector(INTSXP, 0));
			} else {
				PROTECT(ans = allocVector(INTSXP, count - wS + 1));
				PROTECT(pos = allocVector(INTSXP, count - wS + 1));
				rans = INTEGER(ans);
				p = INTEGER(pos);
				for (j = 0; j <= count - wS; j++) {
					if (*(ANS + j) == -1) {
						*(rans + j) = NA_INTEGER;
					} else {
						temp = *(ANS + j);
						
						// determine the k-mer's signature
						bitCount = pop((unsigned int)temp);
						if (temp < 0)
							bitCount += 33;
						
						// set all recent repeats to NA
						if (repeat[bitCount] < 0) {
							// initialize repeat
							repeat[bitCount] = j;
							values[bitCount] = temp;
							*(rans + j) = temp;
						} else { // check for a repeat
							if (values[bitCount]==temp) {
								*(rans + j) = NA_INTEGER;
								*(rans + repeat[bitCount]) = NA_INTEGER;
							} else { // not a recent repeat
								repeat[bitCount] = j;
								values[bitCount] = temp;
								*(rans + j) = temp;
							}
						}
					}
					*(p + j) = *(POS + j);
				}
			}
			Free(POS);
			Free(ANS);
		}
		
		SET_VECTOR_ELT(ret_list, 2*i, ans);
		SET_VECTOR_ELT(ret_list, 2*i + 1, pos);
		UNPROTECT(2);
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ret_list;
}

//ans_start <- .Call("enumerateGappedSequenceAA", myDNAStringSet, wordSize, ordering, PACKAGE="DECIPHER")
SEXP enumerateGappedSequenceAA(SEXP x, SEXP wordSize, SEXP ordering)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int i, j, k, wS, sum, ambiguous, *rans, *p;
	int *o = INTEGER(ordering);
	int l = length(ordering);
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	wS = asInteger(wordSize); // [1 to 7]
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, l*2));
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = 20;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*20;
	}
	
	for (i = 0; i < l; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, *(o + i) - 1);
		SEXP ans, pos;
		if ((x_i.length - wS + 1) < 1) {
			PROTECT(ans = allocVector(INTSXP, 0));
			PROTECT(pos = allocVector(INTSXP, 0));
		} else {
			int *POS = Calloc(x_i.length - wS + 1, int); // initialized to zero
			int *ANS = Calloc(x_i.length - wS + 1, int); // initialized to zero
			int bases[wS];
			int count = 0;
			int repeat[66] = {-1}; // array of previous occurrences
			int values[66];
			int bitCount, temp;
			for (j = 0; j < x_i.length; j++) {
				if (count < wS - 1) {
					if (x_i.ptr[j] != 45 && x_i.ptr[j] != 46) { // not gap
						alphabetFrequencyAA(&x_i, &bases[count], j); // fill initial numbers
						count++;
					}
					continue;
				}
				if (x_i.ptr[j] == 45 || x_i.ptr[j] == 46) // gap
					continue;
				alphabetFrequencyAA(&x_i, &bases[wS - 1], j);
				count++;
				*(POS + count - wS) = j + 1;
				sum = bases[0];
				ambiguous = 0;
				if (bases[0] < 0)
					ambiguous = 1;
				for (k = 1; k < wS; k++) {
					sum += bases[k]*pwv[k];
					if (bases[k] < 0)
						ambiguous = 1;
					bases[k - 1] = bases[k]; // shift numbers left
				}
				if (ambiguous) {
					*(ANS + count - wS) = -1;
				} else {
					*(ANS + count - wS) = sum;
				}
			}
			
			if ((count - wS + 1) < 1) {
				PROTECT(ans = allocVector(INTSXP, 0));
				PROTECT(pos = allocVector(INTSXP, 0));
			} else {
				PROTECT(ans = allocVector(INTSXP, count - wS + 1));
				PROTECT(pos = allocVector(INTSXP, count - wS + 1));
				rans = INTEGER(ans);
				p = INTEGER(pos);
				for (j = 0; j <= count - wS; j++) {
					if (*(ANS + j) == -1) {
						*(rans + j) = NA_INTEGER;
					} else {
						temp = *(ANS + j);
						
						// determine the k-mer's signature
						bitCount = pop((unsigned int)temp);
						if (temp < 0)
							bitCount += 33;
						
						// set all recent repeats to NA
						if (repeat[bitCount] < 0) {
							// initialize repeat
							repeat[bitCount] = j;
							values[bitCount] = temp;
							*(rans + j) = temp;
						} else { // check for a repeat
							if (values[bitCount]==temp) {
								*(rans + j) = NA_INTEGER;
								*(rans + repeat[bitCount]) = NA_INTEGER;
							} else { // not a recent repeat
								repeat[bitCount] = j;
								values[bitCount] = temp;
								*(rans + j) = temp;
							}
						}
					}
					*(p + j) = *(POS + j);
				}
			}
			Free(POS);
			Free(ANS);
		}
		
		SET_VECTOR_ELT(ret_list, 2*i, ans);
		SET_VECTOR_ELT(ret_list, 2*i + 1, pos);
		UNPROTECT(2);
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ret_list;
}

static void alphabetFrequencyReducedAA(const Chars_holder *P, int *bits, int position, int *alpha)
{
	const char *p;
	p = (P->ptr + position);
	
	switch (*p) {
		case 65: // A
			*(bits) = *(alpha);
			break;
		case 82: // R
			*(bits) = *(alpha + 1);
			break;
		case 78: // N
			*(bits) = *(alpha + 2);
			break;
		case 68: // D
			*(bits) = *(alpha + 3);
			break;
		case 67: // C
			*(bits) = *(alpha + 4);
			break;
		case 81: // Q
			*(bits) = *(alpha + 5);
			break;
		case 69: // E
			*(bits) = *(alpha + 6);
			break;
		case 71: // G
			*(bits) = *(alpha + 7);
			break;
		case 72: // H
			*(bits) = *(alpha + 8);
			break;
		case 73: // I
			*(bits) = *(alpha + 9);
			break;
		case 76: // L
			*(bits) = *(alpha + 10);
			break;
		case 75: // K
			*(bits) = *(alpha + 11);
			break;
		case 77: // M
			*(bits) = *(alpha + 12);
			break;
		case 70: // F
			*(bits) = *(alpha + 13);
			break;
		case 80: // P
			*(bits) = *(alpha + 14);
			break;
		case 83: // S
			*(bits) = *(alpha + 15);
			break;
		case 84: // T
			*(bits) = *(alpha + 16);
			break;
		case 87: // W
			*(bits) = *(alpha + 17);
			break;
		case 89: // Y
			*(bits) = *(alpha + 18);
			break;
		case 86: // V
			*(bits) = *(alpha + 19);
			break;
		default: // other
			*(bits) = -1;
			break;
	}
}

SEXP enumerateSequenceReducedAA(SEXP x, SEXP wordSize, SEXP alphabet, SEXP mask)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, wS, maskReps, sum, ambiguous, *rans;
	int *alpha = INTEGER(alphabet);
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	wS = asInteger(wordSize); // [1 to 20]
	maskReps = asInteger(mask);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, x_length));
	
	int m = 0; // hold max of alphabet
	for (i = 0; i < 20; i++) {
		if (*(alpha + i) > m)
			m = *(alpha + i);
	}
	m++; // start at 1
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = m;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*m;
	}
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		SEXP ans;
		if ((x_i.length - wS + 1) < 1) {
			PROTECT(ans = allocVector(INTSXP, 0));
		} else {
			PROTECT(ans = allocVector(INTSXP, x_i.length - wS + 1));
			rans = INTEGER(ans);
			int bases[wS];
			for (j = 0; j < (wS - 1); j++) {
				alphabetFrequencyReducedAA(&x_i, &bases[j], j, alpha); // fill initial numbers
			}
			for (j = wS - 1; j < x_i.length; j++) {
				alphabetFrequencyReducedAA(&x_i, &bases[wS - 1], j, alpha);
				sum = bases[0];
				ambiguous = 0;
				if (bases[0] < 0)
					ambiguous = 1;
				for (k = 1; k < wS; k++) {
					sum += bases[k]*pwv[k];
					if (bases[k] < 0)
						ambiguous = 1;
					bases[k - 1] = bases[k]; // shift numbers left
				}
				if (ambiguous) {
					*(rans + j - wS + 1) = NA_INTEGER;
				} else {
					*(rans + j - wS + 1) = sum;
				}
			}
		}
		
		if (maskReps)
			maskRepeats(ans, wS, 3, 11, 15);
		
		SET_VECTOR_ELT(ret_list, i, ans);
		UNPROTECT(1);
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ret_list;
}

// returns the size of a balanced alphabet with equivalent entropy
SEXP alphabetSize(SEXP x)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, letter;
	double sum = 0;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 1));
	double *rans = REAL(ans);
	rans[0] = 0;
	
	double dist[4] = {0}; // distribution of nucleotides
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (j = 0; j < x_i.length; j++) {
			alphabetFrequency(&x_i, &letter, j);
			if (letter >= 0)
				dist[letter]++;
		}
	}
	
	for (i = 0; i < 4; i++)
		sum += dist[i];
	
	double p; // proportion of each letter
	for (i = 0; i < 4; i++) {
		p = dist[i]/sum;
		if (p > 0)
			rans[0] -= p*log(p); // negative entropy
	}
	rans[0] = exp(rans[0]);
	
	UNPROTECT(1);
	
	return ans;
}

// returns the size of a balanced alphabet with equivalent entropy
SEXP alphabetSizeReducedAA(SEXP x, SEXP alphabet)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, letter;
	int *alpha = INTEGER(alphabet);
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 1));
	double *rans = REAL(ans);
	rans[0] = 0;
	
	int m = 0; // hold max of alphabet
	for (i = 0; i < 20; i++) {
		if (*(alpha + i) > m)
			m = *(alpha + i);
	}
	m++; // start at 1
	
	double dist[m]; // distribution
	for (i = 0; i < m; i++)
		dist[i] = 0;
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (j = 0; j < x_i.length; j++) {
			alphabetFrequencyReducedAA(&x_i, &letter, j, alpha);
			if (letter >= 0)
				dist[letter]++;
		}
	}
	
	double sum = 0;
	for (i = 0; i < m; i++)
		sum += dist[i];
	
	double p; // proportion of each letter
	for (i = 0; i < m; i++) {
		p = dist[i]/sum;
		if (p > 0)
			rans[0] -= p*log(p); // negative entropy
	}
	rans[0] = exp(rans[0]);
	
	UNPROTECT(1);
	
	return ans;
}

// equalivent to t(oligonucleotideFrequency(x, wordSize, fast.moving.side="left"))
SEXP frequencies(SEXP x, SEXP wordSize)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, wS, sum, ambiguous, *rans;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	wS = asInteger(wordSize); // [1 to 15]
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = 4;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*4;
	}
	
	int words = 4;
	for (i = 1; i < wS; i++)
		words *= 4;
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, words, x_length));
	rans = INTEGER(ans);
	for (i = 0; i < words*x_length; i++)
		rans[i] = 0;
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		int bases[wS];
		for (j = 0; j < (wS - 1); j++) {
			alphabetFrequency(&x_i, &bases[j], j); // fill initial numbers
		}
		for (j = wS - 1; j < x_i.length; j++) {
			alphabetFrequency(&x_i, &bases[wS - 1], j);
			sum = bases[0];
			ambiguous = 0;
			if (bases[0] < 0)
				ambiguous = 1;
			for (k = 1; k < wS; k++) {
				sum += bases[k]*pwv[k];
				if (bases[k] < 0)
					ambiguous = 1;
				bases[k - 1] = bases[k]; // shift numbers left
			}
			if (ambiguous == 0)
				(*(rans + words*i + sum))++;
		}
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP frequenciesReducedAA(SEXP x, SEXP wordSize, SEXP alphabet)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, wS, sum, ambiguous, *rans;
	int *alpha = INTEGER(alphabet);
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	wS = asInteger(wordSize); // [1 to 20]
	
	int m = 0; // hold max of alphabet
	for (i = 0; i < 20; i++) {
		if (*(alpha + i) > m)
			m = *(alpha + i);
	}
	m++; // start at 1
	
	// fill the position weight vector
	int pwv[wS]; // wS[0] is ignored
	if (wS > 1)
		pwv[1] = m;
	for (i = 2; i < wS; i++) {
		pwv[i] = pwv[i - 1]*m;
	}
	
	int words = m;
	for (i = 1; i < wS; i++)
		words *= m;
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, words, x_length));
	rans = INTEGER(ans);
	for (i = 0; i < words*x_length; i++)
		rans[i] = 0;
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		int bases[wS];
		for (j = 0; j < (wS - 1); j++) {
			alphabetFrequencyReducedAA(&x_i, &bases[j], j, alpha); // fill initial numbers
		}
		for (j = wS - 1; j < x_i.length; j++) {
			alphabetFrequencyReducedAA(&x_i, &bases[wS - 1], j, alpha);
			sum = bases[0];
			ambiguous = 0;
			if (bases[0] < 0)
				ambiguous = 1;
			for (k = 1; k < wS; k++) {
				sum += bases[k]*pwv[k];
				if (bases[k] < 0)
					ambiguous = 1;
				bases[k - 1] = bases[k]; // shift numbers left
			}
			if (ambiguous == 0)
				(*(rans + words*i + sum))++;
		}
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ans;
}
