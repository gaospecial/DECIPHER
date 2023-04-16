/****************************************************************************
 *                        Manipulating Raw Vectors                          *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// for calloc/free
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

// computes the graph Laplacian from adjacencies between bit vectors
SEXP graphLaplacian(SEXP vecs, SEXP bitMask, SEXP nThreads)
{
	int i, j, k;
	int l = length(vecs); // number of bytes
	int n = 8*l; // number of bits
	int size = n*(n - 1)/2;
	unsigned char *b = RAW(bitMask);
	int nthreads = asInteger(nThreads);
	
	// build a vector of thread-safe pointers
	unsigned char **ptrs = Calloc(l, unsigned char *); // vectors
	for (i = 0; i < l; i++)
		ptrs[i] = RAW(VECTOR_ELT(vecs, i));
	int t = length(VECTOR_ELT(vecs, 0));
	
	// lower triangle of graph laplacian
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, size));
	double *rans;
	rans = REAL(ans);
	
	// count of ones in each bit vector
	SEXP total;
	PROTECT(total = allocVector(INTSXP, n));
	int *tot;
	tot = INTEGER(total);
	
	#pragma omp parallel for private(i,j,k) schedule(guided) num_threads(nthreads)
	for (i = 0; i < n; i++) {
		int byte_i = i/8;
		int bit_i = i % 8;
		unsigned char *v1 = ptrs[byte_i];
		int count = 0;
		
		// count ones
		for (k = 0; k < t; k++)
			if ((v1[k] & b[bit_i]) != 0)
				count++;
		tot[i] = count;
		
		// record indices of ones
		int *ones = (int *) malloc(count*sizeof(int)); // thread-safe on Windows
		int j = 0;
		for (k = 0; j < count; k++)
			if ((v1[k] & b[bit_i]) != 0)
				ones[j++] = k;
		
		for (j = i + 1; j < n; j++) {
			int byte_j = j/8;
			int bit_j = j % 8;
			unsigned char *v2 = ptrs[byte_j];
			
			int index = n*i - i*(i + 1)/2 + j - i - 1;
			rans[index] = 0;
			
			// count cooccurrences of shared ones
			for (k = 0; k < count; k++)
				if ((v2[ones[k]] & b[bit_j]) != 0)
					rans[index]++;
		}
		free(ones);
	}
	Free(ptrs);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret_list, 0, ans);
	SET_VECTOR_ELT(ret_list, 1, total);
	
	UNPROTECT(3);
	
	return ret_list;
}

// assigns clusters by voting on clusters
SEXP assignClusters(SEXP vecs, SEXP bitMask, SEXP clusters, SEXP nThreads)
{
	int i, j, k;
	int l = length(vecs); // number of bytes
	int n = 8*l; // number of bits
	unsigned char *b = RAW(bitMask);
	int *c = INTEGER(clusters);
	int nthreads = asInteger(nThreads);
	
	// build a vector of thread-safe pointers
	unsigned char **ptrs = Calloc(l, unsigned char *); // vectors
	for (i = 0; i < l; i++)
		ptrs[i] = RAW(VECTOR_ELT(vecs, i));
	int t = length(VECTOR_ELT(vecs, 0));
	
	// record the maximum cluster number
	int m = 0;
	for (i = 0; i < n; i++)
		if (c[i] > m)
			m = c[i];
	m++;
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, t));
	int *rans;
	rans = INTEGER(ans);
	
	#pragma omp parallel num_threads(nthreads)
	{
		int *votes = (int *) malloc(m*sizeof(int)); // thread-safe on Windows
		votes[0] = 0; // assume clusters start from 1
		
		#pragma omp for private(i,j,k) schedule(guided)
		for (k = 0; k < t; k++) { // index in bit vectors
			// reset votes
			for (i = 1; i < m; i++)
				votes[i] = 0;
			
			// record votes
			for (i = 0; i < n; i++) {
				int byte_i = i/8;
				int bit_i = i % 8;
				unsigned char *v = ptrs[byte_i];
				if ((v[k] & b[bit_i]) != 0)
					votes[c[i]]++;
			}
			
			// plurality vote
			j = 0; // default vote
			for (i = 1; i < m; i++)
				if (votes[i] > votes[j])
					j = i; // new vote
			rans[k] = j;
		}
		
		free(votes);
	}
	Free(ptrs);
	
	UNPROTECT(1);
	
	return ans;
}
