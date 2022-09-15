/****************************************************************************
 *                       Cluster by Spherical K-means                       *
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

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// for calloc/free
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

SEXP sphericalKmeans(SEXP X, SEXP Y, SEXP Kmeans, SEXP tolerance, SEXP maxIterations, SEXP verbose, SEXP nThreads)
{	
	int i, j, k, it, count, org;
	double temp, best, score;
	double *b, *c; // previous and current centers
	
	SEXP dim = getAttrib(X, R_DimSymbol);
	int d = INTEGER(dim)[0];
	int n = INTEGER(dim)[1];
	double *x = REAL(X);
	int *y = INTEGER(Y);
	int K = asInteger(Kmeans);
	double tol = 1 - asReal(tolerance);
	int maxIt = asInteger(maxIterations);
	int v = asInteger(verbose);
	int nthreads = asInteger(nThreads);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, n));
	int *rans = INTEGER(ans);
	for (i = 0; i < n; i++)
		rans[i] = y[i] - 1; // start from 0
	
	// initialize upper and lower bounds (Spherical Simplified Elkan’s Algorithm)
	double *l = (double *) calloc(n, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *u = (double *) malloc(K*n*sizeof(double)); // thread-safe on Windows
	for (i = 0; i < K*n; i++)
		u[i] = 1; // initialize upper bound to one
	
	it = 0; // iteration
	do {
		it++;
		// determine cluster centroids
		c = (double *) calloc(d*K, sizeof(double)); // initialized to zero (thread-safe on Windows)
		for (j = 0; j < n; j++)
			for (i = 0; i < d; i++)
				c[d*rans[j] + i] += x[d*j + i];
		
		// calculate row and column means
		double *colsums = (double *) calloc(K, sizeof(double)); // initialized to zero (thread-safe on Windows)
		double *rowsums = (double *) calloc(d, sizeof(double)); // initialized to zero (thread-safe on Windows)
		for (j = 0; j < K; j++) {
			for (i = 0; i < d; i++) {
				colsums[j] += c[d*j + i];
				rowsums[i] += c[d*j + i];
			}
		}
		
		// relocate empty clusters
		double *sumsquares = (double *) calloc(K, sizeof(double)); // initialized to zero (thread-safe on Windows)
		for (j = 0; j < K; j++) {
			if (colsums[j] == 0) { // empty cluster
				// reassign to rowmeans
				for (i = 0; i < d; i++)
					c[d*j + i] = rowsums[i]/(double)K;
				// reset upper bound
				for (i = 0; i < n; i++)
					u[K*i + j] = 1;
			}
			for (i = 0; i < d; i++)
				sumsquares[j] += c[d*j + i]*c[d*j + i];
			sumsquares[j] = sqrt(sumsquares[j]);
		}
		free(rowsums);
		free(colsums);
		
		// normalize centroids to unit length
		for (j = 0; j < K; j++)
			for (i = 0; i < d; i++)
				c[d*j + i] /= sumsquares[j];
		free(sumsquares);
		
		if (it > 1) {
			double *p = (double *) calloc(K, sizeof(double)); // initialized to zero (thread-safe on Windows)
			int less = 1;
			for (j = 0; j < K; j++) {
				for (i = 0; i < d; i++)
					p[j] += c[d*j + i]*b[d*j + i]; // dot product with prior center
				if (p[j] <= tol)
					less = 0;
			}
			free(b);
			if (less) { // all within tolerance
				free(p);
				break; // centers stopped moving
			}
			
			for (j = 0; j < n; j++) {
				// adjust lower bound
				temp = (1 - l[j]*l[j])*(1 - p[rans[j]]*p[rans[j]]);
				l[j] *= p[rans[j]];
				if (temp > 0)
					l[j] -= sqrt(temp);
				
				// adjust upper bounds
				for (i = 0; i < K; i++) {
					temp = (1 - u[j*K + i]*u[j*K + i])*(1 - p[i]*p[i]);
					u[j*K + i] *= p[i];
					if (temp > 0)
						u[j*K + i] += sqrt(temp);
				}
			}
			free(p);
		}
		
		count = 0;
		#pragma omp parallel for private(i,j,k,org,best,score) num_threads(nthreads)
		for (j = 0; j < n; j++) {
			best = -1e50;
			org = rans[j];
			for (k = 0; k < K; k++) {
				if (u[j*K + k] > l[j]) {
					score = 0;
					for (i = 0; i < d; i++)
						score += c[d*k + i]*x[d*j + i];
					u[j*K + k] = score;
					if (score > best) {
						best = score;
						rans[j] = k;
						l[j] = score;
					}
				}
			}
			if (rans[j] != org)
				count++;
		}
		
		if (count == 0) // no changes
			break;
		b = c;
		if (v == 1) {
			Rprintf("\riteration %d (%1.1f%% stability)  ", it, 100*(1 - (double)count/(double)n));
		} else {
			R_CheckUserInterrupt();
		}
	} while (it < maxIt);
	free(c);
	free(l);
	free(u);
	
	for (i = 0; i < n; i++)
		rans[i]++; // start from 1
	
	if (v != 0)
		Rprintf("\riteration %d (100%% stability)  \n\n", it);
	
	UNPROTECT(1);
	
	return ans;
}
