/************************************************************************************************
 *																								*
 * Compute simultaneous tolerance bounds for a matrix of random vectors representing the null	*
 * distribution, e.g. N(0,1). 																	*
 *			   																			        *
 * Author:			Dr. André Schützenmeister											        *
 *					Centralised and Point of Care Solutions										*
 *																								*
 *					Roche Diagnostics GmbH														*
 *					DXREBC..6164																*
 *					Nonnenwald 2																*
 *					82377 Penzberg / Germany													*
 *																								*																										*
 *				  																			    *
 * Last modified:	2016-06-07															        *
 *																								*
 * - added OpenMP based parallelization to coverage-computation and STB-computation             *	
 * - removed all unnecessary content before putting it on CRAN          			            *
 *																								*
 ************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <omp.h>

#ifndef PI
#define PI 3.141592653589793115998
#endif


#ifndef EPS
#define EPS 1.e-20		/* precision to be used, e.g. in calls to SASquantile */
#endif


int compare_doubles(const void *X, const void *Y);
double SASquantile(double *v, double p, int len, double tol);
double coverage(double *mat, double *lower, double *upper, int nCol, int nRow, int nCpu);
double abs_double(double x);
void getSTB(double *mat, int *nCol, int *nRow, double *alpha, double *tol, int *max_iter, int *nCpu, double *Q, double *cov);
void getSASquantile(double *v, double *p, int *len, int *numQ, double *tol, double *Q);
void getCoverage(double *mat, double *lower, double *upper, int *nCol, int *nRow, int *nCpu, double *cov);




/*
	Interface for SAS PCTLDEF5 quantile computationen from R.
	
	v		(double *) pointer to the 1st element of the double vector from which 1 or multiple quantile should be computed
	p		(double *) pointer to the 1st element of a double vector containing requestes percentages for quantiles
	len		(int *) pointer to an integer specifying the length of 'v'
	numQ	(int *) pointer to an integer specifying the number of requested quantiles
	Q		(double *) pointer to the 1st element of a double vector to which quantiles will be written
*/

void getSASquantile(double *v, double *p, int *len, int *numQ, double *tol, double *Q)
{
	int i;

	if(*numQ <= 0)
		return;
	
	for(i=0;i<*numQ;i++)
	{
		Q[i] = SASquantile(v, p[i], *len, *tol);
	}
	
	return;
}



/*
	Compute coverage of the 100(1-alpha)% STB accessible from R.
*/

void getCoverage(double *mat, double *lower, double *upper, int *nCol, int *nRow, int *nCpu, double *cov)
{
	cov[0] = coverage(mat, lower, upper, *nCol, *nRow, *nCpu);
	return;
}


int compare_doubles(const void *X, const void *Y)
{
  double x = *((double *)X);
  double y = *((double *)Y);

  if (x > y) return 1;
  else if (x < y) return -1;
  else return 0;
}

double SASquantile(double *v, double p, int len, double tol)
{
	double np, g;
	int j;
	
	qsort(v, len, sizeof(double), compare_doubles);
	np = (len) * p;
	j = floor(np);
	g = np - j;
	if(g < tol)
	{
		return( (v[j-1]+v[j])/2 );		/* indices run from 0...n-1 for vector of length n */
	}
	else
	{
		return( v[j] );
	}			
}



/*
	Compute coverage, i.e. the number of rows completely enclosed by upper bounds 'upper'
	and lower bound 'lower' divided by the number of rows.
*/

double coverage(double *mat, double *lower, double *upper, int nCol, int nRow, int nCpu)
{
	int covered=0, check=0, i=0, j=0;
	
	omp_set_num_threads(nCpu);									/* OpenMP support, set number of CPUs */
	
	#pragma omp parallel for \
			private(i, j) \
			firstprivate(check, lower, upper, nCol, nRow)\
			shared(mat, covered ) 								/* OpenMP parallelization */
	
	for(i=0; i<nRow; i++)										/* over each simulated value of the j-th row */						
	{
		check = 0;
		
		for(j=0; j<nCol; j++)									/* over columns --> j-th order statistics */				
		{
			if(mat[j*nRow+i] < lower[j])
			{
				check = 1;
				break;
			}
			if(mat[j*nRow+i] > upper[j])
			{
				check = 1;
				break;
			}
		}
		if(check == 0) 
		{
			#pragma omp atomic								/* protect access of shared variable (memory) */
			covered += 1;
		} 
	}
	return((double)covered / nRow);
}

/* Absolute value of a double-variable. */

double abs_double(double x)
{
	if(x < 0) return(x * (-1));
	else return(x);
}

void sort(double * vec, int * len)
{
	qsort(vec, *len, sizeof(double), compare_doubles);
}



/* 
	Compute bounds of the 100(1-alpha)% simultaneous tolerance band for matrix 'mat'. 
	
	mat				... (*double) array, representing the matrix of sorted random vectors								
	nCol   			... (*integer) pointer to the variable which contains the number of columns of 'mat'
	nRow			... (*integer) pointer to the variable which contains thenumber of rown of 'mat'
	alpha			... (*double) pointer to the variable of the requested 1-alpha tolerance level of the STB
	tol				... (*double) pointer to the variable of convergence tolerance for the bisection algorithm
	max_iter        ... (*integer) pointer to the variable of the maximum number of iterations of the bisection algorithm
	nCpu			... (*integer) pointer to the variable specifyin the number of cores to be used for parallel processing (currently not used)

*/

void getSTB(double *mat, int *nCol, int *nRow, double *alpha, double *tol, int *max_iter, int *nCpu, double *Q, double *cov)
{
	int check=1, iter=0, i, j;
	double include, alpha_old, tmp, best_cov; // best_alpha;
	double *tmp_col, *lower, *upper, *best_Q;
		
	lower = calloc(*nCol, sizeof(double));
	upper = calloc(*nCol, sizeof(double));
	best_Q = calloc(2*(*nCol), sizeof(double));

	include = 1-(*alpha);											/* requested simultaneous tolerance limit, i.e. coverage */
	alpha_old = *alpha;
	*alpha = (*alpha)/2;
	best_cov=1.0;
	
	omp_set_num_threads(*nCpu);										/* set number of cores */

	while(check==1)
	{
		iter += 1;
		
		#pragma omp parallel for\
			shared(mat, lower, upper, Q)\
			private(i,j,tmp_col)
		for(i=0; i<*nCol; i++)																/* pont-wise tolerance limits */
		{			
			tmp_col = calloc(*nRow, sizeof(double));
		
			for(j=0;j<(*nRow);j++)
			{
				tmp_col[j] = mat[i*(*nRow)+j];										/* copy i-th col of the matrix 'mat'*/
			}
			Q[i*2]=SASquantile(tmp_col, *alpha, *nRow, EPS);					/* 1st row */
			lower[i]=Q[i*2];
			Q[i*2+1]=SASquantile(tmp_col, 1-(*alpha), *nRow, EPS);		/* 2nd row */
			upper[i]=Q[i*2+1];
			
			free(tmp_col);
		}	
		*cov=coverage(mat, lower, upper, *nCol, *nRow, *nCpu);				/* compute coverage of current STB */
		if(*cov >= include)																		/* at least 100(1-alpha)% coverage */
		{
			if(*cov < best_cov)																	/* coverage better than former result */
			{
				//best_alpha = *alpha;
				best_Q = Q;
				best_cov = *cov;
			}
		}
		if( abs_double((*alpha) - alpha_old)/2 == 0 )					/* terminate if alpha cannot become smaller */
		{ 
			check = 0;
		}
		if( iter == *max_iter ) 															/* terminate if max number of iterations reached */
		{
			check = 0;
		}
				
		if( abs_double(*cov - include) <= *tol &&  (*cov - include) >= 0 )
		{
			check = 0;																					/* convergence tolerance reached */
		}
		else																									/* bisection step */
		{
			if( (*cov - include) < 0)
			{
				tmp = *alpha;
				*alpha = *alpha - abs_double(alpha_old-(*alpha))/2;
				alpha_old = tmp;
			}
			else
			{
				tmp = *alpha;
				*alpha = *alpha +  abs_double(alpha_old-(*alpha))/2;
				alpha_old = tmp;
			}
		}
	}
	*cov = best_cov;			/* Use best values from iteration history */
	Q = best_Q;
}

