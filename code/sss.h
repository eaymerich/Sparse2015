/**********************************************************
* Sparse 2015 Project
* Directed by: Euripides Montagne.
*              eurip@cs.ucf.edu
*
* Coded by: Edward Aymerich. 
*           edward.aymerich@knights.ucf.edu
*
* Copyright 2015 - University of Central Florida.
**********************************************************/

#ifndef SP_SSS
#define SP_SSS

#include <stdio.h>
#include <stdlib.h>
#include "coo.h"

/**
 * Sparse matrix structure stored
 * in Symmetric Sparse Skyline (SSS) format.
 * Only the lower triangle is stored.
 */
typedef struct {
	double *dval;		 // Values in main diagonal.
	double *val;       // Value of each non-zero element.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int *ptr; // Initial position of each row.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} SSS_Matrix;

void sss_mvp_sym(SSS_Matrix *m, double *x, double *y);

void sss_init_matrix(SSS_Matrix *m);
void sss_free_matrix(SSS_Matrix *m);
void sss_print_matrix(SSS_Matrix *m);

void coo2sss(COO_Matrix *in, SSS_Matrix *out);
int sss_load_matrix(char* filename, SSS_Matrix *m);

/**
 * Stores on vector y the product of matrix A and vector x:
 *   y = Ax
 */
void sss_mvp_sym(SSS_Matrix *m, double *x, double *y){
	int r,j,c;
	double *dval = m->dval;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Iterate over rows
	for(r = 0; r < m->rows; r++){
		y[r] = dval[r] * x[r];
		for(j = ptr[r]; j < ptr[r+1]; j++){
			c = col[j];
			y[r] += val[j] * x[c];
			y[c] += val[j] * x[r];
		}
	}
}

/**
 * Stores on vector y the product of matrix A and vector x:
 *   y = Ax
 * Slightly optimized version.
 */
void sss_mvp2_sym(SSS_Matrix *m, double *x, double *y){
	register int r,j,c;
	register double y_r, x_r, a_rc;
	double *dval = m->dval;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Iterate over rows.
	for(r = 0; r < m->rows; r++){
		x_r = x[r];
		y_r = dval[r] * x_r;
		// Loop over non-zero elements in row r.
		for(j = ptr[r]; j < ptr[r+1]; j++){
			c = col[j];
			a_rc = val[j];
			y_r += a_rc * x[c];
			y[c] += a_rc * x_r;
		}
		y[r] = y_r;
	}
}

/**
 * Initialize matrix m.
 */
void sss_init_matrix(SSS_Matrix *m){
	m->dval = NULL;
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void sss_free_matrix(SSS_Matrix *m){
	if(m->dval != NULL) free(m->dval);
	if(m->val != NULL) free(m->val);
	if(m->col != NULL) free(m->col);
	if(m->ptr != NULL) free(m->ptr);
	sss_init_matrix(m);
}

/**
 * Prints matrix m data on console.
 */
void sss_print_matrix(SSS_Matrix *m){
	unsigned int i,j;
	unsigned int nzod = m->nz - m->rows;
	printf("dval= ");
	for(i=0; i < m->rows; i++){
		printf("%.2g ", m->dval[i]);
	}
	printf("\nval= ");
	j = 1;
	while(0 == m->ptr[j]) j++;
	for(i=0; i < nzod; i++){
		if(i == m->ptr[j]){
			printf("|");
			while(i == m->ptr[j]) j++;
		}
		printf("%.2g ", m->val[i]);
	}
	j = 1;
	while(0 == m->ptr[j]) j++;
	printf("\ncol= ");
	for(i=0; i < nzod; i++){
		if(i == m->ptr[j]){
			printf("|");
			while(i == m->ptr[j]) j++;
		}
		printf("%d ", m->col[i]);
	}
	printf("\nptr= ");
	for(i=0; i < m->rows+1; i++){
		printf("%d ", m->ptr[i]);
	}
	printf("\n");
	printf("\tnz= %d\n", m->nz);
	printf("\trows= %d\n", m->rows);
	printf("\tcols= %d\n", m->cols);
}

/**
 * Copies the COO matrix 'in' into SSS matrix 'out'.
 */
void coo2sss(COO_Matrix *in, SSS_Matrix *out){
	unsigned int i,nzod,ndcount,tot;
	COO_Matrix coo;
	
	// Init mem for SSS matrix
	nzod = in->nz - in->rows; // Non zero elements outside the main diagonal.
	out->dval = (double *) malloc(in->rows * sizeof(double));
	out->val = (double *) malloc(nzod * sizeof(double));
	out->col = (unsigned int *) malloc(nzod * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy in matrix
	coo_copy(in,&coo);
	coo_reorder_by_rows(&coo);
	
	// Copy elements row by row
	tot = 0;
	ndcount = 0;
	out->ptr[0] = ndcount;
	for(i = 0; i < coo.rows; i++){
		// Copy full row (minus last element)
		while(tot < coo.nz && coo.row[tot+1] == i){
			out->val[ndcount] = coo.val[tot];
			out->col[ndcount] = coo.col[tot];
			ndcount++;
			tot++;
		}
		// Copy last element in row to diagonal
		out->dval[i] = coo.val[tot++];
		
		out->ptr[i+1] = ndcount;
	}
	
	coo_free_matrix(&coo);
}

/**
 * Load SSS matrix from file filename.
 * File must be in Matrix Market format.
 */
int sss_load_matrix(char* filename, SSS_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	if(!sym){
		fprintf(stderr, "SSS format doesn't support non-symmetric matrices.\n");
		exit(EXIT_FAILURE);
	}
	coo2sss(&temp,m);
	coo_free_matrix(&temp);
	return sym;
}

#endif // SP_SSS