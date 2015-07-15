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

#ifndef SP_CSR
#define SP_CSR

#include "coo.h"

/**
 * Sparse matrix structure stored
 * in Compressed Sparse Column (CSC) format.
 * In case of symmetric matrices, only
 * the lower triangle is stored.
 */
typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *row; // Corresponding row coordinate.
	unsigned int *ptr; // Initial position of each column.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} CSC_Matrix;

void csc_mvp(CSC_Matrix *m, double *x, double *y);
void csc_mvp_sym(CSC_Matrix *m, double *x, double *y);

void csc_init_matrix(CSC_Matrix *m);
void csc_free_matrix(CSC_Matrix *m);
void csc_print_matrix(CSC_Matrix *m);
void coo2csc(COO_Matrix *in, CSC_Matrix *out);
int csc_load_matrix(char* filename, CSC_Matrix *m);

void csc_mvp(CSC_Matrix *m, double *x, double *y){
	unsigned int i,j;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	for(j = 0; j < m->cols; j++){
		for(i = ptr[j] ; i < ptr[j+1]; i++){
			y[row[i]] += val[i] * x[j];
		}
	}
}

void csc_mvp_sym(CSC_Matrix *m, double *x, double *y){
	unsigned int i,j,mrow;
	double x_j, a_ij;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	for(j = 0; j < m->cols; j++){
		x_j = x[j];
		
		// Handle diagonal
		i = ptr[j];
		mrow = row[i];
		y[mrow] += val[i] * x_j;

		// Handle everything else
		i++;
		for(; i < ptr[j+1]; i++){
			mrow = row[i];
			a_ij = val[i];
			y[mrow] += a_ij * x_j;
			y[j] += a_ij * x[mrow];
		}
	}
}

void csc_mvp_sym2(CSC_Matrix *m, double *x, double *y){
	unsigned int i,j,mrow;
	double *val = m->val;
	double y_j,x_j,mval;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	for(j = 0; j < m->cols; j++){
		y_j = 0.0;
		x_j = x[j];
		i = ptr[j];
		
		// Handle diagonal
		y_j += (*val++) * x_j;
		row++;
		i++;
		
		// Handle everything else
		for(; i < ptr[j+1]; i++){
			mrow = (*row++);
			mval = (*val++);
			y[mrow] += mval * x_j;
			y_j += mval * x[mrow];
		}
		y[j] += y_j;
	}
}

void csc_init_matrix(CSC_Matrix *m){
	m->val = NULL;
	m->row = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

void csc_free_matrix(CSC_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->row != NULL) free(m->row);
	if(m->ptr != NULL) free(m->ptr);
	csc_init_matrix(m);
}

void csc_print_matrix(CSC_Matrix *m){
	// unsigned int i,j;
	// printf("val= ");
	// j = 1;
	// for(i=0; i < m->nz; i++){
		// if(i == m->ptr[j]){
			// printf("|");
			// j++;
		// }
		// printf("%.2g ", m->val[i]);
	// }
	// j = 1;
	// printf("\nrow= ");
	// for(i=0; i < m->nz; i++){
		// if(i == m->ptr[j]){
			// printf("|");
			// j++;
		// }
		// printf("%d ", m->row[i]);
	// }
	// printf("\nptr= ");
	// for(i=0; i < m->rows+1; i++){
		// printf("%d ", m->ptr[i]);
	// }
	printf("\tnz= %d\n", m->nz);
	printf("\trows= %d\n", m->rows);
	printf("\tcols= %d\n", m->cols);	
}

void coo2csc(COO_Matrix *in, CSC_Matrix *out){
	unsigned int i;
	unsigned int tot = 0;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->row = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
		
	// Copy elements row by row
	out->ptr[0] = tot;
	for(i = 0; i < in->cols; i++){
		while(tot < in->nz && in->col[tot] == i){
			out->val[tot] = in->val[tot];
			out->row[tot] = in->row[tot];
			tot++;
		}
		out->ptr[i+1] = tot;
	}
}

int csc_load_matrix(char* filename, CSC_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	coo2csc(&temp,m);
	coo_free_matrix(&temp);
	return sym;
}

#endif // SP_CSR
