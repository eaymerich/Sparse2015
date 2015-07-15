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

#ifndef SP_OSKI
#define SP_OSKI

#include "coo.h"

/**
 * Sparse matrix structure stored
 * in compressed row storage (CSR) format for OSKI
 * algorithms. In case of symmetric matrices, only
 * the upper triangle is stored.
 */
typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int *ptr; // Initial position of each row.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} OSKI_Matrix;

void oski_mvp(OSKI_Matrix *m, double *x, double *y);
void oski_mvp_sym(OSKI_Matrix *m, double *x, double *y);
void oski_mvp_sym2(OSKI_Matrix *m, double *x, double *y);

void oski_init_matrix(OSKI_Matrix *m);
void oski_free_matrix(OSKI_Matrix *m);
void oski_print_matrix(OSKI_Matrix *m);

void coo2oski(COO_Matrix *in, OSKI_Matrix *out);
void coo2oski_sym(COO_Matrix *in, OSKI_Matrix *out);
int oski_load_matrix(char* filename, OSKI_Matrix *m);

/**
 * Stores on vector y the product of matrix A and vector x:
 *   y = Ax
 * Uses exactly Berkeley implementation.
 */
void oski_mvp(OSKI_Matrix *m, double *x, double *y){
	unsigned int i,j;
	double y_i;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Loop over rows.
	for(i = 0; i < m->rows; i++){
		y_i = 0.0;
		// Loop over non-zero elements in row i.
		for(j = ptr[i]; j < ptr[i+1]; j++){
			y_i += (*val++) * x[*col++];
		}
		y[i] = y_i;
	}
}

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Uses exactly Berkeley implementation.
 */
void oski_mvp_sym(OSKI_Matrix *m, double *x, double *y){
	unsigned int i,j,jj;
	double y_i, x_i, a_ij;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	// Loop over rows.
	for(i = 0; i < m->rows; i++){
		y_i = y[i];
		x_i = x[i];
		jj = ptr[i];
		
		// Special handling of the diagonal
		if(i == *col){
			y_i += x[i] * (*val++);
			col++;
			jj++;
		}
		else y_i = 0.0;
		
		// Loop over non-zeroes in row i
		for( ; jj < ptr[i+1]; jj++){
			j = *col++;			// column index j
			a_ij = *val++;		// matrix element A(i,j)
			y_i += a_ij * x[j];
			y[j] += a_ij * x_i;
		}
		
		y[i] = y_i;
	}
}

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Uses exactly Berkeley implementation.
 */
void oski_mvp_sym2(OSKI_Matrix *m, double *x, double *y){
	unsigned int i,j,jj,end;
	double y_i, x_i, a_ij;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	// Loop over rows.
	end = 0;
	for(i = 0; i < m->rows; i++){
		x_i = x[i];
		jj = end;
		
		// Special handling of the diagonal
		if(i == *col){
			y_i = x[i] * (*val++);
			col++;
			jj++;
		}
		else y_i = 0.0;
		
		// Loop over non-zeroes in row i
		end = ptr[i+1];
		for( ; jj < end; jj++){
			j = *col++;			// column index j
			a_ij = *val++;		// matrix element M(i,j)
			y_i += a_ij * x[j];
			y[j] += a_ij * x_i;
		}
		
		y[i] += y_i;
	}
}

/**
 * Initialize matrix m.
 */
void oski_init_matrix(OSKI_Matrix *m){
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void oski_free_matrix(OSKI_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->col != NULL) free(m->col);
	if(m->ptr != NULL) free(m->ptr);
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Prints matrix m data on console.
 */
void oski_print_matrix(OSKI_Matrix *m){
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
	// printf("\ncol= ");
	// for(i=0; i < m->nz; i++){
		// if(i == m->ptr[j]){
			// printf("|");
			// j++;
		// }
		// printf("%d ", m->col[i]);
	// }
	// printf("\nptr= ");
	// for(i=0; i < m->rows+1; i++){
		// printf("%d ", m->ptr[i]);
	// }
	printf("\tnz= %d\n", m->nz);
	printf("\trows= %d\n", m->rows);
	printf("\tcols= %d\n", m->cols);
}

/**
 * Copies the COO matrix 'in' into OSKI matrix 'out'.
 */
void coo2oski(COO_Matrix *in, OSKI_Matrix *out){
	unsigned int i,j;
	unsigned int tot = 0;
	COO_Matrix coo;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy in matrix
	coo_copy(in,&coo);
	coo_reorder_by_rows(&coo);
	
	// Copy elements row by row
	out->ptr[0] = tot;
	j = 0;
	for(i = 0; i < coo.rows; i++){
		while(tot < coo.nz && coo.row[tot] == i){
			out->val[tot] = coo.val[tot];
			out->col[tot] = coo.col[tot];
			tot++;
		}
		out->ptr[i+1] = tot;
	}
	
	coo_free_matrix(&coo);
}

/**
 * Copies the COO matrix 'in' into OSKI matrix 'out'.
 * If 'in' is symmetric, it is converted into an upper triangular matrix.
 */
void coo2oski_sym(COO_Matrix *in, OSKI_Matrix *out){
	unsigned int *t;
	unsigned int i,j,tot = 0;
	COO_Matrix coo;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy in matrix
	coo_copy(in,&coo);
	
	// Change COO from lower triangle to upper triangle
	t = coo.row; 
	coo.row = coo.col;
	coo.col = t;
	coo_reorder_by_rows(&coo);
	
	// Copy elements row by row
	out->ptr[0] = tot;
	j = 0;
	for(i = 0; i < coo.rows; i++){
		while(tot < coo.nz && coo.row[tot] == i){
			out->val[tot] = coo.val[tot];
			out->col[tot] = coo.col[tot];
			tot++;
		}
		out->ptr[i+1] = tot;
	}
	
	coo_free_matrix(&coo);
}

/**
 * Load OSKI matrix from file filename.
 * File must be in Matrix Market format.
 */
int oski_load_matrix(char* filename, OSKI_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	if(sym){
		coo2oski_sym(&temp,m);
	}else{
		coo2oski(&temp,m);
	}
	coo_free_matrix(&temp);
	return sym;
}

#endif // SP_OSKI
