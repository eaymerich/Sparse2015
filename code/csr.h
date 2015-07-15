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
 * in Compressed Sparse Row (CSR) format.
 * In case of symmetric matrices, only
 * the lower triangle is stored.
 */
typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int *ptr; // Initial position of each row.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} CSR_Matrix;

void csr_mvp(CSR_Matrix *m, double *x, double *y);
void csr_mvp2(CSR_Matrix *m, double *x, double *y);
void csr_mvp_oski(CSR_Matrix *m, double *x, double *y);
void csr_mvp_oski2(CSR_Matrix *m, double *x, double *y);
void csr_mvp_sym(CSR_Matrix *m, double *x, double *y);
void csr_mvp_sym2(CSR_Matrix *m, double *x, double *y);
void csr_mvp_sym_oski_lo(CSR_Matrix *m, double *x, double *y);
void csr_mvp_sym_oski2(CSR_Matrix *m, double *x, double *y);

void csr_init_matrix(CSR_Matrix *m);
void csr_free_matrix(CSR_Matrix *m);
void csr_print_matrix(CSR_Matrix *m);
void coo2csr(COO_Matrix *in, CSR_Matrix *out);
void coo22csr(COO_Matrix *in, CSR_Matrix *out);
void coo222csr(COO_Matrix *in, CSR_Matrix *out);
int csr_load_matrix(char* filename, CSR_Matrix *m);

/**
 * Stores on vector y the product of matrix M and vector x:
 *   y = Mx
 */
void csr_mvp(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,end;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	end = 0;
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
		j = end;
		end = ptr[i+1];
		for( ; j < end; j++){
			y[i] += val[j] * x[col[j]];
		}
	}
}

/**
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 * Uses Berkeley optimization.
 */
void csr_mvp2(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,end;
	double tempy;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	end = 0;
	
	// Loop over rows.
	for(i = 0; i < m->rows; i++){
		tempy = 0.0;
		j = end;
		end = ptr[i+1];
		// Loop over non-zero elements in row i.
		for( ; j < end; j++){
			tempy += val[j] * x[col[j]];
		}
		y[i] = tempy;
	}
}

/**
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 * Uses exactly Berkeley implementation.
 */
void csr_mvp_oski(CSR_Matrix *m, double *x, double *y){
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
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 * Uses a bit improved Berkeley implementation.
 */
void csr_mvp_oski2(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,end;
	double y_i;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Loop over rows.
	end = 0;
	for(i = 0; i < m->rows; i++){
		y_i = 0.0;
		j = end;
		end = ptr[i+1];
		// Loop over non-zero elements in row i.
		for(; j < end; j++){
			y_i += (*val++) * x[*col++];
		}
		y[i] = y_i;
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 */
void csr_mvp_sym(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,mcol;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
		for(j = ptr[i]; j < ptr[i+1]; j++){
			mcol = col[j];
			y[i] += val[j] * x[mcol];
			if(i != mcol){
				y[mcol] += val[j] * x[i];
			}
		}
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 * Uses Berkeley optimization.
 */
void csr_mvp_sym2(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,end,mcol;
	double tempy;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	end = 0;
	for(i = 0; i < m->rows; i++){
		tempy = 0.0;
		j = end;
		end = ptr[i+1];
		for( ; j < end; j++){
			mcol = col[j];
			tempy += val[j] * x[mcol];
			if(i != mcol){
				y[mcol] += val[j] * x[i];
			}
		}
		y[i] = tempy;
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 * Uses Berkeley implementation, except that this one assumes that only
 * the lower triangle is stored
 */
void csr_mvp_sym_oski_lo(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,jj;
	double y_i, x_i, a_ij;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Loop over rows.
	for(i = 0; i < m->rows; i++){
		y_i = 0.0;
		x_i = x[i];
		
		// Loop over non-zeros in row i.
		jj = ptr[i];
		for(; jj < ptr[i+1]-1; jj++){
			j = *col++; // column index j
			a_ij = *val++;  // matrix element A(i,j)
			y_i += a_ij * x[j];
			y[j] += a_ij * x_i;
		}
		
		// Special handling of the diagonal
		if(jj < ptr[i+1]){ // Check if this row is not empty
			if(i == *col){
				// If this is diagonal
				y_i += (*val++) * x_i;
				col++;
			}else{
				j = *col++; // column index j
				a_ij = *val++;  // matrix element A(i,j)
				y_i += a_ij * x[j];
				y[j] += a_ij * x_i;
			}
		}
		
		y[i] = y_i;
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 * Uses Berkeley implementation variation, except that this one assumes that only
 * the lower triangle is stored
 */
void csr_mvp_sym_oski2(CSR_Matrix *m, double *x, double *y){
	unsigned int i,j,jj;
	double y_i, x_i, a_ij;
	double *val = m->val;
	unsigned int *col = m->col;
	unsigned int *ptr = m->ptr;
	
	// Loop over rows.
	for(i = 0; i < m->rows; i++){
		y_i = 0.0;
		x_i = x[i];
		
		// Loop over non-zeros in row i.
		for(jj = ptr[i]; jj < ptr[i+1]; jj++){
			j = *col++; // column index j
			a_ij = *val++;  // matrix element A(i,j)
			y_i += a_ij * x[j];
			y[j] += a_ij * x_i;
		}
		
		// No special handling for diagonal
		// The correct value is written over the incorrect value.
		
		y[i] = y_i;
	}
}

/**
 * Initialize matrix m.
 */
void csr_init_matrix(CSR_Matrix *m){
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void csr_free_matrix(CSR_Matrix *m){
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
void csr_print_matrix(CSR_Matrix *m){
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
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 */
void coo2csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i,j;
	unsigned int tot = 0;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy elements row by row
	out->ptr[0] = tot;
	for(i = 0; i < out->rows; i++){
		for(j = 0; j < in->nz; j++){
			if(in->row[j] == i){
				out->val[tot] = in->val[j];
				out->col[tot] = in->col[j];
				tot++;
			}
		}
		out->ptr[i+1] = tot;
	}
}

/**
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 */
void coo222csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i;
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
	//coo_print_matrix(&coo);
	coo_reorder_by_rows(&coo);
	//coo_print_matrix(&coo);
	
	// Copy elements row by row
	out->ptr[0] = tot;
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
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 * Faster than previous function. (I hope...)
 * TODO: FIX this function doesn't work
 */
void coo22csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i,j;
	unsigned int tot;
	unsigned int *col_start = (unsigned int *)malloc(in->cols * sizeof(unsigned int));
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Find column starts.
	//printf("col_start= ");
	j = 0;
	for(i = 0; i < in->cols; i++){
		col_start[i] = j;
		//printf("%d ", col_start[i]);
		while(j < in->nz && in->col[j] == i){
			j++;
		}		
	}
	//printf("\n");
	
	// Copy elements row by row
	tot = 0;
	out->ptr[0] = tot;
	for(i = 0; i < in->rows; i++){
		//printf("For row %d\n",i);
		for(j = 0; j < in->cols; j++){
			//printf("\trow: %d\n",in->row[col_start[j]] );
			if(in->row[col_start[j]] == i){
				// Copy element
				out->val[tot] = in->val[col_start[j]];
				out->col[tot] = in->col[col_start[j]];
				tot++;
				// Advance column start
				col_start[j]++;
			}
		}
		out->ptr[i+1] = tot;
	}
	
	free(col_start);
}

/**
 * Load CSR matrix from file filename.
 * File must be in Matrix Market format.
 */
int csr_load_matrix(char* filename, CSR_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	coo222csr(&temp,m);
	coo_free_matrix(&temp);
	return sym;
}

#endif // SP_CSR