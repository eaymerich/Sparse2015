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

#ifndef SP_TJDS
#define SP_TJDS

#include "coo.h"

/**
 * Sparse matrix structure stored
 * in transposed jagged diagonal storage (TJDS) format.
 */
typedef struct {
	double *val;        // Value of each non-zero element.
	unsigned int *row;  // Corresponding row coordinate.
	unsigned int *ptr;  // Initial position of each transposed jagged diagonal.
	unsigned int *perm; // Permutation vector. (Not needed if we store the diagonal as the first tjd. TODO: fix later.)
	unsigned int nz;    // Total non-zero elements.
	unsigned int rows;  // Total rows of matrix.
	unsigned int cols;  // Total columns of matrix.
	unsigned int diags; // Total transposed jagged diagonals of matrix.
} TJDS_Matrix;

void tjds_mvp(TJDS_Matrix *m, double *x, double *y);
void tjds_mvp_np(TJDS_Matrix *m, double *x, double *y);
void tjds_mvp_np_sym(TJDS_Matrix *m, double *x, double *y);
void tjds_init_matrix(TJDS_Matrix *m);
void tjds_free_matrix(TJDS_Matrix *m);
void tjds_print_matrix(TJDS_Matrix *m);
void tjds_permute_vector(TJDS_Matrix *m, double *x);
void tjds_permute_back_vector(TJDS_Matrix *m, double *x);
void coo2tjds(COO_Matrix *in, TJDS_Matrix *out);
int tjds_load_matrix(char* filename, TJDS_Matrix *m);

// Utility functions used to create TJDS matrix data.
void MergeSort(unsigned int *A, unsigned int *A2, unsigned int *A3, unsigned int n);
unsigned int MergeMin(unsigned int a, unsigned int b);
void CopyArray(unsigned int *B, unsigned int *A, unsigned int n);
void BottomUpMerge(
	unsigned int *A, unsigned int *A2, unsigned int *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, unsigned int *B3 );

/**
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 */
void tjds_mvp(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	// Creates permuted vector x.
	double *px = (double *)malloc(m->cols * sizeof(double));
	for(i = 0; i < m->cols; i++){
		px[i] = x[m->perm[i]];
	}
	
	for(i = 0; i < m->cols; i++){
		y[i] = 0.0;
	}
	
	// Multiply
	for(i = 0; i < m->diags; i++){
		j = 0;
		for(k = ptr[i]; k < ptr[i+1]; k++){
			y[row[k]] += val[k] * px[j];
			j++;
		}
	}
	
	free(px);
}

/**
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 * Assumes that x is already permuted to work with TJDS matrix m.
 */
void tjds_mvp_np(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k,end;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	for(i = 0; i < m->cols; i++){
		y[i] = 0.0;
	}
	
	// Multiply
	end = 0;
	for(i = 0; i < m->diags; i++){
		j = 0;
		k = end;
		end = ptr[i+1];
		for( ; k < end; k++){
			y[row[k]] += val[k] * x[j];
			j++;
		}
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 * Assumes that x is already permuted to work with TJDS matrix m.
 * This function works with matrices that may have zero elements
 * in the diagonal.
 */
void tjds_mvp_np_sym(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k,r,end;
	double x_r, a_ij;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	// Clear vector y.
	for(i = 0; i < m->cols; i++){
		y[i] = 0.0;
	}
	
	// Multiply
	end = 0;
	for(i = 0; i < 1; i++){ // Handle diagonal
		j = 0;
		k = end;
		end = ptr[i+1];
		for( ; k < end; k++){
			r = row[k];
			y[r] = val[k] * x[j];
			if(r != j){
				y[j] += val[k] * x[r];
			}
			j++;
		}
	}
	for(i = 1; i < m->diags; i++){ // Handle everything else
		j = 0;
		k = end;
		end = ptr[i+1];
		for( ; k < end; k++){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y[j] += a_ij * x[r];
			j++;
		}
	}
}
// void tjds_mvp_np_sym(TJDS_Matrix *m, double *x, double *y){
	// unsigned int i,j,k,r,end;
	// double *val = m->val;
	// unsigned int *row = m->row;
	// unsigned int *ptr = m->ptr;
	
	// // Clear vector y.
	// for(i = 0; i < m->cols; i++){
		// y[i] = 0.0;
	// }
	
	// // Multiply
	// end = 0;
	// for(i = 0; i < m->diags; i++){
		// j = 0;
		// k = end;
		// end = ptr[i+1];
		// for( ; k < end; k++){
			// r = row[k];
			// y[r] += val[k] * x[j];
			// if(r != j){
				// y[j] += val[k] * x[r];
			// }
			// j++;
		// }
	// }
// }

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Assumes that x is already permuted to work with TJDS matrix m.
 * Assumes that all elements in the diagonal of matrix m are non-zero.
 * Makes light use of C pointer arithmetic.
 */
void tjds_mvp2_np_sym(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k,r;
	double a_ij;
	double *val = m->val;
	unsigned int *row;
	unsigned int *ptr = m->ptr;
	
	// First multiply all elements in diagonal.
	for(i = 0; i < m->cols; i++){
		y[i] = (*val++)*x[i];
	}

	// Now multiply other elements.
	row = &(m->row[1]);
	for(i = 1; i < m->diags; i++){
		j = 0;
		for(k = ptr[i]; k < ptr[i+1]; k++){
			r = *row++;				// row index r
			a_ij = *val++;			// matrix element A(i,j)
			y[j] += a_ij * x[r];
			y[r] += a_ij * x[j];
			j++;
		}
	}
}
// void tjds_mvp2_np_sym(TJDS_Matrix *m, double *x, double *y){
	// unsigned int i,j,k,r;
	// double a_ij;
	// double *val = m->val;
	// unsigned int *row = m->row;
	// unsigned int *ptr = m->ptr;
	
	// // First multiply all elements in diagonal.
	// for(i = 0; i < m->cols; i++){
		// y[i] = val[i]*x[i];
	// }

	// // Now multiply other elements.
	// for(i = 1; i < m->diags; i++){
		// j = 0;
		// for(k = ptr[i]; k < ptr[i+1]; k++){
			// r = row[k];
			// a_ij = val[k];
			// y[j] += a_ij * x[r];
			// y[r] += a_ij * x[j];
			// j++;
		// }
	// }
// }

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Assumes that x is already permuted to work with TJDS matrix m.
 * Assumes that all elements in the diagonal of matrix m are non-zero.
 * Makes strong use of C pointer arithmetic.
 */
void tjds_mvpP_np_sym(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,k,r,end;
	double a_ij;
	double *val = m->val,*xx=x,*yy=y;
	unsigned int *row;
	unsigned int *ptr = &(m->ptr[1]);
	
	// First multiply all elements in diagonal.
	for(i = 0; i < m->cols; i++){
		(*yy++) = (*val++)*(*xx++);
	}
	
	// Now multiply other elements.
	end = (*ptr++);
	row = &(m->row[end]);
	for(i = 1; i < m->diags; i++){
		xx = x;
		yy = y;
		k = end;
		end = (*ptr++);
		for( ; k < end; k++){
			r = (*row++);
			a_ij = (*val++);
			(*yy++) += a_ij * x[r];
			y[r] += a_ij * (*xx++);
		}
	}
}

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Assumes that x is already permuted to work with TJDS matrix m.
 * Assumes that all elements in the diagonal of matrix m are non-zero.
 * Caches first element.
 */
void tjds_mvpA_np_sym(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k,r,end;
	double a_ij,y_0=0.0;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	// First multiply all elements in diagonal.
	for(i = 0; i < m->cols; i++){
		y[i] = val[i]*x[i];
	}
	
	// Now multiply other elements.
	end = ptr[1];
	for(i = 1; i < m->diags; i++){
		j = 0;
		k = end;
		end = ptr[i+1];
		
		// Cache first element
		r = row[k];
		a_ij = val[k];
		y[r] += a_ij * x[j];
		y_0 += a_ij * x[r];
		k++;
		j++;
		
		// Handle all others
		for( ; k < end; k++){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y[j] += a_ij * x[r];
			j++;
		}
	}
	
	y[0] += y_0;
}

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Assumes that x is already permuted to work with TJDS matrix m.
 * Assumes that all elements in the diagonal of matrix m are non-zero.
 * Caches first two elements.
 */
void tjds_mvpB_np_sym(TJDS_Matrix *m, double *x, double *y){
	unsigned int i,j,k,r,end;
	double a_ij,y_0=0.0,y_1=0.0;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	// First multiply all elements in diagonal.
	for(i = 0; i < m->cols; i++){
		y[i] = val[i]*x[i];
	}
	
	// Now multiply other elements.
	end = ptr[1];
	for(i = 1; i < m->diags; i++){
		j = 0;
		k = end;
		end = ptr[i+1];
		
		// Cache first element
		r = row[k];
		a_ij = val[k];
		y[r] += a_ij * x[j];
		y_0 += a_ij * x[r];
		k++;
		j++;
		
		// Cache second element (if any)
		if(k < end){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y_1 += a_ij * x[r];
			k++;
			j++;
		}
		
		// Handle all others
		for( ; k < end; k++){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y[j] += a_ij * x[r];
			j++;
		}
	}
	
	y[0] += y_0;
	y[1] += y_1;
}

/**
 * Stores on vector y the product of symmetric matrix A and vector x:
 *   y = Ax
 * Assumes that x is already permuted to work with TJDS matrix m.
 * Assumes that all elements in the diagonal of matrix m are non-zero.
 * Caches first two elements.
 */
void tjds_mvpC_np_sym(TJDS_Matrix *m, double *x, double *y){
	register unsigned int i,j,k,r,end;
	register double a_ij,y_0=0.0,y_1=0.0,y_2=0.0;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *ptr = m->ptr;
	
	// First multiply all elements in diagonal.
	for(i = 0; i < m->cols; i++){
		y[i] = val[i]*x[i];
	}
	
	// Now multiply other elements.
	end = ptr[1];
	for(i = 1; i < m->diags; i++){
		j = 0;
		k = end;
		end = ptr[i+1];
		
		// Cache first element
		r = row[k];
		a_ij = val[k];
		y[r] += a_ij * x[j];
		y_0 += a_ij * x[r];
		k++;
		j++;
		
		// Cache second element (if any)
		if(k < end){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y_1 += a_ij * x[r];
			k++;
			j++;
		}
		
		// Cache second element (if any)
		if(k < end){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y_2 += a_ij * x[r];
			k++;
			j++;
		}
		
		// Handle all others
		for( ; k < end; k++){
			r = row[k];
			a_ij = val[k];
			y[r] += a_ij * x[j];
			y[j] += a_ij * x[r];
			j++;
		}
	}
	
	y[0] += y_0;
	y[1] += y_1;
	y[2] += y_2;
}

/**
 * Initialize matrix m.
 */
void tjds_init_matrix(TJDS_Matrix *m){
	m->val = NULL;
	m->row = NULL;
	m->ptr = NULL;
	m->perm = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void tjds_free_matrix(TJDS_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->row != NULL) free(m->row);
	if(m->ptr != NULL) free(m->ptr);
	if(m->perm != NULL) free(m->perm);
	m->val = NULL;
	m->row = NULL;
	m->ptr = NULL;
	m->perm = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Prints matrix m data on console.
 */
void tjds_print_matrix(TJDS_Matrix *m){
	unsigned int i,j;
	j = 1;
	// printf("val= ");
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
	// for(i=0; i < m->diags+1; i++){
		// printf("%d ", m->ptr[i]);
	// }
	// printf("\nperm= ");
	// for(i=0; i < m->cols; i++){
		// printf("%d ", m->perm[i]);
	// }
	// printf("\n");
	printf("\tnz= %d\n", m->nz);
	printf("\trows= %d\n", m->rows);
	printf("\tcols= %d\n", m->cols);
	printf("\tdiags= %d\n", m->diags);
}

/**
 * Permutes vector x so it works correctly 
 * when multiplying it with TJDS matrix m.
 */
void tjds_permute_vector(TJDS_Matrix *m, double *x){
	unsigned int i;
	double *t = (double *)malloc(m->cols * sizeof(double));
	for(i = 0; i < m->cols; i++){
		t[i] = x[m->perm[i]];
	}
	for(i = 0; i < m->cols; i++){
		x[i] = t[i];
	}
	free(t);
}

/**
 * Reverse the work of the previous function.
 * Permutes back vector x to its original order.
 */
void tjds_permute_back_vector(TJDS_Matrix *m, double *x){
	unsigned int i;
	double *t = (double *)malloc(m->cols * sizeof(double));
	for(i = 0; i < m->cols; i++){
		t[m->perm[i]] = x[i];
	}
	for(i = 0; i < m->cols; i++){
		x[i] = t[i];
	}
	free(t);
}

/**
 * Copies the COO matrix 'in' into TJDS matrix 'out'.
 */
void coo2tjds(COO_Matrix *in, TJDS_Matrix *out){
	unsigned int i,j,count,diags,end, col_begin;
	unsigned int *col_index = (unsigned int *)malloc(in->cols * sizeof(unsigned int));
	unsigned int *col_count = (unsigned int *)malloc(in->cols * sizeof(unsigned int));
	unsigned int *col_start = (unsigned int *)malloc(in->cols * sizeof(unsigned int));
	for(i = 0; i < in->cols; i++){
		col_index[i] = i;
	}
	
	// Count elements per column and
	// record column start.
	j = 0;
	for(i = 0; i < in->cols; i++){
		col_start[i] = j;
		count=0;
		while(j < in->nz && in->col[j] == i){
			count++;
			j++;
		}
		col_count[i] = count;
	}
	
	// printf("col_index= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_index[i]);
	// }
	// printf("\n");
	
	// printf("col_count= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_count[i]);
	// }
	// printf("\n");
	
	// printf("col_start= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_start[i]);
	// }
	// printf("\n");
	
	// Sort columns from max to min number of elements
	MergeSort(col_count, col_index, col_start, in->cols);
	
	// printf("col_index= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_index[i]);
	// }
	// printf("\n");
	
	// printf("col_count= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_count[i]);
	// }
	// printf("\n");
	
	// printf("col_start= ");
	// for(i = 0; i < in->cols; i++){
		// printf("%d ", col_start[i]);
	// }
	// printf("\n");
	
	diags = col_count[0];
	// printf("Number of diags=%d\n", diags);
	
	// Allocate memory for TJDS data
	out->val = (double *)malloc(in->nz * sizeof(double));
	out->row = (unsigned int *)malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *)malloc((diags+1) * sizeof(unsigned int));
	out->perm = col_index;
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	out->diags = diags;
	
	// Copy elements
	end = in->cols-1;
	count = 0;
	for(i = 0; i < diags; i++){
		out->ptr[i] = count;
		while(i >= col_count[end]){
			end--;
		}
		// printf("Segment %d has %d elements= ", i, end+1);
		for(j = 0; j <= end; j++){
			col_begin = col_start[j];
			// printf("%.1f ", in->val[col_begin+i]);
			// Save value and row
			out->val[count] = in->val[col_begin+i];
			out->row[count] = in->row[col_begin+i];
			count++;
		}
		// printf("\n");
	}
	out->ptr[i] = count;
	
	// Permutate rows using segment 1
	if(diags > 1){
		for(i = 0; i < in->cols; i++){
			col_count[out->row[i]] = i; // reuse col_count to save position of each row
		}
		for(i = out->ptr[1]; i < out->nz; i++){
			out->row[i] = col_count[out->row[i]];
		}
	}
	
	free(col_count);
	free(col_start);
	//free(col_index);
}

/**
 * Load TJDS matrix from file filename.
 * File must be in Matrix Market format.
 */
int tjds_load_matrix(char* filename, TJDS_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	coo2tjds(&temp,m);
	coo_free_matrix(&temp);
	return sym;
}

/**************************************
 * Merge sort utility functions.
 * Code adapted from:
 * http://en.wikipedia.org/wiki/Merge_sort
 *************************************/
void BottomUpMerge(
	unsigned int *A, unsigned int *A2, unsigned int *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, unsigned int *B3 ) {
	
	unsigned int i0 = iLeft;
	unsigned int i1 = iRight;
	unsigned int j;
 
	/* While there are elements in the left or right runs */
	for (j = iLeft; j < iEnd; j++){
		/* If left run head exists and is <= existing right run head */
		if (i0 < iRight && (i1 >= iEnd || A[i0] >= A[i1])){
			B[j] = A[i0];
			B2[j] = A2[i0];
			B3[j] = A3[i0];
			i0++;
		} else {
			B[j] = A[i1];
			B2[j] = A2[i1];
			B3[j] = A3[i1];
			i1++;
		}
	}
}
 
void CopyArray(unsigned int *B, unsigned int *A, unsigned int n) {
	unsigned int i;
   for(i = 0; i < n; i++){
        A[i] = B[i];
	}
}

inline unsigned int MergeMin(unsigned int a, unsigned int b){
	return (a<b)?a:b;
}

/**
 * Use MergeSort algorithm to order the elements of A, from high to low.
 * Elements of A2 and A3 are moved correspondingly to elements of A.
 */
void MergeSort(unsigned int *A, unsigned int *A2, unsigned int *A3, unsigned int n) {
	/* array A[] has the items to sort; array B[] is a work array */
	unsigned int width,i;
	unsigned int *t;
	unsigned int swaps = 0;
	unsigned int *B = (unsigned int *)malloc(n * sizeof(unsigned int));
	unsigned int *B2 = (unsigned int *)malloc(n * sizeof(unsigned int));
	unsigned int *B3 = (unsigned int *)malloc(n * sizeof(unsigned int));
	/* Each 1-element run in A is already "sorted". */
	/* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted. */
	for (width = 1; width < n; width = 2 * width){
      /* Array A is full of runs of length width. */
      for (i = 0; i < n; i = i + 2 * width){
         /* Merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
         /* or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
         BottomUpMerge(A,A2,A3, i, MergeMin(i+width, n), MergeMin(i+2*width, n), B,B2,B3);
      }
      /* Now work array B is full of runs of length 2*width. */
      /* Copy array B to array A for next iteration. */
      /* A more efficient implementation would swap the roles of A and B */
      //CopyArray(B, A, n);
		//CopyArray(B2, A2, n);
		//CopyArray(B3, A3, n);
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		t = B3; B3 = A3; A3 = t;
		swaps++;
      /* Now array A is full of runs of length 2*width. */
   }
	if(swaps%2 == 1){
		//printf("Swap A and B\n");
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		t = B3; B3 = A3; A3 = t;
		CopyArray(B, A, n);
		CopyArray(B2, A2, n);
		CopyArray(B3, A3, n);
	}
	free(B);
	free(B2);
	free(B3);
}

/**************************************
 * End of merge sort utility functions.
 *************************************/

#endif // SP_TJDS
