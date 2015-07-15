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

#ifndef SP_COO
#define SP_COO

#include "mmio.h"

/**
 * Sparse matrix structure stored
 * in coordinate list (COO) format.
 */
typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *row; // Corresponding row coordinate.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} COO_Matrix;

void coo_mvp(COO_Matrix *m, double *x, double *y);
void coo_mvp_sym(COO_Matrix *m, double *x, double *y);
void coo_init_matrix(COO_Matrix *m);
void coo_free_matrix(COO_Matrix *m);
void coo_print_matrix(COO_Matrix *m);
int coo_load_matrix(char* filename, COO_Matrix *coo);
void coo_copy(COO_Matrix *in, COO_Matrix *out);
void coo_reorder_by_rows(COO_Matrix *m);

void cooBottomUpMerge(
	unsigned int *A, unsigned int *A2, double *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, double *B3 );
void cooCopyArray(unsigned int *B, unsigned int *A, unsigned int n);
void cooCopyArrayDouble(double *B, double *A, unsigned int n);
unsigned int cooMergeMin(unsigned int a, unsigned int b);
void cooMergeSort(unsigned int *A, unsigned int *A2, double *A3, unsigned int n);

/**
 * Stores on vector y the product of matrix m and vector x:
 *   y = mx
 */
void coo_mvp(COO_Matrix *m, double *x, double *y){
	unsigned int i;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *col = m->col;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	for(i = 0; i < m->nz; i++){
		y[row[i]] += val[i]*x[col[i]];
	}
}

/**
 * Stores on vector y the product of symmetric matrix m and vector x:
 *   y = mx
 */
void coo_mvp_sym(COO_Matrix *m, double *x, double *y){
	unsigned int i,r,c;
	double *val = m->val;
	unsigned int *row = m->row;
	unsigned int *col = m->col;
	
	for(i = 0; i < m->rows; i++){
		y[i] = 0.0;
	}
	
	for(i = 0; i < m->nz; i++){
		r = row[i];
		c = col[i];
		y[r] += val[i]*x[c];
		if(r != c){
			y[c] += val[i]*x[r];
		}
	}
}

/**
 * Initialize matrix m.
 */
void coo_init_matrix(COO_Matrix *m){
	m->val = NULL;
	m->row = NULL;
	m->col = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void coo_free_matrix(COO_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->row != NULL) free(m->row);
	if(m->col != NULL) free(m->col);
	m->val = NULL;
	m->row = NULL;
	m->col = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Prints matrix m data on console.
 */
void coo_print_matrix(COO_Matrix *m){
	unsigned int i;
	printf("val= ");
	for(i=0; i < m->nz; i++){
		printf("%.2g ", m->val[i]);
	}
	printf("\nrow= ");
	for(i=0; i < m->nz; i++){
		printf("%d ", m->row[i]);
	}
	printf("\ncol= ");
	for(i=0; i < m->nz; i++){
		printf("%d ", m->col[i]);
	}
	printf("\nnz= %d\n", m->nz);
	printf("rows= %d\n", m->rows);
	printf("cols= %d\n", m->cols);
}

/**
 * Load COO matrix from file filename.
 * File must be in Matrix Market format.
 */
int coo_load_matrix(char* filename, COO_Matrix *coo){
	FILE *file;
	MM_typecode code;
	int m,n,nz,i,ival;
	unsigned int dc = 0;
	
	// Open file
	file = fopen(filename, "r");
	if(file == NULL){
		fprintf(stderr, "ERROR: can't open file \"%s\"\n", filename);
		exit(EXIT_FAILURE);
	}
	
	if (mm_read_banner(file, &code) != 0){
		fprintf(stderr, "ERROR: Could not process Matrix Market banner.\n");
		exit(EXIT_FAILURE);
	}
	
	// Check if input matrix is supported
	if(!mm_is_matrix(code) || !mm_is_sparse(code) || !mm_is_coordinate(code) ){
		fprintf(stderr, "Sorry, this application does not support ");
		fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(code));
		exit(EXIT_FAILURE);
	}
	
	// Read Matrix size
	if ((mm_read_mtx_crd_size(file, &m, &n, &nz)) != 0){
		fprintf(stderr, "ERROR: Could not read matrix size.\n");
		exit(EXIT_FAILURE);
	}
	coo->rows = m;
	coo->cols = n;
	coo->nz = nz;
	
	// Reserve space for matrix data
	coo->val = (double *) malloc(nz * sizeof(double));
	coo->row = (unsigned int *) malloc(nz * sizeof(unsigned int));
	coo->col = (unsigned int *) malloc(nz * sizeof(unsigned int));
	
	// Read matrix
	for(i = 0; i < nz; i++){
		if(mm_is_pattern(code)){
			if(fscanf(file, "%d %d\n", &(coo->row[i]), &(coo->col[i]) ) < 2){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
			coo->val[i] = 1.0;
		}else if (mm_is_real(code)){
			if(fscanf(file, "%d %d %lg\n", &(coo->row[i]), &(coo->col[i]), &(coo->val[i])) < 3){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
		}else if (mm_is_integer(code)){
			if(fscanf(file, "%d %d %d\n", &(coo->row[i]), &(coo->col[i]), &ival) < 3){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
			coo->val[i] = (double)ival;
		}
		coo->row[i]--;
		coo->col[i]--;
		if(coo->row[i] == coo->col[i]) { ++dc; }
	}
	
	fclose(file);
	
	return mm_is_symmetric(code);
}

void coo_copy(COO_Matrix *in, COO_Matrix *out){
	unsigned int i;
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->row = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	for(i = 0; i < in->nz; i++){
		out->val[i] = in->val[i]; 
	}
	for(i = 0; i < in->nz; i++){
		out->row[i] = in->row[i]; 
	}
	for(i = 0; i < in->nz; i++){
		out->col[i] = in->col[i]; 
	}
}

void coo_reorder_by_rows(COO_Matrix *m){
	cooMergeSort(m->row, m->col, m->val, m->nz);
}

/**************************************
 * Merge sort utility functions.
 * Code adapted from:
 * http://en.wikipedia.org/wiki/Merge_sort
 *************************************/
void cooBottomUpMerge(
	unsigned int *A, unsigned int *A2, double *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, double *B3 ) {
	
	unsigned int i0 = iLeft;
	unsigned int i1 = iRight;
	unsigned int j;
 
	/* While there are elements in the left or right runs */
	for (j = iLeft; j < iEnd; j++){
		/* If left run head exists and is <= existing right run head */
		if (i0 < iRight && (i1 >= iEnd || A[i0] <= A[i1])){
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
 
void cooCopyArray(unsigned int *B, unsigned int *A, unsigned int n) {
	unsigned int i;
   for(i = 0; i < n; i++){
        A[i] = B[i];
	}
}

void cooCopyArrayDouble(double *B, double *A, unsigned int n) {
	unsigned int i;
   for(i = 0; i < n; i++){
        A[i] = B[i];
	}
}

inline unsigned int cooMergeMin(unsigned int a, unsigned int b){
	return (a<b)?a:b;
}

/**
 * Use MergeSort algorithm to order the elements of A, from high to low.
 * Elements of A2 and A3 are moved correspondingly to elements of A.
 */
void cooMergeSort(unsigned int *A, unsigned int *A2, double *A3, unsigned int n) {
	/* array A[] has the items to sort; array B[] is a work array */
	unsigned int width,i;
	unsigned int *t;
	double *td;
	unsigned int swaps = 0;
	unsigned int *B = (unsigned int *)malloc(n * sizeof(unsigned int));
	unsigned int *B2 = (unsigned int *)malloc(n * sizeof(unsigned int));
	double *B3 = (double *)malloc(n * sizeof(double));
	/* Each 1-element run in A is already "sorted". */
	/* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted. */
	for (width = 1; width < n; width = 2 * width){
      /* Array A is full of runs of length width. */
      for (i = 0; i < n; i = i + 2 * width){
         /* Merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
         /* or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
         cooBottomUpMerge(A,A2,A3, i, cooMergeMin(i+width, n), cooMergeMin(i+2*width, n), B,B2,B3);
      }
      /* Now work array B is full of runs of length 2*width. */
      /* Copy array B to array A for next iteration. */
      /* A more efficient implementation would swap the roles of A and B */
      //cooCopyArray(B, A, n);
		//cooCopyArray(B2, A2, n);
		//cooCopyArrayDouble(B3, A3, n);
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		td = B3; B3 = A3; A3 = td;
		swaps++;
      /* Now array A is full of runs of length 2*width. */
   }
	if(swaps%2 == 1){
		//printf("Swap A and B\n");
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		td = B3; B3 = A3; A3 = td;
		cooCopyArray(B, A, n);
		cooCopyArray(B2, A2, n);
		cooCopyArrayDouble(B3, A3, n);
	}
	
	free(B);
	free(B2);
	free(B3);
}

/**************************************
 * End of merge sort utility functions.
 *************************************/

#endif // SP_COO
