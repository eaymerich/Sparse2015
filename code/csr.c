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

#include <stdio.h>
#include <stdlib.h>

#include "mmio.h"
#include "clock.h"
#include "coo.h"
#include "csr.h"

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);

int main(int argc, char *argv[]){
	unsigned int i;
	int sym;
	CSR_Matrix csr;
	Clock clock;
	char* filename;
	double *x;
	double *y;
	
	parse_args(argc, argv);
	filename = argv[1];
	
	// Initialize matrices.
	csr_init_matrix(&csr);
	
	// Load matrix from file into COO.
	printf("Loading matrix \"%s\"\n", filename);
	sym = csr_load_matrix(filename, &csr);
	
	if(sym) printf("Matrix is symmetric\n");
	else	printf("Matrix is general (non-symmetric)\n");
	
	// Print Matrix data
	//printf("CSR matrix data:\n");
	//csr_print_matrix(&csr);
	
	// Initialize vectors for multiplication
	x = (double *)malloc(csr.cols * sizeof(double));
	y = (double *)malloc(csr.rows * sizeof(double));
	for(i = 0; i < csr.cols; i++){
		x[i] = i+1;
	}
	//print_vector("\nx= ", x, coo.cols);

	printf("Calculating matrix-vector product for %u iterations...\n", iterations);

	if(sym){
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp_sym(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp time: %.2fus \n", clock.sec);
				
		// Time CSR matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			csr_mvp_sym2(&csr,x,y);
		}
		clock_stop(&clock);
		printf("CSR mvp2 time: %.2fus \n\n", clock.sec);
				
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp_sym_oski_lo(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp oski time: %.2fus \n\n", clock.sec);
		
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp_sym_oski2(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp oski2 time: %.2gs \n", clock.sec);
		
		//print_vector("y= ", y, csr.rows);
	}else{
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp time: %.2gs \n", clock.sec);
		//print_vector("y= ", y, csr.cols);
		
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp2(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp2 time: %.2gs \n", clock.sec);
		//print_vector("y= ", y, csr.cols);
		
		// Time CSR matrix-vector product
		// clock_start(&clock);
		// for(i=0; i < iterations; i++){
			// csr_mvp_oski(&csr,x,y);
		// }
		// clock_stop(&clock);
		// printf("CSR mvp oski time: %.2gs \n", clock.sec);
		//print_vector("y= ", y, csr.cols);
		
		// Time CSR matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			csr_mvp_oski2(&csr,x,y);
		}
		clock_stop(&clock);
		printf("CSR mvp oski2 time: %.2gs \n\n", clock.sec);
		//print_vector("y= ", y, csr.cols);
	}
	
	// Free resources
	free(x);
	free(y);
	csr_free_matrix(&csr);
	
	exit(EXIT_SUCCESS);
}

void print_vector(char* pre, double *v, unsigned int size){
	unsigned int i;
	printf("%s", pre);
	for(i = 0; i < size; i++){
		//printf("%.1f ", v[i]);
		printf("%e \n", v[i]);
	}
	printf("\n");
}

void parse_args(int argc, char *argv[]){
	int i;
	if (argc < 2) {
		printf("Usage: %s input_matrix [iterations]\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	if(argc >= 3){
		i = atoi(argv[2]);
		if (i <= 0){
			printf("Invalid number of iterations.\n");
			exit(EXIT_FAILURE);
		}
		iterations = i;
	}
}
