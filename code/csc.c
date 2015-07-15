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

#include "clock.h"
#include "csc.h"

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);

int main(int argc, char *argv[]){
	unsigned int i;
	int sym;
	CSC_Matrix m;
	Clock clock;
	char* filename;
	double *x;
	double *y;
	
	parse_args(argc, argv);
	filename = argv[1];
	
	// Initialize matrices.
	csc_init_matrix(&m);
	
	// Load matrix from file
	printf("Loading matrix \"%s\"\n", filename);
	sym = csc_load_matrix(filename, &m);
	
	if(sym) printf("Matrix is symmetric\n");
	else	printf("Matrix is general (non-symmetric)\n");
	
	// Print Matrix data
	//printf("CSC matrix data:\n");
	//csc_print_matrix(&m);
	
	// Initialize vectors for multiplication
	x = (double *)malloc(m.cols * sizeof(double));
	y = (double *)malloc(m.rows * sizeof(double));
	for(i = 0; i < m.cols; i++){
		x[i] = i+1;
	}
	//print_vector("\nx= ", x, m.cols);

	printf("Calculating matrix-vector product for %u iterations...\n", iterations);

	if(sym){
		// Time matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			csc_mvp_sym(&m,x,y);
		}
		clock_stop(&clock);	
	}else{
		// Time matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			csc_mvp(&m,x,y);
		}
		clock_stop(&clock);
	}
	printf("CSC mvp time: %.2fus \n\n", clock.sec);
	//print_vector("y= \n", y, m.rows);
	
	// Free resources
	free(x);
	free(y);
	csc_free_matrix(&m);
	
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
