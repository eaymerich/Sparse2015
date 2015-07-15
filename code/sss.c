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
#include "sss.h"

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);

int main(int argc, char *argv[]){
	unsigned int i;
	int sym;
	SSS_Matrix m;
	Clock clock;
	char* filename;
	double *x;
	double *y;
	
	parse_args(argc, argv);
	filename = argv[1];
	
	// Initialize matrices.
	sss_init_matrix(&m);
	
	// Load matrix from file
	//printf("Loading matrix \"%s\"\n", filename);
	sym = sss_load_matrix(filename, &m); // If matrix is non-symmetric, this call fails.
	
	//printf("Matrix is symmetric\n");
	
	// Print Matrix data
	//printf("SSS matrix data:\n");
	//sss_print_matrix(&m);
	
	// Initialize vectors for multiplication
	x = (double *)malloc(m.cols * sizeof(double));
	y = (double *)malloc(m.rows * sizeof(double));
	for(i = 0; i < m.cols; i++){
		x[i] = i+1;
	}
	//print_vector("\nx= ", x, m.cols);

	// Time matrix-vector product
	//printf("Calculating matrix-vector product for %u iterations...\n", iterations);
	clock_start(&clock);
	for(i=0; i < iterations; i++){
		sss_mvp2_sym(&m,x,y);
	}
	clock_stop(&clock);	
	//printf("SSS mvp time: %.2fus \n\n", clock.sec);
	//print_vector("y= \n", y, m.rows);
	
	// Print minimal info for testing purposes
	printf("%s\t%.2f\n", filename, clock.sec);
	
	// Free resources
	free(x);
	free(y);
	sss_free_matrix(&m);
	
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
