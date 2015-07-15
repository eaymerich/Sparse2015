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
#include "tjds.h"

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);

int main(int argc, char *argv[]){
	unsigned int i;
	int sym;
	TJDS_Matrix tjds;
	Clock clock;
	char* filename;
	double *x;
	double *y;
	
	parse_args(argc, argv);
	filename = argv[1];
	
	// Initialize matrices.
	tjds_init_matrix(&tjds);
	
	// Load matrix from file into COO.
	//printf("Loading matrix \"%s\"\n", filename);
	sym = tjds_load_matrix(filename, &tjds);
	
	//if(sym) printf("Matrix is symmetric\n");
	//else	printf("Matrix is general (non-symmetric)\n");
	
	// Print Matrix data
	//printf("TJDS matrix data:\n");
	//tjds_print_matrix(&tjds);
	
	// Initialize vectors for multiplication
	x = (double *)malloc(tjds.cols * sizeof(double));
	y = (double *)malloc(tjds.rows * sizeof(double));
	for(i = 0; i < tjds.cols; i++){
		x[i] = i+1;
	}
	//print_vector("\nx= ", x, coo.cols);
	
	// Permute X vector before doing TJDS product
	tjds_permute_vector(&tjds, x);

	//printf("Calculating matrix-vector product for %u iterations...\n", iterations);

	if(sym){
		// Time TJDS matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			tjds_mvpP_np_sym(&tjds,x,y);
		}
		clock_stop(&clock);
	}else{
		// Time TJDS matrix-vector product
		clock_start(&clock);
		for(i=0; i < iterations; i++){
			tjds_mvp_np(&tjds,x,y);
		}
		clock_stop(&clock);
	}
	//printf("TJDS mvp time: %.2fus \n\n", clock.sec);
	//tjds_permute_back_vector(&tjds, y);
	//print_vector("y= \n", y, tjds.rows);
	
	// Print minimal info for testing purposes
	printf("%s\t%.2f\n", filename, clock.sec);
	
	// Free resources
	free(x);
	free(y);
	tjds_free_matrix(&tjds);
	
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
