#include <stdio.h>
#include <stdlib.h>

#include "coo.h"

int main(int argc, char *argv[]){
	unsigned int j,col,dfound;
	COO_Matrix coo;
	
	if (argc < 2) {
		printf("Usage: %s input_matrix\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	
	coo_init_matrix(&coo);
	coo_load_matrix(argv[1], &coo);
	
	j = 0;
	for(col = 0; col < coo.cols; col++){
		dfound = 0;
		while(j < coo.nz && coo.col[j] == col){
			if(coo.col[j] == coo.row[j]){
				dfound = 1;
			}
			j++;
		}
		if(dfound == 0) break;
	}
	
	coo_free_matrix(&coo);
	
	if(dfound == 1) printf("Matrix %s has a complete diagonal.\n", argv[1]);
	else printf("Matrix %s does not have a complete diagonal.\n", argv[1]);
	
	exit(EXIT_SUCCESS);
}
