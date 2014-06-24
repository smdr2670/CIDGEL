#include <stdio.h>
#include <stdlib.h>
#include <time.h>




int main(int argc, char* argv[]){

    srand(time(NULL));

    int rows ;
    int cols ;
    int dim ;
    char *ofname;

    sscanf (argv[1], "%i", &cols);
    sscanf (argv[2], "%i", &rows);
    sscanf (argv[3], "%i", &dim);
    ofname = argv[4];
 

    /*
    if(argc != 4){
        fprintf(stderr, "l2print");
    }
    */
    
    int i;
    int j;
    int k;

    int zerocounter;
    int border;
    int counter;

    // Matrix as 1 dimensional array
    int *mat = (int *)malloc(rows * cols * sizeof(int));

    int offset;
    // now mat[offset] corresponds to m(i, j)

	/* open outfile */
    /*
    if (name!=0 && (outfile=fopen(name,"w"))==0){
        fprintf(stderr,"  couldn't open %s for output\n",ofname);
        exit(1);
    }
    */

    //m(i,j)
    for(i=0;i<rows;i++){
    	for(j=0;j<cols;j++){
            offset = i * cols + j;
            /*Identity Matrix */
    		if(rows > j){
                if(i == j){
                    mat[offset] = 1;
    			}else{
                    mat[offset] = 0;
    			}
            /* Generate random */
    		}else{
                mat[offset] = (rand()%(dim));
            }
    	}
    }
    
    /* checks if the column has a nonzero column vector */
    for(j=rows+1;j<cols;j++ ){
        counter = 0;
        for(i=0;i<rows;i++){
            offset = i * cols + j;
            if(mat[offset]==0){
                zerocounter++;
            }
        }
        if(counter == rows){
            border = ((rand()%rows)+1)%rows;
            for(k=0;k<border;k++){
                offset = i * cols +k;
                mat[offset] = dim-1;
            }
        }
    }


    
    /* Print matrix */
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            offset = i * cols + j;
            printf("%d", mat[offset]);
        }
        printf("\n");
    }


    free(mat);
}


void printMatrix(){


}