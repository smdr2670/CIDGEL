#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void checkCols(int *mat,int rows, int cols, int dim);
void checkRows(int *mat,int rows, int cols, int dim);
void printMatrix(char *ofname,int *mat,int rows, int cols, int dim );


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
 

    int i;
    int j;
    int k;

    int zerocounter;
    int border;
    int counter;

    // Matrix as 1 dimensional array
    int *mat = (int *)malloc(rows * cols * sizeof(int));

    // now mat[offset] corresponds to m(i, j)
    int offset;
    



    /* Beschreibe Matrix */
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

    checkCols(mat,rows,cols,dim);
    checkRows(mat,rows,cols,dim);



    
    /* Print matrix */
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            offset = i * cols + j;
            printf("%d", mat[offset]);
        }
        printf("\n");
    }
    printf("\n");

    printMatrix(ofname,mat,rows,cols, dim);


    free(mat);
}


void printMatrix(char *ofname,int *mat,int rows, int cols, int dim ){
    FILE *outfile;
    int offset,i,j;
    /* open outfile */
    if (ofname!=0 && (outfile=fopen(ofname,"w"))==0){
        fprintf(stderr,"  couldn't open %s for output\n",ofname);
        exit(1);
    }

    fprintf(outfile, "M: { %d %d %d :\n", rows, cols,dim );

        /* Print matrix */
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            offset = i * cols + j;
            fprintf(outfile, "%d ", mat[offset]);
        }
        fprintf(outfile,"\n");
    }

    fprintf(outfile, "}" );
}

    




void checkCols(int *mat,int rows, int cols, int dim){
    srand(time(NULL));
    int k,border,zerocounter,offset,i,j;
    /* checks if the column has a nonzero column vector */
    for(j=rows;j<cols;j++ ){
        zerocounter = 0;
        for(i=0;i<rows;i++){
            offset = i * cols + j;
            if(mat[offset]==0){
                zerocounter++;
            }
        }
        if(zerocounter == rows){
            border = ((rand()%rows)+1)%rows;
            for(k=0;k<border;k++){
                offset = k * cols + j;
                mat[offset] = dim-1;
            }
        }
    }

}

void checkRows(int *mat,int rows, int cols, int dim){
    srand(time(NULL));
    int k,border,zerocounter,offset,i,j;
        /* Check if A has a nonzero row vector*/
    for(i = 0; i<rows ; i++){
        zerocounter = 0;
       for(j=rows;j<cols;j++){
            offset = i * cols + j;
            if(mat[offset]== 0){
                zerocounter++;
            }
       }
       // TODO zufalls zahl neu generieren
       if(zerocounter == cols-rows){
            border = (rand()%(zerocounter)+1) + rows;
            printf("row : %d  border %d \n",i, border );
            for(k=(rows);k<border;k++){
                offset = i * cols +k;
                mat[offset] = dim-1;
            }

       }
    }


}