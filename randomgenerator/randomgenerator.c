#include <stdio.h>
#include <stdlib.h>




int main(){

    srand(time(NULL));

    int i;
    int j;
    int result;
	char *ofname=0;
	FILE *outfile;
	int n;
	int k;
	int dim;
	char name[80];
	printf("This program generates a random Generator matrix for a linear [n,k] Code.\n");
	printf("Please enter n:");
	scanf("%d",&n);
	printf("Please enter k:");
	scanf("%d",&k);
	if(n <= k){
		printf("Error! n must be greater than k");
		exit(1);
	}
	printf("Please enter dimension:");
	scanf("%d",&dim);
	printf("Now the outputname: ");
	scanf("%s",name);

	/* open outfile */
    if (name!=0 && (outfile=fopen(name,"w"))==0){
        fprintf(stderr,"  couldn't open %s for output\n",ofname);
        exit(1);
    }

    fprintf(outfile, "M: { %d %d %d :\n",k,n,dim);
    for(i=0;i<k;i++){
    	for(j=0;j<n;j++){
    		if(k > j){
    			if(i == j){
    				fprintf(outfile, "1 ");
    			}else{
    				fprintf(outfile, "0 ");
    			}
    		}else{
                
                result = (rand()%(dim));
    			fprintf(outfile, "%d ", result );
    		}
    	}
    	fprintf(outfile,"\n");
    }
    fprintf(outfile, "}" );
}