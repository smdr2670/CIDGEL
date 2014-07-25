/*
** toric.c                                 Birk Huber, 4/99 
** -- support for computing first grobner bases of toric ideals.
** -- includes functions gset_rlex_rgb() and gset_compute_colon()
**    which once tested and working will move back to the gset.c file. 
**
**
** TiGERS,  Toric Groebner Basis Enumeration by Reverse Search 
** copyright (c) 1999  Birk Huber
**
*/
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "gset.h"
#include "matrices.h"
#include "Rsimp.h"

FILE *outfile;



/*
** void gset_compute_colon(gset g,int lv):
**    Compute deg rev lex rgb (lv least) for colon of g with variable lv.
**    - set ring_lv=lv.
**    - convert g to rgb wrt deg_rev_lex as implemented in binomial.h.
**    - divide each binomial by highest power of var lv possible.
**    - autoreduce.
**
*/
#define min(m,n) (((m)<(n))?(m):(n))



/**
 * @brief gset_code_ideal Given an mxn integer matrix M -- compute an rgb for the code ideal I_C.
 *        Expects generator matrix in standard form, reads off equations for a generating system, then
 *        uses Buchberger algorithm.
 * @param M input Matrix
 * @param m number of rows
 * @param n number of columns
 */
gset gset_code_ideal(int **M,int m, int n){
    int i;
    int j;
    int jk;
    binomial b;
    int **new_imatrix();
    gset g;

    /* set ring dimension*/
    if (ring_N!=n){
        fprintf(stderr,"ERROR gset_toric_ideal():");
        fprintf(stderr," matrix dimensions incompatible with ring\n");
        return 0;
    }



    

    /*
    fprintf(outfile,"Debug Matrixes, Corank: %d \n", crk);
    print_imatrix(outfile,debug_str,S,n,m);
    fprintf(outfile,"\n");
    print_imatrix(outfile,debug_str,V,n,n);
   */

    /* create new gset and read off equations corresponding to rows of V*/
    /*
    g=gset_new();
    for(i=crk;i<n;i++){
        b=binomial_new();
        for(j=0;j<n;j++){
            jk=IMref(V,i,j);
            if (jk>0){
                binomial_lead(b)[j]=jk;
                binomial_trail(b)[j]=0;
            }else {
                binomial_trail(b)[j]=-1*jk;
                binomial_lead(b)[j]=0;
            }
        }
        gset_insert(g,b);
    }
    */

    // read equations of generator matrix */

    g = gset_new();
    for(i=0;i<m;i++){
        b = binomial_new();
        for(j=0;j<n;j++){
            jk = IMref(M,i,j);
            if(j <= m-1){
                binomial_lead(b)[j] = jk;
                binomial_trail(b)[j] = 0;
            }else{
                binomial_lead(b)[j] = 0;
                binomial_trail(b)[j] = (code_dim-jk)%code_dim;
            }
            gset_insert(g,b);
        }
    }

    /* Adding non prime-part */
    gset_add_nonprime(g);


    //fprintf(outfile,"debug of Gset \n ");
    //gset_print(outfile, g);


    /* compute grobner basis wrt to this grade lex term order */
    // FIXME: lex basis can read off of the matrix with standard-form
    //        With exhaustive search, lex basis is not good enough
    //if(degree_comp == TRUE){
    //gset_rgb(g,monomial_grlexcomp);
    //}

    /* result is a grobner basis, make sure to reduce it to make it an rgb */
    gset_autoreduce(g);



    /* compute successive colon ideals */
    /*
    for(i=0;i<n;i++){
        fprintf(outfile,"debug compute colon : ");
        gset_print(outfile, g);
        gset_compute_colon(g,i);
        fprintf(outfile,"\n\n");

    }
    */

    /*unset weighting */
    for(i=0;i<n;i++) {
        ring_weight[i]=1;
    }

    /* free space used by matrices */
    free_imatrix(M);
 


    /* make sure result is a reduced gradedlex gbasis */
    // FIXME : lex basis  standard form
    //gset_rgb(g,monomial_grlexcomp);

    return g;
}












