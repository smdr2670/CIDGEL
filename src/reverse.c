/** 
*   @file reverse.c
*   @brief Reverse search loop and main calling program for CIDGEL.
*          
*   
*   This file is reused from TiGERS,
*   Toric Groebner Basis Enumeration by Reverse Search 
*   copyright (c) 1999  Birk Huber
*
*   @author Birk Huber, 4/99 
*   @author Daniel Rembold
*   @bug No known bugs
*
*/

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "gset.h"
#include "reverse.h"
#include "match.h"

int rsearch_cache=TRUE;
int degree_comp=FALSE;



/*
** Output initial ideals and numbers of facets at each stage.
*/
FILE *outfile;
int print_init=TRUE;
int print_tree=TRUE;

int stats_minfacets=-1;
int stats_minelts=-1;
int stats_mindeg=-1;
int stats_maxfacets=0;
int stats_maxelts=0;
int stats_maxdeg=0;
int stats_tdepth=0;
int stats_ecounter=0;

int match_fan = 0;
//int no_print = 0;

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

struct node* list = NULL;
struct node* list2 = NULL;



void vertex_print(gset g1,int no_print ){



    stats_maxfacets=max(gset_nfacets(g1),stats_maxfacets);
    stats_maxelts=max(gset_nelts(g1),stats_maxelts);
    stats_maxdeg=max(gset_deg(g1),stats_maxdeg);

    if (stats_minfacets<0){
        stats_minfacets=gset_nfacets(g1);
    }else {
        stats_minfacets=min(stats_minfacets,gset_nfacets(g1));
    }

    if (stats_minelts<0) {
        stats_minelts=gset_nelts(g1);
    }else{
        stats_minelts=min(stats_minelts,gset_nelts(g1));
    }

    if (stats_mindeg<0){
        stats_mindeg=gset_deg(g1);
    }else {
        stats_mindeg=min(stats_mindeg,gset_deg(g1));
    }

    if(no_print == TRUE){
        return;
    }

    fprintf(outfile,"Vtx: %d (%d facets/%d among them degree compatible /%d binomials/degree %d)\n",
            gset_id(g1),gset_nfacets(g1),countDC(g1), gset_nelts(g1),gset_deg(g1));

    if (print_init==TRUE){
        fprintf(outfile,"    Initial ideal:");gset_init_print(outfile,g1);
        fprintf(outfile,"\n  Facet Binomials:");gset_facet_print(outfile,g1);
    } else{
        gset_print(outfile,g1);
    }
    if (print_tree==TRUE && rsearch_cache==TRUE && gset_cache_vtx(g1)!=0){
        fprintf(outfile,"\n  Reducing Edge: [%d,%d]",
                gset_id(g1),gset_id(gset_cache_vtx(g1)));
    }
    fprintf(outfile,"\n\n");
}



int flip_condition_tig(binomial b){
    if(degree_comp == TRUE){
        return (binomial_degree_compatible(b) == 0 && binomial_grlexordered(b)==TRUE);
    }else{
        return (binomial_lexordered(b)==TRUE);
    }
}



int rsearch(gset g1, int number, int no_print){

    gset G1=0;
    gset G2=0;
    binomial b=0;
    binomial btmp;
    int counter=0,done=FALSE, depth=0;

    stats_tdepth=0;
    stats_ecounter=0;

    struct node *new_node;

    if(DEBUG){
    fprintf(stderr,"CHECKPOINT 1\n");
    }
    /* make copy of first grobner basis */
    G1=gset_new();
    for(b=gset_first(g1);b!=0;b=binomial_next(b)){
        btmp=binomial_new();
        binomial_copy(b,btmp);
        gset_insert(G1,btmp);
    }

    if(DEBUG){
    fprintf(stderr,"CHECKPOINT 2\n");
    }



    /* Do groebner walk to grobner basis wrt default term order*/

    if(degree_comp == TRUE){
        while((b=gset_downedge_grlex(G1))!=0){
            fprintf(stderr,"Warning: rsearch doing walk to get to root\n");
            G2=gset_flip(G1,b);
            gset_free(G1);
            G1=G2;
        }

    }else{
        while((b=gset_downedge_lex(G1))!=0){
        fprintf(stderr,"Warning: rsearch doing walk to get to root\n");
        G2=gset_flip(G1,b);
        gset_free(G1);
        G1=G2;
        }
    }

    fprintf(stderr,"doing walk finished\n");

    if(DEBUG){
    fprintf(stderr,"CHECKPOINT 3\n");
    }

    btmp=(b=gset_first(G1));

    /* output first groebner basis*/
    //FIXME testing warum lange laufzeit
    gset_setfacets(G1);

    if(DEBUG){
    fprintf(stderr,"CHECKPOINT 4\n");
    }
    gset_id(G1)=++counter;

    
    vertex_print(G1,no_print);

    if(match_fan == TRUE){
        new_node=newNode(gset_id(G1),gset_nfacets(G1),gset_nelts(G1),gset_deg(G1));
        if(number == 1){
            sortedInsert(&list, new_node);   
        }else{
            sortedInsert(&list2, new_node);
        }
    }

     if(DEBUG){
    fprintf(stderr,"CHECKPOINT 5\n");
    }
    /* exits directly if Groebner Basis is the only degree compatible */
    if(gset_only_degreecompatible(G1) && degree_comp == TRUE){
        fprintf(outfile, "\nThe only degree compatible Groebner Basis found!\n");
        vertex_print(G1,no_print);

        if(match_fan == TRUE){
            new_node=newNode(gset_id(G1),gset_nfacets(G1),gset_nelts(G1),gset_deg(G1));
            if(number == 1){
                sortedInsert(&list, new_node);   
            }else{
                sortedInsert(&list2, new_node);
            }
        }
        gset_free(G1);
        return counter;
    }

     if(DEBUG){
    fprintf(stderr, "entering while loop" );
    }
    
    while(done==FALSE){

        /* go through all binomials of G1 */
        while(b!=0){

            //FIXME : lange laufzeit in den Griff bekommen
            //if (binomial_grlexordered(b)==TRUE && flip_condition_tig(b)   ){
            if (gset_isfacet(G1,b)==TRUE && flip_condition_tig(b) ){
                stats_ecounter++;

                G2=gset_flip(G1,b);
                 if(DEBUG){
                    fprintf(stderr,"CHECKPOINT 6\n");
                 }   
                if(degree_comp == TRUE){
                    btmp=gset_downedge_grlex(G2);
                }else{
                    btmp=gset_downedge_lex(G2);
                }
                

                /* check if G1 and G2 are adjacent in the reverse search tree*/
                if (monomial_equal(binomial_lead(b),binomial_trail(btmp))==TRUE &&
                        monomial_equal(binomial_trail(b),binomial_lead(btmp))==TRUE){

                    if (rsearch_cache==TRUE){
                        gset_cache_vtx(G2)=G1;
                        gset_cache_edge(G2)=b->next;
                    }
                    else{
                        gset_free(G1);
                    }

                    G1=G2;
                    depth++;
                    if (stats_tdepth<depth){
                        stats_tdepth=depth;
                    }
                    b=gset_first(G1);
                    if(DEBUG){
                        fprintf(stderr,"CHECKPOINT 7\n");
                    }
                    
                    gset_setfacets(G1);
                    if(DEBUG){
                         fprintf(stderr,"CHECKPOINT 8\n");
                    }               
                    
                    gset_id(G1)=++counter;
                    if(counter%1000 == 0){
                        fprintf(stderr,"%d Groebner bases counted yet...\n" , counter);
                    }
                    vertex_print(G1,no_print);

                    /*Sets the linked list for matching 2 Groebner fans */
                    if(match_fan == TRUE){
                        new_node=newNode(gset_id(G1),gset_nfacets(G1),gset_nelts(G1),gset_deg(G1));
                        if(number == 1){
                            sortedInsert(&list, new_node);   
                        }else{
                            sortedInsert(&list2, new_node);
                        }
                    }
                }
                else {
                    gset_free(G2);
                    b=binomial_next(b);
                }
            }
            else b=binomial_next(b);
        }

        depth--;
        if (rsearch_cache==TRUE){
            b=gset_cache_edge(G1);
            if ((G2=gset_cache_vtx(G1))==0){
                done=TRUE;
            }else{
                gset_free(G1);
                G1=G2;
            }
        }else{
            if(degree_comp == TRUE){

                
                if ((btmp=gset_downedge_grlex(G1))==0){
                     done=TRUE;
                }else{
                    G2=gset_flip(G1,btmp);
                    b=gset_first(G2);
                    while(monomial_equal(binomial_lead(b),binomial_trail(btmp))==FALSE){
                        b=binomial_next(b);
                    }
                    b=binomial_next(b);
                    gset_free(G1);
                    G1=G2;
                }
            }else{

                if ((btmp=gset_downedge_lex(G1))==0){
                     done=TRUE;
                }else{
                    G2=gset_flip(G1,btmp);
                    b=gset_first(G2);
                    while(monomial_equal(binomial_lead(b),binomial_trail(btmp))==FALSE){
                        b=binomial_next(b);
                    }
                    b=binomial_next(b);
                    gset_free(G1);
                    G1=G2;
                }

            }
        }
    }
    gset_free(G1);

    if(match_fan == TRUE && number == 2){

        if(0==match(&list,&list2)){
           fprintf(outfile,"\nGroebnerfan of the Codes are compatible\n");
         }else{
           fprintf(outfile,"\nGroebnerfan of the Codes are NOT compatible\n");
         }
    }  
    return counter;
}











