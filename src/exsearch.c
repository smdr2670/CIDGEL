/*
** exsearch.c                                 Birk Huber, 4/99 
**   -- exhaustive search main loop.
*/
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "gset.h"
#include "exsearch.h"


/**
 * @brief insert use next pointer in Gsets to implement a simple linked list.
 *               for now stored in no particular order.
 * @param g1
 * @param L
 */
void insert(gset g1, gset *L){
    g1->next=*L;
    *L=g1;
}

/**
 * @brief gset_equiv determine if two gsets contain equivalent marked binomials
 * @param g1 first binomial
 * @param g2 second binomial
 * @return TRUE(=0) of FALSE(=1)
 */
int gset_equiv(gset g1, gset g2){
    binomial p1,p2;

    p1=gset_first(g1);
    p2=gset_first(g2);

    while(p1!=0 && p2!=0){
        if(monomial_equal(binomial_lead(p1),binomial_lead(p2))!=TRUE){
            return FALSE;    
        } 
        if(monomial_equal(binomial_trail(p1),binomial_trail(p2))!=TRUE) {
            return FALSE;
        }
        p1=binomial_next(p1);
        p2=binomial_next(p2);
    }
    if (p1!=p2) {
        return FALSE;
    }
    return TRUE;
}



/**
 * @brief find check a gset against those on list untill an equivalent one is found.
 * @param g1 gset which shall be found in list l
 * @param L list that will be searched after g1
 *
 *
 */
gset find(gset g1, gset L){
    gset g2=0;
    for(g2=L;g2!=0;g2=g2->next){
        if (gset_equiv(g1,g2)==TRUE){
            break;
        }
    }
    return g2;
}


/**
 * @brief flip_condition extra condition for degree compatible groebner base for facet binomial to be flipped
 * @param b given binomial
 * @return 1 for fullfilled condition, 0 for not fullfilled
 */
int flip_ex_condition(binomial b,int degree_comp){
    if(degree_comp == TRUE){
        return (binomial_degree_compatible(b) == 0 );
    }else{
        return 1;
    }
}

extern int stats_ecounter;

/**
 * @brief exsearch breadth first search to find all degree compatible groebner bases
 *                 of a code ideal
 * @param g1 starting groebner base
 * @return number of groebner bases found
 */
int exsearch(gset g1, int degree_comp,int no_print){
    gset G1;
    gset G2;
    gset todo=0;
    gset seen=0;
    //gset tmp;

    binomial b=0,btmp;
    int counter=0;
    stats_ecounter=0;

    /* make copy of first grobner basis */
    G1=gset_new();
    for(b=gset_first(g1);b!=0;b=binomial_next(b)){
        btmp=binomial_new();
        binomial_copy(b,btmp);
        gset_insert(G1,btmp);
    }

    insert(G1,&todo);
    gset_setfacets(G1);
    gset_id(G1)=++counter;
    vertex_print(G1,no_print);

    while(todo!=0){
        G1=todo; todo=todo->next;
        insert(G1,&seen);
        for(b=gset_first(G1);b!=0;b=binomial_next(b)){
            if (gset_isfacet(G1,b)==TRUE && binomial_grlexordered(b)==TRUE  && flip_condition(b,degree_comp)  ){
                G2=gset_flip(G1,b);
                stats_ecounter++;
                if (find(G2,todo)==0 && find(G2,seen)==0){
                    insert(G2,&todo);
                    gset_setfacets(G2);
                    gset_id(G2)=++counter;
                    vertex_print(G2,no_print);
                }
                else gset_free(G2);
            }
        }
    }

    while(seen!=0){
        G2=seen;
        seen=G2->next;
        gset_free(G2);
    }

    return counter;
}







