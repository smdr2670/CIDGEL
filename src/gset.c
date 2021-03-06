/** 
*   @file gset.c
*   @brief definitions and basic operations on gsets
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
#include "matrices.h"
#include "Rsimp.h"


gset gset_new(){
    gset g=0;
    if ((g=(gset)malloc(sizeof(struct gset_tag)))==0){
        fprintf(stderr,"gset_new(): malloc failed\n");
        return 0;
    }
    gset_first(g)=0;
    gset_cache_vtx(g)=0;
    gset_cache_edge(g)=0;
    return g;
}


void gset_free(gset g){
    binomial b;
    while((b=gset_first(g))!=0){
        gset_first(g)=binomial_next(b);
        binomial_free(b);
    }
    free((void *)g);
}


void gset_eliminatezero(gset g){
    binomial ptr;
    ptr = gset_first(g);
    for(ptr = gset_first(g); ptr != 0 ; ptr = binomial_next(ptr) ){
        if( monomial_stddegree(binomial_lead(ptr)) == 0 && monomial_stddegree(binomial_trail(ptr)) == 0){
            gset_delete_binomial(g,ptr);
        }
    }
}


void gset_delete_binomial(gset g, binomial b){
    binomial tmp;
    binomial prev;

    tmp = gset_first(g);

    if(0 == binomial_equal(b,tmp) ){
        binomial tmp2 = binomial_next(tmp);
        gset_first(g) = tmp2;
        binomial_free(tmp);
        return;
    }

    
    for(prev = gset_first(g); prev!=0 ; prev= binomial_next(prev)){
        tmp = binomial_next(prev);
        if( 0 == binomial_equal(b,tmp) ){
            binomial_next(prev)=binomial_next(tmp);
            binomial_free(tmp);
            return;
        }
    }
}


int gset_read(FILE *is,gset g){
    char c;
    binomial b;

    eatwhite(is);
    c=getc(is);
    if (c!='{') return FALSE;
    c=',';
    while(c==','){
        b=binomial_new();
        binomial_read(is,b);
        gset_insert(g,b);
        eatwhite(is);
        c=getc(is);
    }
    if (c!='}'){
        fprintf(stderr,"gset_read(): gset should end with a '}'\n");
        return FALSE;
    }
    ungetc(c,is);
    return TRUE;
}


void gset_print(FILE *of,gset g){
    binomial b;

    fprintf(of,"{");
    for(b=g->bottom;b!=0;b=b->next){
        binomial_print(of,b);
        if (b->next!=0){
            fprintf(of,", ");
        }
    }
    fprintf(of,"}");
}



void gset_init_print(FILE *of,gset g){
    binomial b;
    fprintf(of,"{");
    for(b=g->bottom;b!=0;b=b->next){
        print_monomial(of,binomial_lead(b));
        if (b->next!=0){
            fprintf(of,", ");
        }
    }
    fprintf(of,"}");
}


void gset_facet_print(FILE *of,gset g){
    binomial b;
    int tog=0;
    fprintf(of,"{");
    for(b=g->bottom;b!=0;b=b->next){
        if (b->ff==FACET){
            if (tog==1){
                fprintf(of,",");
            }else{
                tog=1;
            }
            binomial_print(of,b);
        }
    }
    fprintf(of,"}");
}



void gset_insert(gset g, binomial b){
    binomial t1,t2;
    if (gset_first(g)==0 || binomial_compair(gset_first(g),b)==FALSE){
        binomial_next(b)=gset_first(g);
        gset_first(g)=b;
    }else{
        t1=gset_first(g);
        while(t1!=0){
            t2=binomial_next(t1);
            if (t2==0 || binomial_compair(t2,b)==FALSE){
                binomial_next(t1)=b;
                binomial_next(b)=t2;
                break;
            }
            t1=t2;
        }
    }
}




binomial gset_downedge_lex(gset g){       
    binomial ptr=0;
    for(ptr=gset_first(g);ptr!=0;ptr=binomial_next(ptr)){
        if (binomial_lexordered(ptr)==FALSE && gset_isfacet(g,ptr)==TRUE ){
            break;
        }
    }
    return ptr;
}

binomial gset_downedge_grlex(gset g){       
    binomial ptr=0;
    for(ptr=gset_first(g);ptr!=0;ptr=binomial_next(ptr)){
        if (binomial_grlexordered(ptr)==FALSE && gset_isfacet(g,ptr)==TRUE ){
            break;
        }
    }
    return ptr;
}



binomial bmleadreduce(binomial m1, binomial B, binomial mlist){
    binomial res,ptr;
    /* case 1: the binomial is used in the decomp:*/
    if (monomial_divides(binomial_lead(B),binomial_lead(m1))==TRUE){
        res=binomial_new();
        binomial_copy(B,res);
        binomial_bumpto(res,m1);

        /* now reduce by lt(b) as many times as possible */
        while(monomial_divides(binomial_lead(B),binomial_trail(res))==TRUE){
            reducetrail(res,B);
        }

        /* now find a monomial dividing whats left and reduce by it*/
        for(ptr=mlist;ptr!=0;ptr=ptr->next){
            if (monomial_divides(binomial_lead(ptr),binomial_trail(res))==TRUE){
                reducetrail(res,ptr);
                return res;
            }
        }
    }else{ /* case 2 -- our monomial is a multiple of some monomial*/
        for(ptr=mlist;ptr!=0;ptr=ptr->next){
            if (monomial_divides(binomial_lead(ptr),binomial_lead(m1))==TRUE){
                res=binomial_new();
                binomial_copy(ptr,res);
                binomial_set(res);
                binomial_bumpto(res,m1);
                return res;
            }
        }
    }
    fprintf(stderr, "Warning BMreduce does not get to zero\n");
    return(binomial_new());
}




binomial find_divisor(binomial S,binomial L){
    binomial ptr;
    for(ptr=L;ptr!=0;ptr=ptr->next){
        if(ptr!=S&&monomial_divides(binomial_lead(ptr),binomial_lead(S))==TRUE){
            break;
        }
    }
    return ptr;
}



void remove_multiples(binomial S, binomial *L){
    binomial tlist=*L,ptr;
    *L=0;
    while(tlist!=0){
        ptr=tlist;  tlist=tlist->next; ptr->next=0;       /* pop tlist*/
        if (S==ptr || monomial_divides(binomial_lead(S),binomial_lead(ptr))==FALSE) {
            ptr->next=*L;  /* push L */
            *L=ptr;
        }
        else binomial_free(ptr);
    }
}



void bmrgb(binomial B, binomial *ms){
    //int i,jk;
    binomial DL;
    binomial SL;
    binomial S;
    binomial tmp;

    //binomial tlist;
    //binomial ptr;
    //int szero;

    DL=*ms;
    SL=0;

    while(DL!=0){
        tmp=DL; DL=DL->next;     /*  pop DL */
        tmp->next=SL; SL=tmp;    /*  push SL */

        /* calculate normal form of S pair with respect to {B}*/
        S=monomial_spair(B,binomial_lead(tmp));

        /* check if anything on DL or on SL divides S*/
        if (find_divisor(S,SL)!=0||find_divisor(S,DL)!=0)binomial_free(S);
        else {
            remove_multiples(S,&SL); /* remove any multiples of S from SL*/
            remove_multiples(S,&DL); /* remove any multiples of S from DL*/
            S->next=DL; DL=S;        /* finally push S onto DL */
        }
    }

    remove_multiples(B,&SL); /* now that we have a grobner basis
                               we can remove multiples of lead(B)
                               from the monomial list */

    // replace list bottom pointer
    *ms=SL;

}


gset gset_flip(gset g1, binomial b){
    binomial old=0,new=0,tmp,ptr,B,Bold;
    gset gres=0;
    //gset gtest=0;

    /* sainity check -- remove eventually*/
#ifdef GTEST
    if (lp_isfacet(g1,b)!=FACET){
        fprintf(stdout,"WARNING flip asked to flip non-facetbinomial\n");
    }
#endif

    /* copy all binomials of g's list ecept b to old and new lists*/
    for(ptr=gset_first(g1);ptr!=0;ptr=ptr->next){
        if (ptr!=b){
            tmp=binomial_new();
            binomial_copy(ptr,tmp);
            monomial_set(tmp);
            binomial_next(tmp)=new;
            new=tmp;
            tmp=binomial_new();
            binomial_copy(ptr,tmp);
            monomial_set(tmp);
            binomial_next(tmp)=old;
            old=tmp;
        }
    }
    Bold=binomial_new();
    binomial_copy(b,Bold);

    /* make flipped copy of b*/
    B=binomial_new();
    binomial_copy(b,B);
    binomial_flip(B);

    /* calculate the reduced GB of the binomial b and the monomials in M */
   
    bmrgb(B,&new);
    

    /* create new gset and insert flipped copy of b */
    gres=gset_new();
    gset_insert(gres,B);

   

    /* do lifting of elements in mlist and add them to new gset*/
    while((tmp=new)!=0){
        new=new->next;
        ptr=bmleadreduce(tmp,Bold,old);
        gset_insert(gres,ptr);
        binomial_free(tmp);
    }
   

    binomial_free(Bold);
    while(old!=0){
        ptr=old;
        old=old->next;
        binomial_free(ptr);
    }

    gset_autoreduce(gres);

    return gres;
}




void gset_autoreduce(gset g){
    binomial b1,b2;
    int changed;

    /* First remove redundant elements from the grobner basis*/
    while(find_divisor(gset_first(g),gset_first(g))!=0){
        gset_first(g)=(b1=gset_first(g))->next;
        binomial_free(b1);
    }
    b1=gset_first(g);
    while(b1!=0 && binomial_next(b1)!=0){
        if (find_divisor(binomial_next(b1),gset_first(g))!=0) {
            b2=binomial_next(b1);
            binomial_next(b1)=binomial_next(b2);
            binomial_free(b2);
        }
        else b1=binomial_next(b1);
    }

    /* now reduce trailing terms */
    b1=gset_first(g);
    while(b1!=0){
        changed=1;
        while(changed==1){
            changed=0;
            b2=gset_first(g);
            while(b2!=0){
                if (b2!=b1 &&
                        monomial_divides(binomial_lead(b2),binomial_trail(b1))==TRUE){
                    reducetrail(b1,b2);
                    changed=1;
                }
                b2=binomial_next(b2);
            }
        }
        b1=binomial_next(b1);
    }

    /* now reinsert terms to ensure list is in order */
    b1=gset_first(g);
    gset_first(g)=0;
    while((b2=b1)!=0) {
        b1=binomial_next(b1);
        gset_insert(g,b2);
    }

}




#define ordered(b) ((comp(binomial_lead(b),binomial_trail(b))>=0)?TRUE:FALSE)
int changed;

void gset_rgb(gset g, int (*comp)(monomial,monomial)){
    void process_pair(binomial,binomial,binomial,binomial,binomial *,int (*comp)());
    binomial ptr,b,DL=0,SL=0,NL=0;
    int passes=0;

    do{
        DL=0; SL=0; NL=0;
        changed=FALSE;

        /* reorient g and put binomial list on DL*/
        ptr=gset_first(g);
        while((b=ptr)!=0){
            ptr=binomial_next(ptr);
            if (ordered(b)==FALSE) binomial_flip(b);
            binomial_next(b)=DL;
            DL=b;
        }

        /*  */
        while(DL!=0){
            NL=0;
            b=DL; DL=DL->next;  /* pop DL*/
            b->next=SL; SL=b; /* push ptr onto SL*/
            for(ptr=SL;ptr!=0;ptr=binomial_next(ptr)) process_pair(b,ptr,DL,SL,&NL,comp);
            for(ptr=DL;ptr!=0;ptr=binomial_next(ptr)){
                process_pair(b,ptr,DL,SL,&NL,comp);
                if (binomial_next(ptr)==0){ binomial_next(ptr)=NL; NL=0;}
            }
        }

        /* reinsert elements in SL into gset */
        ptr=SL;
        gset_first(g)=0;
        while(ptr!=0){
            ptr=(b=ptr)->next;
            gset_insert(g,b);
        }


        gset_autoreduce(g);

        if (passes>1){
            if (changed==TRUE ){fprintf(stderr,"looping pass (%d)\n",passes);}
            else {fprintf(stderr,"stopping after %d passes\n",passes);}
        }
        passes++;
    }while(changed==TRUE);

}


void process_pair(binomial b1,binomial b2,binomial DL,binomial SL,binomial *NL,
                  int (*comp)(monomial,monomial)){
    binomial S,ptr,find_divisor(binomial,binomial);
    //int i,tmp;

    /* test is lead terms are rel prime [buchberger's first criterion]*/
    if (monomial_rel_prime(binomial_lead(b1),binomial_lead(b2))==TRUE)return;

    S=binomial_new();
    binomial_copy(b1,S);
    if (binomial_spair(S,b2)==TRUE){
        binomial_free(S);
        return;
    }
    if (ordered(S)==FALSE) binomial_flip(S);
    while(((ptr=find_divisor(S,DL))!=0)||
          ((ptr=find_divisor(S,SL))!=0)||
          ((ptr=find_divisor(S,*NL))!=0)){
        if (binomial_spair(S,ptr)==TRUE){
            binomial_free(S);
            return;
        }
        if (ordered(S)==FALSE) binomial_flip(S);
    }
    changed=TRUE;
    binomial_next(S)=*NL;
    *NL=S;
}


#define max(a,b) (((a)>(b))?(a):(b))


void gset_setfacets(gset g){
    int debug = 0;
    binomial b;
    monomial ml,mt;
    int i,d=0,dl,dt;
    gset_nfacets(g)=0;
    gset_nelts(g)=0;
    for(b=gset_first(g);b!=0;b=binomial_next(b)){
        //fprintf(stderr, "%d ", debug );
        debug++;
        if (gset_isfacet(g,b)==FACET){
         gset_nfacets(g)++;
        }
        gset_nelts(g)++;
        ml=binomial_lead(b); dl=0;
        mt=binomial_trail(b); dt=0;
        for(i=0;i<ring_N;i++){
            dl+=ml[i]; dt+=mt[i];
        }
        d=max(d,max(dl,dt));
    }
    gset_deg(g)=d;
}


int lptest=1;   /* 1 test only LP, 2 test only flipability, 3 test both */
int gset_isflippable(gset g,binomial b); 

int gset_isfacet(gset g,binomial b){
    
    int rval=UNKNOWN;
    rval=b->ff;

    /* This cost */
    if(monomial_stddegree(binomial_lead(b))== 2 && monomial_stddegree(binomial_trail(b))==0 ){
        rval = NONFACET;
    }

    /* No considerable speedup reached */
    //if(binomial_degree_compatible(b) != 0 ){
    //    rval = NONFACET;
    //}else{
    //    rval = FACET;
    //}


    if (rval!=FACET && rval!=NONFACET){
        switch(lptest){
        // use linear programming
        case 1: rval=(b->ff=lp_isfacet(g,b));
            break;
        // use flippability test
        case 2: if (gset_isflippable(g,b)==FALSE) rval=(b->ff=NONFACET);
            else rval=(b->ff=FACET);
            break;
        // use both
        case 3: if (gset_isflippable(g,b)==FALSE) rval=(b->ff=NONFACET);
            else rval=(b->ff=lp_isfacet(g,b));
            break;
        default: fprintf(stderr,"Warning lptest has invalid value in gset_isflippable()\n");
            exit(4);
        }
    }

    return rval;
}


int gset_isflippable(gset g,binomial b){
    int ctg=0,ctn=1,val=TRUE;
    binomial new=0,tmp,ptrg,ptrn,btmp;
    /* copy all binomials of g's list except b to old and new lists*/
    for(ptrg=gset_first(g);ptrg!=0;ptrg=ptrg->next){
        ctg++;
        if (ptrg!=b){
            tmp=binomial_new();
            binomial_copy(ptrg,tmp);
            monomial_set(tmp);
            binomial_next(tmp)=new;
            new=tmp;
        }
    }
    btmp=binomial_new();
    binomial_copy(b,btmp);


    /* calculate the reduced GB of the binomial b and the monomials in M */
    bmrgb(btmp,&new);

    /* test if new and old have same leading terms */
    for(ptrn=new           ;ptrn!=0;ptrn=ptrn->next) ctn++;
    if (ptrn!=ptrg) val=FALSE;  /* initial ideals have diferent numbers of generators */
    else {
        for(ptrn=new;ptrn!=0;ptrn=ptrn->next){
            for(ptrg=gset_first(g);ptrg!=0;ptrg=ptrg->next){
                if (monomial_equal(binomial_lead(ptrn),binomial_lead(ptrg))==TRUE)break;
            }
            if (ptrg==0)break;  /* no match found */
        }
        if (ptrn!=0)val=FALSE; /* not all terms had matches found */
    }

    /* clean up */
    binomial_free(btmp);
    while(new!=0){
        ptrn=new;
        new=new->next;
        binomial_free(ptrn);
    }

    return val;
}


int lp_isfacet(gset g,binomial b){
    binomial ptr;
    int D=ring_N,M=0;
    int idx=0,i;
    //double dtmp0,dtmp;

    /* need cheeper way to determine the size of g */
    for(ptr=gset_first(g);ptr!=0;ptr=binomial_next(ptr)){
        if (ptr->ff!=NONFACET){
            M++;
        }
    }

    /* make sure the lp structures have enough space for the program*/
    LP_get_space(M,D+2*M);

    /*  Set up linear program. */
    ptr=gset_first(g);
    idx=0;
    while(ptr!=0){
        if (ptr->ff!=NONFACET){
            for(i=0;i<D;i++){  /* make sure supports of trail and lead disjoint? */
                LP_A(idx,i)=binomial_lead(ptr)[i]-binomial_trail(ptr)[i];
            }
            for(i=0;i<M;i++){
                LP_A(idx,D+i)=0;
                LP_A(idx,D+M+i)=0;
            }
            if (ptr==b) {for(i=0;i<D;i++) LP_A(idx,i)*=-1;}
            LP_A(idx,D+M+idx)=1;
            LP_A(idx,D+idx)=-1;
            LP_B[idx]=1;
            idx++;
        }
        ptr=ptr->next;
    }
    for(i=0;i<D+M;i++){
        LP_C[i]=0;
        LP_X[i]=0;
        LP_NonBasis[i]=i;
    }
    for(i=0;i<M;i++){
        LP_Basis[i]=D+M+i;
        LP_X[D+M+i]=1;
        LP_C[D+M+i]=1;
    }

    /* Run simplex algorithm */
    Rsimp(LP_M,LP_N,LP_A,LP_B,LP_C,LP_X,LP_Basis,LP_NonBasis,
          LP_R,LP_Q,LP_t1,LP_t2);

    /* use results to determine if the b is a facet*/
    b->ff=FACET;
    for(i=0;i<M;i++) if (LP_Basis[i]>=M+D) {b->ff=NONFACET;break;}

    return b->ff;
}


void gset_add_nonprime(gset g){
    binomial ptr;
    int *tmp;
    int i;
    tmp= malloc(ring_N*sizeof(int));
    /*check first term of binomial */
    for(ptr=gset_first(g);ptr!=0;ptr = binomial_next(ptr)){
        if(binomial_first_term_degree(ptr)==1){
            //write function which returns the position, which shall NOT be added
            *(tmp+binomial_variable_position(ptr)) = 1;
        }
    }

    /* all non-prime binomials which have to be added are now calculated */

    for(i=0;i<ring_N;i++){
        /* non-prime binomal has to be added */
        if( *(tmp+i) != 1) {
            binomial b = binomial_new();
            *((b->exps1)+i) = code_dim;
            gset_insert(g,b);
        }
    }

}


int gset_only_degreecompatible(gset g){
    binomial ptr;
    for(ptr=gset_first(g);ptr!=0;ptr = binomial_next(ptr)){
        if( (binomial_degree_compatible(ptr)) != 1){

            return 0;
        }
    }
    return 1;

}


int countDC(gset g1){
    int d=0;
    binomial ptr;
    for(ptr = gset_first(g1); ptr!= 0 ; ptr =binomial_next(ptr) ){
        if(binomial_degree_compatible(ptr)== TRUE && ptr->ff==FACET){
            d++;
        }
    }

    return d;
}


void gset_copy(gset src, gset dest){
    binomial ptr;
    for(ptr = gset_first(src); ptr!=0; ptr = binomial_next(ptr)){
        binomial ptr2 = binomial_new();
        binomial_copy(ptr,ptr2);
        gset_insert(dest,ptr2);
        
    }
}




void gset_puncture(gset G1, int pos){
    binomial ptr;
    for(ptr = gset_first(G1); ptr!= NULL ; ptr=binomial_next(ptr) ){
        binomial_puncture(ptr,pos);
    }
}


