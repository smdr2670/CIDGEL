/*
** exsearch.h                                 Daniel Rembold 
** -- header file for the exhaustive search
**  
**
**
*/


int flip_ex_condition(binomial b,int degree_comp);
void insert(gset g1, gset *L);
int gset_equiv(gset g1, gset g2);
gset find(gset g1, gset L);
int exsearch(gset g1,int degree_comp, int no_print);