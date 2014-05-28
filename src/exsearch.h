/*
** exsearch.h                                 Daniel Rembold 
** -- header file for the exhaustive search
**  
**
**
*/

void insert(gset g1, gset *L);
int gset_equiv(gset g1, gset g2);
gset find(gset g1, gset L);
int exsearch(gset g1);