/*
** binomial.h                                 Birk Huber, 4/99 
** -- header file for binomal type and main operations on binomials
**  
**
** TiGERS,  Toric Groebner Basis Enumeration by Reverse Search 
** copyright (c) 1999  Birk Huber
**
** Commented and slightly edited by Daniel Rembold
**
*/



#ifndef BINOMIAL_H 
extern int ring_N;
extern int ring_lv;
extern int *ring_weight;
extern int code_dim;
#endif

void code_dim_set(int n);
int ring_set(int n);
int ring_getvar(FILE *ifile);
#define ring_putvar(ofile,v) fprintf(ofile,"%c",'a'+v)


typedef int *monomial;
int monomial_divides(monomial m1, monomial m2);
int monomial_equal(monomial m1, monomial m2);
int monomial_lexcomp(monomial m1, monomial m2);
int monomial_rlexcomp(monomial m1, monomial m2);
int monomial_grlexcomp(monomial m1, monomial m2);
int monomial_rel_prime(monomial m1, monomial m2);
void print_monomial(FILE *of, int *exps);

typedef struct bin_tag *binomial;
struct bin_tag{
    int *exps1;
    int *exps2;
    int *E;
    int ff;
    int bf;
    binomial next;
};


#define BINOMIAL 1
#define MONOMIAL 0
#define FACET    0
#define NONFACET 1
#define UNKNOWN -1

#define binomial_next(b) (b->next)
#define binomial_lead(b) (b->exps1)
#define binomial_trail(b) (b->exps2)
#define monomial_set(b)   (b->bf=MONOMIAL)
#define binomial_set(b)   (b->bf=BINOMIAL)
#define binomial_facet(b) (b->ff)


binomial binomial_new();
void binomial_free(binomial m);
void binomial_read(FILE *is, binomial b);
void binomial_print(FILE *of, binomial b);
binomial monomial_spair(binomial b,monomial m);
void binomial_flip(binomial b);
int binomial_spair(binomial b1, binomial b2);
int binomial_first_term_degree(binomial b);
int binomial_degree_compatible(binomial b);
void binomial_copy(binomial src,binomial dest);
int binomial_variable_position(binomial b);
void reducetrail(binomial b1, binomial b2);
void binomial_bumpto(binomial b1, binomial b2);
int binomial_compair(binomial b1,binomial b2);
int ring_read(FILE *infile);


/*
** return TRUE if binomials marking agrees with respect to the
** lexicographic term order
*/
#define binomial_lexordered(b)\
    ((monomial_lexcomp(binomial_lead(b),binomial_trail(b))>=0) ? TRUE : FALSE)

/*
** return TRUE if binomials marking agrees with respect to the
** graded lexicograpic term order
*/
#define binomial_grlexordered(b)\
    ((monomial_grlexcomp(binomial_lead(b),binomial_trail(b))>=0) ? TRUE : FALSE)
