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

/**
 * @brief Determines if m1 divides m2.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_divides(monomial m1, monomial m2);

/**
 * @brief Determines if monomials m1 and m2 have same exponents.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_equal(monomial m1, monomial m2);


int monomial_lexcomp(monomial m1, monomial m2);


int monomial_rlexcomp(monomial m1, monomial m2);


int monomial_grlexcomp(monomial m1, monomial m2);

/**
 * @brief Determines if m1 is reletively prime to m2.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_rel_prime(monomial m1, monomial m2);

/**
 * @brief Prints monomial in a file.
 * @param of Outputfile where monomial will be printed.
 * @param exps Vector which contains the exponents.
 */
void print_monomial(FILE *of, int *exps);


/**
 * @brief Reads ascii representation of monomial as power product
 *         onvert to exponent vector stored in *exps.
 * @param is The input stream.
 * @param exps Store the exponent vector in exps.s
 *
 */
void get_monomial(FILE *is,int *exps);

/**
 * @brief monomial_lcm copy least common multiple of m1 and m2 into m3.
 * @param m1 first monomial
 * @param m2 second monomial
 * @param m3 lcm monomial
 * @return
 */
void monomial_lcm(monomial m1, monomial m2, monomial m3);


/**
 * @brief Determines the degree of a monomial.
 * @param m  The monomial.
 * @return The degree as an integer.
 */
int monomial_stddegree(monomial m);



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

/**
 * @brief Sets the ring_N.
 * @param n Number of Dimensions.
 * @return TRUE(=0) if function was successfull.
 */
int ring_set(int n);

/**
 * @brief Reads the number of variables out of a file and sets the ring dimension.
 * @param infile File/Stream consisting the number of dimensions.
 * @return TRUE(=0) if function was successfull.
 */
int ring_read(FILE *infile);

/**
* @brief Reads the variable of the inputfile.
* @param infile File/Stream consisting the variables.
* @return ASCII code of the variable.
*/
int ring_getvar(FILE *infile);



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
void binomial_puncture(binomial b,int position);
int ring_read(FILE *infile);
int binomial_equal(binomial b1, binomial b2);


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
