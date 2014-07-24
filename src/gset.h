/** 
*   @file gset.h
*   @brief header file with definitions and basic operations on generating sets
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
#include "binomial.h"
#define DEBUG 0

typedef struct gset_tag *gset;

/* Linked List of gset_tag which contiains the binomial and the
   caching informations*/
struct gset_tag{
    int id;
    int nfacets;
    int nelts;
    int deg;
    binomial bottom;
    binomial cache_edge;
    struct gset_tag *cache_vtx;
    struct gset_tag *next;
};

#define gset_first(g) (g->bottom)                            
#define gset_cache_vtx(g) (g->cache_vtx)
#define gset_cache_edge(g) (g->cache_edge)
#define gset_nelts(g) (g->nelts)
#define gset_nfacets(g) (g->nfacets)
#define gset_id(g) (g->id)
#define gset_deg(g) (g->deg)
#define gset_next(g) (g->next)


/**
 * @brief Creates an empty gset, initializes all its pointers to null.
 * @return The new gset.
 */
gset gset_new();

/**
 * @brief Deletes a gset and all binomials in it.
 * @param g The gset to be deleted.
 */
void gset_free(gset g);


/**
 * @brief Reads a gset from an input stream.
 * @param is The input stream.
 * @param g  All informaion stored in gset.
 * @return TRUE(=0) if gset is read in successfully.
 */
int gset_read(FILE *is,gset g);

/**
 * @brief Prints a gset.
 * @param of The output stream.
 * @param g The gset to be printed.
 */
void gset_print(FILE *of,gset g);

/**
 * @brief Displays the initial/leading ideal of Gset.
            
            See definition 2.5 of the thesis for defintion
            of a leading/initial ideal. 
 * @param of The output stream.
 * @param g The initial/leading ideal of the gset g to be printed.
 */
void gset_init_print(FILE *of,gset g);

/**
 * @brief Places binomials in gset in decreasing lexicographic order.
 *        (using binomial_compair(,) to determine relative order or binomials)
 * @param g The gset which gets the new binomial.
 * @param b The binomial to be inserted.
 */
void gset_insert(gset g, binomial b);



/**
 * @brief Returns the first facet binomial which is mis-marked with respect to the 
 *		  lexicographic order.
 *               A mis-marked binomial is a binomial that does not suit
 *               with the term order, but it is a facet binomial.
 * @param g The gset g to be  examined.
 * @return The first facet binomial which is mis-marked.
 */
binomial gset_downedge_lex(gset g);

/**
 * @brief Returns the first facet binomial which is mis-marked with respect to the 
 *		  graded-lexicographic order.
 * @param g The gset g to be  examined.
 * @return The first facet binomial which is mis-marked.
 */
binomial gset_downedge_grlex(gset g);


/*
 * @brief bmleadreduce "lift" monomial m1 to a binomial using bmlist B,mlist.
            This function represents line 6 and 7 from algorithm 5 of this thesis.
 * @param m1 The binomial to be lifted.
 * @param B The binomial which may divide m1.
 * @param mlist The list of monomials to compare with.
 * @return The new generated binomial.
 *
 */
binomial bmleadreduce(binomial m1, binomial B, binomial mlist);


/**
* @brief Checks if a binomial is divisible by a list on binomials.
* @param S The binomial which is the divident.
* @param L List of binimials which may divide S.
* @return The binomial which divides S.
*/
binomial find_divisor(binomial S,binomial L);

/**
* @brief Passes through a list L, removing all multiples of S encountered (except S if it is on the list).
* @param S Binomial with multiples to be removed.
* @param L List of binomials.
*
*/
void remove_multiples(binomial S, binomial *L);

/**
* @brief Transforms a marked binomial and a list of monomials into the
*        reduced grobner basis for the ideal to generate.
* @param B  The marked binomial
* @param ms List of binomials.
*/
void bmrgb(binomial B, binomial *ms)


/**
* @brief Take a gset and one of its facet binomials and produce the new gset
*        corresponding to the grobner cone which shares the given facet binomial.
*        This function represents algorithm 5 in this thesis.
* @param g1 A reduced Gröbner basis.
* @param b A prescribed facet binomial.
*
*/
gset gset_flip(gset g1, binomial b);


/**
* @brief Autoreduces a Gröbner basis with respect to its marking.
* @param g Given Gröbner basis.
*
*/
void gset_autoreduce();



/**
* @brief Checks whether b is known to be a facet or not. if not call lp_isfacet()
*        to use a linear program to set b's facet flag.
* @param g Gröbner basis to be researched.
* @param b Binomial to be checked if it is a facet.
* @return  Interger-flag for facets.
*/
int gset_isfacet(gset g,binomial b);

/**
* @brief Uses repeated calls to gset_isfacet to ensure that all facet flags are set.
*        Also sets the flags nfacets, nelts and deg counters. 
*        These contain garbage except right after a call to set facets.
* @param g Gröbner basis to be set.
*/
void gset_setfacets(gset g);

/**
* @brief Checks if a binomial is flippable without linear programming.
* @param g Gröbner basis to be researched.
* @param b Binomial to be checked if it is a facet.
* @return  Interger-flag for facets.
*/
int gset_isflippable(gset g,binomial b);

gset gset_toric_ideal(int **,int,int);

/**
* @brief Use buchberger's algorithm to transform g into groebner basis
*        wrt the term order implemented by comp.
* @param g The given ideal.
* @param comp Function pointer with standard monomial orders.
*/
void gset_rgb(gset g, int (*comp)(monomial,monomial));

/**
 * @brief Displays facet binomials.
 * @param of The output stream.
 * @param g The given gset.
 */
void gset_facet_print(FILE *of,gset g);


/**
* @brief Processes S-Pair defined by binomials b1 and b2.
*        Taking S-Pair with respect to union of DL SL and *NL lists.
*        If result is non-zero add result to NL list.  
* @param b1 Binomial.
* @param b2 Binomial.
* @param DL Binomial.
* @param SL Binomial.
* @param NL List of binomials.  
* @param comp Function pointer with standard monomial orders.
*/
void process_pair(binomial b1,binomial b2,binomial DL,binomial SL,binomial *NL,
                  int (*comp)(monomial,monomial));


/**
 * @brief Adds the nonprime Terms to a generating set in order to get a CodeIdeal,
 *        check all first terms of all binomials.
 *        If the exponent vector has the degree of 1, then it does not appear in
 *        the non-prime set.
 *        example: if no leading term divides x1^2, add x1^2-1 to the set
 * @param g gset which gets the nonprime binomials
 */
void gset_add_nonprime(gset g);


/**
 * @brief gset_only_degreecompatible Determines whether a Groebner Basis is the only degree compatible Groebner Basis or not.
 *        Example: In the 3 dimensional vector space, the facet of the Groebner Fan is the only
 *                 one which contains the all-one vector.     
 * @param g Gröbner basis to be researched.
 * @return 1 if it is the only degree compatible Gröbner basis and 0 if not.
 */
int gset_only_degreecompatible(gset g);

/**
* @brief Counts the degree compatible binomials in a Gröbner basis
* @param g1 Gröbner basis to be researched.
* @return Number of d.c. binomials
*/
int countDC(gset g1);

/**
* @brief Sets up the linear program to check if a binomial is a facet.
* @param g Gröbner basis consisting the binomial.
* @param b Binomial to be researched.
* @return Interger-flag for facets.
*
*/
int lp_isfacet(gset g,binomial b);


/**
* @brief Copies a Gröbner basis.
* @param src  Source Gröbner basis
* @param dest Destination Gröbner basis
*/
void gset_copy(gset G1, gset G2);


/**
* @brief Punctures a linear code.
*        Simply sets the all exponent-vectors from the gset to 0 at the certain position.
* @param G1 Gröbner basis to be punctured.
* @param pos Variable that will be punctured.
*
*/
void gset_puncture(gset G1, int pos);

/**
 * @brief Deletes a binimoal from a gset.
 * @param g The gset to be researched.
 * @param b The binomial which shall be resarched.
 */
void gset_delete_binomial(gset g, binomial b);

/**
 * @brief Eliminates '1-1' terms in a gset.
 * @param g The gset to be researched.
 */
void gset_eliminatezero(gset g);







