#ifndef TIGERS_H
#define TIGERS_H


/**
 * Main algorithm using reverse search mentioned in algorithm 2.8
 * of this paper.
 *
 */



/**
 * @brief Prints vertex with all its information
 * @param g1 The given groebner base.
 * @param no_print Flag for printing.
 */
void vertex_print(gset g1,int no_print);

/**
 * @brief Extra condition for degree compatible groebner base for facet binomial to be flipped.
 * @param b The given binomial.
 * @return 1 for fullfilled condition, 0 for not fullfilled.
 */
int flip_condition_tig(binomial b);


int rsearch(gset g1,int number,int no_print);
#endif // TIGERS_H
