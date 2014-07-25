/** 
*   @file tigers.h
*   @brief Main algorithm using reverse search mentioned in algorithm 2.8 
*          of the TiGERS paper.
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


#ifndef TIGERS_H
#define TIGERS_H

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


/**
 *
 * @brief Main algorithm using the reverse search tree.
 * @param g1 Given Gröbner base to start of.
 * @param number Important in order to match 2 groebner fans.
 * @return Nummber of counted Gröbner bases.
 *
 */
int rsearch(gset g1,int number,int no_print);
#endif // TIGERS_H
