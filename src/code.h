/** 
*   @file code.h
*   @brief Function prototype for computing code ideals.
*
*   @author Daniel Rembold
*   @bug No known bugs
*
*/


#ifndef CODE_H
#define CODE_H

/**
 * @brief Given an mxn integer matrix M, compute an rgb for the code ideal I_C.
 *        Expects generator matrix in standard form, reads off equations for a generating system, then
 *        uses Buchberger algorithm.
 * @param M input Matrix
 * @param m number of rows
 * @param n number of columns
 */
gset gset_code_ideal(int **M,int m, int n);

#endif