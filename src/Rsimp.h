/** 
*   @file Rsimp.h
*   @brief header file for Rsimp.c with  definitions of linear programming  
*          data structure and basic implementation of revised simplex method.
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
#ifndef RSIMP_H
extern int LP_MAX_N;
extern int LP_MAX_M;
extern int LP_N;
extern int LP_M;
extern double **LP_A;
extern double *LP_B;
extern double *LP_C;
extern double *LP_X;
extern int *LP_Basis;
extern int *LP_NonBasis;
extern double **LP_Q;
extern double **LP_R;
extern double *LP_t1; 
extern double *LP_t2;
#endif

#define LP_A(i,j) LP_A[j][i]
#define LP_OPT 0
#define LP_UNBD 1
#define LP_FAIL -1



/**
* @brief Allocate space for LP data structures.
* @param M Number of rows.
* @param N Number of columns.
*/
void LP_get_space(int M, int N);

/**
* @brief Deallocates space for LP data structures.
*        Sets all LP globals to 0.
*
*/
void LP_free_space();

//void LP_print();


/**
* @brief print LP data structures to stdout.
*
*/
void Print_LP();

/**
* @brief revised simplex method (Using Bland's rule) 
*        and a qr factorization to solve the linear equations
* 
*     Adapted from algorithms presented in 
*              Linear Approximations and Extensions                  
*              (theory and algorithms)
*              Fang & Puthenpura
*              Prentice Hall, Engelwood Cliffs NJ (1993)
*       and 
*             Linear Programming
*             Chvatal 
*             Freeman and Company, New York, 1983
*  
*       (developed first in Octave, many thanks to the author)
*  
* 
*   Solve the problem 
*        minimize C'x, 
*        subject to A*x=b,  x>=0
*        for x,c,b n-vectors, and A an m,n matrix with full row rank 
*  
*  Assumptions:
*     A mxn matrix with full row rank.
*     b an m matrix. 
*     c an n-vector.
*     x an n-vector holding a basic feasible solution, 
*     basis m-vector holding indices of the basic variables in x
*     nonbasis n-m vector holding the indices not appearing in x.
*  
*  Returns: 
*       LP_FAIL if algorithm doesn't terminate.
*       LP_UNBD if problem is unbounded
*       LP_OPT  if optimum found
*   efects:
*     A,b,c unchanged.
*     x basis, nonbasis, hold info describing last basic feasible 
*                        solution.
*     Q,R hold qrdecomp of last basis matrix.
*     t1,t2 undefined.
*
* @param m Number of rows of A.
* @param n Nubmer of colums of A.
* @param A mxn matrix with full row rank.
* @param b m matrix.
* @param c an n-vector.
* @param x an n-vector holding a basic feasible solution, 
*          basis m-vector holding indices of the basic variables in x
*          nonbasis n-m vector holding the indices not appearing in x.
* @param R Matrix for the QR factorization.
* @param Q Matrix for the QR factorization.
* @param t1  Help-vector.
* @param t2  Help-vector.
* @return LP_FAIL if algorithm doesn't terminate.
*         P_UNBD if problem is unbounded
*         LP_OPT  if optimum found  
*/
int Rsimp(int m, int n, double **A, double *b, double *c,
          double *x, int *basis, int *nonbasis,
          double **R, double **Q, double *t1, double *t2);



/**
* @brief Use givens rotations on R to bring it into triangular form.
*        Store orthogonal matrix needed to bring R to triangular form in Q.
*        Assume R is an rxc matrix and Q is rxr.
* @param Q Matrix Q for the QR factorization.  
* @param R Matrix R for the QR factorization. 
* @param r Dimension of the R matrix.
* @param c Additional Dimension for the Q matrix.     
*/
void GQR(int r, int c, double **Q, double **R);






