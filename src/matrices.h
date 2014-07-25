/** 
*   @file matrices.h
*   @brief Definition and Minipulation of integer and double matrices and vectors.
*   	   vectors are simply C-vectors (with indices starting at zero) and matrices
*   	   are stored Numerical Recepies style -- as a vector of pointers to the
*   	   rows of the matrix [which are themselves just regular C vectors].
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
#define IMref(M,i,j) ((M)[i][j])


/**
* @brief read in description of imatrix from is, create matrix and fill it.
*      -- format: { m n : entry_1 .... entry_mn } where m=number of rows n=number of columns
*     e.g. 
*     { 2 4 10 : 1 0 3 5 0 1 8 9 } describes the matrix [ 1, 0, 3, 5] 
*                                                       [ 0, 1, 8, 9]   in F_10
* @param is  Input stream.    
* @param m Code dimension.
* @param n Code length.
* @param f Dimension of the ring.
*/
int **imatrix_read(FILE *is,int *m, int *n, int *f);


/**
* @brief copy prefex string to output file of, then write ascii representation
*        of matrix with m-rows, n-cols and entrees in M to of.
* @param Input stream.
* @param Character prefix.
* @param M Matrix M.
* @param m Number of rows.
* @param n Number of columns.
*/
void print_imatrix(FILE *of, char *prefix,int **M, int m, int n);


/**
* @brief Allocates space for a integer matrix.
* @param r Number of rows.
* @param c Number of colums
* @return Pointer to the matrix.
*/
int **new_imatrix(int r, int c);


/**
* @brief Deallocates the space of a matrix.
* @param M Matrix to be deleted.
*/
void free_imatrix(int **M);

/**
* @brief Allocates space for a new integer vector.
* @param n Length of the new vector.
*/
int *new_ivector(int c);

/**
* @brief Deletes a vector.
* @param M The chosen vector.
*/
void free_ivector(int *M);


/**
* @brief Allocates space for a double matrix.
* @param r Number of rows.
* @param c Number of colums
* @return Pointer to the matrix.
*/
double **new_matrix(int r, int c);

/**
* @brief Deallocates the space of a matrix.
* @param M Matrix to be deleted.
*/
void free_matrix(double **M);

/**
* @brief Allocates space for a new double vector.
* @param n Length of the new vector.
*/
double *new_vector(int c);

/**
* @brief Deletes a vector.
* @param M The chosen vector.
*/
void free_vector(double *V);




