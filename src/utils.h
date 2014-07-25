/** 
*   @file utils.h
*   @brief A place for general global definitions which must be accessible 
*          througout the program.
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
#define TRUE 0
#define FALSE 1


/**
* @brief Skips all whitespaces when reading a file.
* @param is Input stream.
* @return 0 for being successful
*/ 
int eatwhite(FILE *is);
