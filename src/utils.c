/** 
*   @file utils.c
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
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


int eatwhite(FILE *is){
    int c;
    while ((c=fgetc(is))!=EOF){
        if (c=='%'){
        	while ((c=fgetc(is))!=EOF && c!='\n');
        }
        if (isspace(c)==0){
        	break;
        }
    }
    ungetc(c,is);
    return 0;
}

