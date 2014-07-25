/** 
*   @file utime.c
*   @brief Function definitions for time measuring.
*          
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
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "utime.h"



int set_mark(){
    return clock();
}


int read_mark(int timeset){
    return (clock()-timeset) / CLOCKS_PER_SEC;
}





