/** 
*   @file utime.h
*   @brief Function prototypes for time measuring.
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


/**
* @brief Sets the mark for time measuring, like a tic-command. 
* @return The time stamp.
*/
int set_mark(void);

/**
* @brief Gets the mark for time measuring, like a toc-command. 
* @return The absoulute time in seconds.
*/
int read_mark(int);
