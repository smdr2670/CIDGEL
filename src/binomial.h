/** 
*   @file binomial.h   
*   @brief header file for binomal type and main operations on binomials
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


#ifndef BINOMIAL_H 
extern int ring_N;
extern int ring_lv;
extern int *ring_weight;
extern int code_dim;
#endif

void code_dim_set(int n);
int ring_set(int n);
int ring_getvar(FILE *ifile);

#define ring_putvar(ofile,v) fprintf(ofile,"%c",'a'+v)


typedef int *monomial;

/**
 * @brief Determines if m1 divides m2.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_divides(monomial m1, monomial m2);


/**
 * @brief Determines if monomials m1 and m2 have same exponents.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_equal(monomial m1, monomial m2);


/**
 * @brief Determines if m1 is lexicographically greater then m2.
 * @param m1 The first monomial.
 * @param m2 The second monomial.
 * @return >0 if m1 is lexicographically greater then m2.
 */
int monomial_lexcomp(monomial m1, monomial m2);


/**
 * @brief Determines if m1 is degree reverse lexicographically greater than m2
 *        with variables taken in order a>b>..., except that variable ring_lv is taken last (if it is >=0).
 * @param m1 The first monomial.
 * @param m2 The second binomial.
 * @return >0 if m1 is degree reverse lexicographically greater than m2.
 */
int monomial_rlexcomp(monomial m1, monomial m2);


/**
 * @brief Determines if m1 is degree lexicographically greater than m2.
 * @param m1 The first monomial.
 * @param m2 The second monomial.
 * @return >0 if m1 is degree lexicographically greater than m2.
 */
int monomial_grlexcomp(monomial m1, monomial m2);


/**
 * @brief Determines if m1 is reletively prime to m2.
 * @param m1 First monomial.
 * @param m2 Second monomial.
 * @return TRUE(=0) or FALSE(=1).
 */
int monomial_rel_prime(monomial m1, monomial m2);


/**
 * @brief Prints monomial in a file.
 * @param of Outputfile where monomial will be printed.
 * @param exps Vector which contains the exponents.
 */
void print_monomial(FILE *of, int *exps);


/**
 * @brief Reads ascii representation of monomial as power product
 *         onvert to exponent vector stored in *exps.
 * @param is The input stream.
 * @param exps Store the exponent vector in exps.s
 *
 */
void get_monomial(FILE *is,int *exps);


/**
 * @brief monomial_lcm copy least common multiple of m1 and m2 into m3.
 * @param m1 first monomial
 * @param m2 second monomial
 * @param m3 lcm monomial
 * @return
 */
void monomial_lcm(monomial m1, monomial m2, monomial m3);


/**
 * @brief Determines the degree of a monomial.
 * @param m  The monomial.
 * @return The degree as an integer.
 */
int monomial_stddegree(monomial m);



typedef struct bin_tag *binomial;

/**
* @brief A structure to represent binomials
*/
struct bin_tag{
	int *exps1; 	/**< Stores the exponent vector of the first monomial #exps1.	*/
    int *exps2; 	/**< Stores the exponent vector of the second monomial #exps2.	*/
    int *E; 		/**< Points to exps1 and exps2 and is used for allocating and deallocating #E. */
    int ff; 		/**< Flag which shows if a binomial is a facet binomial or not #ff.	*/
    int bf; 		/**< Tells if there is a monomial or binomial #bf.	*/
    binomial next;  /**< Pointer to next binomial #next.  */
};


#define BINOMIAL 1
#define MONOMIAL 0
#define FACET    0
#define NONFACET 1
#define UNKNOWN -1

#define binomial_next(b) (b->next)
#define binomial_lead(b) (b->exps1)
#define binomial_trail(b) (b->exps2)
#define monomial_set(b)   (b->bf=MONOMIAL)
#define binomial_set(b)   (b->bf=BINOMIAL)
#define binomial_facet(b) (b->ff)


/**
 * @brief Allocate space for a new binomial.
 *        Initializes flags and set exponents to zero.
 * @return New allocated binomial.
 */
binomial binomial_new();

/**
 * @brief Sets the ring_N.
 * @param n Number of Dimensions.
 * @return TRUE(=0) if function was successfull.
 */
int ring_set(int n);

/**
 * @brief Reads the number of variables out of a file and sets the ring dimension.
 * @param infile File/Stream consisting the number of dimensions.
 * @return TRUE(=0) if function was successfull.
 */
int ring_read(FILE *infile);

/**
* @brief Reads the variable of the inputfile.
* @param infile File/Stream consisting the variables.
* @return ASCII code of the variable.
*/
int ring_getvar(FILE *infile);


/**
 * @brief Reclaims space allocated by binomial_new().
 * @param m Binomial which shall be deleted.
 */
void binomial_free(binomial m);


/**
 * @brief Reads file and stores a binomial.
 * @param is The input file stream.
 * @param b Read binomial will be stored in this parameter.
 */
void binomial_read(FILE *is, binomial b);



/**
 * @brief Prints the binomial.
 * @param of The output stream.
 * @param b Binomial which shall be printed.
 */
void binomial_print(FILE *of, binomial b);


/**
 * @brief Forms S-Pair between marked binomial and marked monomial
 *        and reduce result by binomial as many times as possible.
 * @param b First binomial.
 * @param m The Mononomial.
 * @return S-Pair between the parameters.
 */
binomial monomial_spair(binomial b,monomial m);


/**
 * @brief Interchange leading and trailing terms of a binomial.
 * @param b Binomial which will be flipped.
 */
void binomial_flip(binomial b);


/**
 * @brief Replaces b1, by spair of b1 and b2 -- the spair is taken with respect
 *        to the markings of b1 and b2, but it's marking is not meaningfull.
 * @param b1 The first binomial.
 * @param b2 The second binomial.
 * @return TRUE(=0) if result is equivalent to zero monomial.
 */
int binomial_spair(binomial b1, binomial b2);


/**
 * @brief Calculates the degree of the first term.
 * @param b The given binomial.
 * @return Degree of first term.
 */
int binomial_first_term_degree(binomial b);


/**
 * @brief Returns whether a binomial is degree compatible or not.
 * @param b The given binomial.
 * @return 1 if first term has a higher degree, 0 if degrees are equal, -1 if second term has a higher degree.
 */
int binomial_degree_compatible(binomial b);


/**
 * @brief Copies source binomials exponents and flags to those of the destination binomials.
 * @param src Source binomial.
 * @param dest Destination binomial.
 */
void binomial_copy(binomial src,binomial dest);


/**
 * @brief Returns the position which will NOT belong to the non-prime ideal.
 *        Note that the function is only called when the leading term of the binomial
 *        has the degree 1.
 * @param b The given binomial.
 * @return Index i of exponent vector with i != 0.
 */
int binomial_variable_position(binomial b);


/**
 * @brief Assuming that lead(b2) divides trail(b1) change b1 into binomial
 *        resulting from bumping b2 so that lead(b2)=trail(b1) and adding
 *        result to b1.
 * @param b1 The first binomial.
 * @param b2 The second binomial.
 */
void reducetrail(binomial b1, binomial b2);


/**
 * @brief Assuming that the leading term of b1 divides that of b2
 *        multiply both sides of b1 so that the leading terms of b1
 *        and b2 are equal.
 * @param b1 The first binomial.
 * @param b2 The second binomial.
 */
void binomial_bumpto(binomial b1, binomial b2);


/**
 * @brief Determines if the lead term of b1 is lexicographically greater than
 *        that of b2 or if they tie return true if the trail term b1 is greater
 *        or equal to that of b2 otherwise return false.
 *
 * @param b1 The first binomial.
 * @param b2 The second binomial.
 * @return TRUE(=0) or FALSE (=1).
 */
int binomial_compair(binomial b1,binomial b2);

/**
* @brief Punctures the binomial with cancelling a certain variable.
* @param b The binomial to be punctured.
* @param position Position of the exp. vector which will be set to zero.
*/
void binomial_puncture(binomial b,int position);



/**
* @brief Checks if two binomials are identical.
* @param b1 First binomial.
* @param b2 Second binomial.
*/
int binomial_equal(binomial b1, binomial b2);


/*
** return TRUE if binomials marking agrees with respect to the
** lexicographic term order
*/
#define binomial_lexordered(b)\
    ((monomial_lexcomp(binomial_lead(b),binomial_trail(b))>=0) ? TRUE : FALSE)

/*
** return TRUE if binomials marking agrees with respect to the
** graded lexicograpic term order
*/
#define binomial_grlexordered(b)\
    ((monomial_grlexcomp(binomial_lead(b),binomial_trail(b))>=0) ? TRUE : FALSE)
