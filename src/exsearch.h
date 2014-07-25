/** 
*   @file exsearch.h 
*   @brief header file for the exhaustive search
*
*   @author Daniel Rembold
*   @bug No known bugs
*
*/

/**
* @brief An extra condition for degree compatible Gröbner basis for facet binomial to be flipped.
* @param b The given binomial.
* @return 1 for fullfilled condition, 0 for not fullfilled.
*/
int flip_ex_condition(binomial b,int degree_comp);

/**
 * @brief insert use next pointer in Gsets to implement a simple linked list.
 *               for now stored in no particular order.
 * @param g1 The gset to be inserted.
 * @param L List of gsets.
 */
void insert(gset g1, gset *L);

/**
 * @brief Determines if two gsets contain equivalent marked binomials.
 * @param g1 First binomial.
 * @param g2 Second binomial.
 * @return TRUE(=0) of FALSE(=1)
 */
int gset_equiv(gset g1, gset g2);

/**
 * @brief Checks a gset against those on list until an equivalent one is found.
 * @param g1 Gset which shall be found in list L.
 * @param L List that will be searched after g1.
 * @return Equivalent gset in L
 */
gset find(gset g1, gset L);

/**
 * @brief Breadth first search to find all degree compatible groebner bases
 *        of a code ideal.
 * @param g1 The starting Gröbner basis.
 * @return The number of groebner bases found.
 */
int exsearch(gset g1,int degree_comp, int no_print);