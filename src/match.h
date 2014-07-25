/** 
*   @file match.h 
*   @brief Header files with function declarations for Gröbner fan matching.
*
*   @author Daniel Rembold
*   @bug No known bugs
*
*/


/**
* @brief Struct which contains all imporant informations of a Gröbner basis
*
*/
struct node{
  int id;			 /**< Number of the vertex of the edge graph #id.  */
  int nfacet;		 /**< The number of facet binomials, #nfacets.   */
  int nbinomials;    /**< Number of binomials #nelts.   */
  int degree;		 /**< Highest degree #deg.   */
  struct node *next; /**< Pointer to next node #next.   */
};

/**
 *@brief Inserts a node in the correct position in the linked list.
 *@param head Linked list.
 *@param newnode New node to be inserted.
 */
void sortedInsert(struct node** head, struct node* newnode);



/**
* @brief Creates a new Node (Vertex without binomials).
* @param id ID of the Vertex.
* @param facet Number of facets of the Vertex.
* @param binomials Number of binomials of the Vertex.
* @param degree Highest degree of the binomials.
* @return Newnode allocated space of a node.
*/
struct node *newNode(int id, int facet, int binomials , int degree);


//void printList(struct node *head);

  /**
  * @brief Searches and deletes a certain element of the list listP.
  * @param listP Linked List.
  * @param nfacet Number of facet.
  * @param nbinomials Number of binomials.
  * @param degree Highest degree.
  */
int ListDelete(struct node **listP, int nfacet, int nbinomials , int degree);


/**
 * @brief Checks if two Lists contains the same elements.
 * @param list1 First list.
 * @param list2 Second list.
 * @return 0 if lists are equal, -1 if not.
 *
 */
int match( struct node **list1, struct node **list2);