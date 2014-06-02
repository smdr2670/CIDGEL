/*
** match.h                                 Daniel Rembold
** -- header file with definitions 
**  
**
**
** Commented and slightly edited by Daniel Rembold
**
*/


struct node{
  int id;
  int nfacet;
  int nbinomials;
  int degree;
  struct node *next;
};


void sortedInsert(struct node** head, struct node* newnode);
struct node *newNode(int id, int facet, int binomials , int degree);
void printList(struct node *head);
int ListDelete(struct node **listP, int nfacet, int nbinomials , int degree);
int match( struct node **list1, struct node **list2);