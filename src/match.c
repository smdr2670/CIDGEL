/** 
*   @file match.c 
*   @brief Function definitions for Gröbner fan matching.
*
*   @author Daniel Rembold
*   @bug No known bugs
*
*/

#include <stdio.h>
#include <stdlib.h>
#include "match.h"


void sortedInsert(struct node** head, struct node* newnode){
    struct node* tmp;
    /* Special case for the head end */
    if (*head == NULL || (*head)->nfacet >= newnode->nfacet)
    {
        newnode->next = *head;
        *head = newnode;
    }
    else
    {
        /* Locate the node before the point of insertion */
        tmp = *head;
        while (tmp->next!=NULL && tmp->next->nfacet < newnode->nfacet){
            tmp = tmp->next;
        }
        newnode->next = tmp->next;
        tmp->next = newnode;
    }
}





struct node *newNode(int id, int facet, int binomials , int degree){
                            /* allocate node */
    struct node* newnode = (struct node*) malloc(sizeof(struct node));

                            /* put in the data  */
    newnode->id  = id;
    newnode->nfacet = facet;
    newnode->nbinomials = binomials;
    newnode->degree = degree;
    newnode->next =  NULL;

    return newnode;
  }

/* Function to print linked list */
void printList(struct node *head){
    struct node *temp = head;
    while(temp != NULL)
    {
      printf("(%d,%d,%d)\n",  temp->nfacet, temp->nbinomials, temp->degree);
      temp = temp->next;
    }
  }


  int ListDelete(struct node **listP, int nfacet, int nbinomials , int degree){
      struct node *currP, *prevP;

      /* For 1st node, indicate there is no previous. */
      prevP = NULL;

      /*
       * Visit each node, maintaining a pointer to
       * the previous node we just visited.
       */
      for (currP = *listP; currP != NULL; prevP = currP, currP = currP->next) {

        if (currP->nfacet == nfacet && currP->nbinomials == nbinomials && currP->degree == degree){ /* Found it. */
              if (prevP == NULL) {
                /* Fix beginning pointer. */
                *listP = currP->next;
              } else {
                /*
                 * Fix previous node's next to
                 * skip over the removed node.
                 */
                prevP->next = currP->next;
              }

              /* Deallocate the node. */
              free(currP);

              /* Done searching. */
              return 0;
        }
      }

      return -1;

}


int match( struct node **list1, struct node **list2){

  if(*list1 == NULL || *list2== NULL){
    printf("EMPTY LIST!");
    return(-1);
  }

  struct node *tmp;

  for( tmp = *list1; tmp!= NULL; tmp=tmp->next){
    if(ListDelete(list2,tmp->nfacet,tmp->nbinomials,tmp->degree) == 0){
      // Found and deleted in second list, delete in first list too and move on
      ListDelete(list1,tmp->nfacet,tmp->nbinomials,tmp->degree);
    }else{
      //Element not found => codes are not compatible => returning
      return (-1);
    }
  }

  return 0;
}



