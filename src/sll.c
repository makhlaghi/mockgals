/*********************************************************************
MockGals - Make mock galaxies and stars from a catalog.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

MockGals is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MockGals is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MockGals. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sll.h"



/****************************************************************
 *****************            Float          ********************
 ****************************************************************/
void
printarrayoffsll(struct fsll **afsll, size_t num)
{
  size_t i;
  struct fsll *tmp;
  for(i=0;i<num;i++)
    {
      printf(" %lu:\n", i);
      for(tmp=afsll[i];tmp!=NULL;tmp=tmp->next)
	printf("%f, ", tmp->v);
      printf("\n");
    }
}





void
add_to_fsll(struct fsll **list, float value)
{
  struct fsll *newnode;

  newnode=malloc(sizeof *newnode);
  assert(newnode!=NULL);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void
pop_from_fsll(struct fsll **list, float *value)
{
  struct fsll *tmp;
  tmp=*list;
  *value=tmp->v;
  *list=tmp->next;
  free(tmp);
}





size_t
numinfsll(struct fsll *list)
{
  size_t num=0;
  struct fsll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    num++;
  return num;
}





void
fslltoarray(struct fsll *list, float **f, size_t *num)
{
  float *tf;
  size_t i=0;
  struct fsll *tmp;
 
  if(*num==0) *num=numinfsll(list);
  *f=malloc(*num*sizeof(float));
  assert(*f!=NULL);
  tf=*f;
    
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    tf[i++]=tmp->v;
}





void
freefsll(struct fsll *list)
{
  struct fsll *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}





void
freearrayoffsll(struct fsll **afsll, size_t num)
{
  size_t i;
  for(i=0;i<num;i++)
    freefsll(afsll[i]);
  free(afsll);
}




















/****************************************************************
 *****************           size_t          ********************
 ****************************************************************/
void
add_to_ssll(struct ssll **list, size_t value)
{
  struct ssll *newnode;

  newnode=malloc(sizeof *newnode);
  assert(newnode!=NULL);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void
pop_from_ssll(struct ssll **list, size_t *value)
{
  struct ssll *tmp;
  tmp=*list;
  *value=tmp->v;
  *list=tmp->next;
  free(tmp);
}





size_t
numinssll(struct ssll *list)
{
  size_t num=0;
  struct ssll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    num++;
  return num;
}





void
sslltoarray(struct ssll *list, size_t **f, size_t *num)
{
  size_t i=0, *tf;
  struct ssll *tmp;
 
  *num=numinssll(list);
  *f=malloc(*num*sizeof(size_t));
  assert(*f!=NULL);
  tf=*f;
    
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    tf[i++]=tmp->v;
}





void
freessll(struct ssll *list)
{
  struct ssll *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 ******************        Two way SLL    ***********************
 *****************           size_t          ********************
 ****************************************************************/

void
add_to_tssll_end(struct tssll **last, size_t value)
{
  struct tssll *newnode;

  assert(( newnode=malloc(sizeof *newnode) )!=NULL);

  newnode->v=value;
  newnode->next=*last;
  newnode->prev=NULL;
  if(*last)			/* If *list is not NULL */
    (*last)->prev=newnode;
  *last=newnode;
}





/* Note that start has to be initialized. */
void
pop_from_tssll_start(struct tssll **first,  size_t *value)
{
  struct tssll *tmp;
  tmp=*first;
  *value=tmp->v;
  *first=tmp->prev;
  free(tmp);
  if(*first)
    (*first)->next=NULL;
}




















/****************************************************************
 ******************        Ordered SLL       ********************
 *****************           size_t          ********************
 ****************************************************************/
/* We want to put the nodes in order based on the 'tosort' value of
each node. The top element should always have the smallest radius. */
void
add_to_ossll(struct ossll **list, size_t value, float tosort)
{
  struct ossll *newnode, *tmp=*list, *prev=NULL;

  assert(( newnode=malloc(sizeof *newnode) )!=NULL);

  newnode->v=value;
  newnode->s=tosort;

  /* *list points to the smallest value in the queue!*/
  while(tmp!=NULL)
    {
      if(tosort<tmp->s) break;
      /* No need for else, it will only come here if the condition
	 above is not satisfied. */
      prev=tmp;
      tmp=tmp->next;
    }

  if(tmp==NULL)	     /* This is the largest value so far. */
    {		     /* '*list' only changes if it is NULL. */
      newnode->next=NULL;
      if(prev) prev->next=newnode;   /* 'prev' is not NULL! */
      else     *list=newnode;	     /* Only for initial node. */
    }
  else
    {
      if(prev) prev->next=newnode;
      else     *list=newnode;	/* 'tosort' is smaller than all. */
      newnode->next=tmp;
    }
}





/* Note that the popped element is the smallest! */
void
pop_from_ossll(struct ossll **list,  size_t *value, float *sortvalue)
{
  struct ossll *tmp;
  tmp=*list;
  *value=tmp->v;
  *sortvalue=tmp->s;
  *list=tmp->next;
  free(tmp);
}





/* Add the elements of an ossll to a ssll. */
void
ossll_into_ssll(struct ossll *in, struct ssll **out)
{
  struct ossll *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      add_to_ssll(out, in->v);
      free(in);
      in=tmp;
    }
}




















/****************************************************************
 ******************   Two way, Ordered SLL   ********************
 *****************           size_t          ********************
 ****************************************************************/
void
print_tossll(struct tossll *l, struct tossll *s)
{
  size_t counter=1;   /* We are not counting array elements :-D ! */
  while(l!=NULL)
    {
      printf("\t%-5lu (%lu, %.4f) \n", counter++, 
	     l->v, l->s);
      l=l->next;
      printf("\t\t\t\t(%lu, %.4f)\n", s->v, s->s);
      s=s->prev;
    }
  printf("\n");
}





/* Very similar to Ordered SLL, but now it is two way. */
void
add_to_tossll_end(struct tossll **largest, struct tossll **smallest, 
		  size_t value, float tosort)
{
  struct tossll *newnode, *tmp=*largest;

  assert(( newnode=malloc(sizeof *newnode) )!=NULL);

  newnode->v=value;
  newnode->s=tosort;
  newnode->prev=NULL;

  while(tmp!=NULL)
    {
      if(tosort >= tmp->s) break;
      /* No need for else, it will only come here if the condition
	 above is not satisfied. */
      newnode->prev=tmp;
      tmp=tmp->next;
    }

  if(tmp==NULL)	     /* This is the smallest value so far.     */
    {		     /* '*largest' only changes if it is NULL. */
      newnode->next=NULL;
      *smallest=newnode;
      if(newnode->prev)		/* 'prev' is not NULL! */
	newnode->prev->next=newnode;   
      else			/* 'prev is NULL, Only first. */
	*largest=newnode;
    }
  else
    {
      if(newnode->prev)
	{ 
	  newnode->prev->next->prev=newnode;
	  newnode->prev->next=newnode;
	}
      else
	{
	  (*largest)->prev=newnode;
	  *largest=newnode;       /* 'tosort' is larger than all. */
	}
      newnode->next=tmp;
    }
}





/* Note that start has to be initialized. */
void
pop_from_tossll_start(struct tossll **largest, struct tossll **smallest,  
		      size_t *value, float *tosort)
{
  struct tossll *tmp=*smallest;

  *value=tmp->v;
  *tosort=tmp->s;

  *smallest=tmp->prev;
  free(tmp);
  if(*smallest)
    (*smallest)->next=NULL;
  else
    *largest=NULL;

  /*printf("Popped v: %lu, s: %f\n", *value, *tosort);*/
}





void
smallest_tossll(struct tossll *largest, struct tossll **smallest)
{
  struct tossll *tmp=largest;

  while(tmp!=NULL)
    {
      if(tmp->next==NULL)
	{
	  *smallest=tmp;
	  break;
	}
      tmp=tmp->next;
    }

  /* If *largest wasn't NULL initially, tmp should not be NULL because
     the loop terminated before it becomes null. But if it was
     initiall NULL, it will never enter the loop, and so smallest
     should also be NULL. */
  if(tmp==NULL) *smallest=NULL;
}




void
tossll_into_ssll(struct tossll *in, struct ssll **out)
{
  struct tossll *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      add_to_ssll(out, in->v);
      free(in);
      in=tmp;
    }
}
