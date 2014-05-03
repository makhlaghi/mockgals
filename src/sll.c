/*********************************************************************
sll - simple linked list operations.

Copyright (C) 2013 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

sll is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sll is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sll. If not, see <http://www.gnu.org/licenses/>.

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
