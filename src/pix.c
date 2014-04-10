/*********************************************************************
pix - Define, find and manipulate a 2D pixel structure in images.

Each pixel will have a structuring element, keeping its index 
(in the 1D array), its value and its neighbors. This file
contains functions to define, create and manipulate such structures
which will become very handy in all sorts of image processing jobs.
Like sorting, segmentation and so forth.

Copyright (C) 2013 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

pix is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pix is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pix. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "pix.h"
#include "fitsarrayvv.h"




















/****************************************************************
 *****************     pixlist functions     ********************
 ****************************************************************/
void
addtopixlist(struct pixlist **queue, struct pix *p)
{
  struct pixlist *newnode;

  newnode=malloc(sizeof *newnode);
  assert(newnode!=NULL);

  newnode->p=p;
  newnode->ignore=0;
  newnode->next=*queue;
  *queue=newnode;
}





/* This is very similar to branching in Git!  We have a simple linked
   list: we want two entry points into this list for two sets, one of
   the sets completely over laps the other.  For example 8
   connectivity and 4 connectivity.  All 4-connected pixels are also
   8-connected.  So first we make a list with 4-connectivity, then we
   take the starting point of that and add it to the 8 connected
   pixels.  */
void
addtoanotherpixlist(struct pixlist **newq, 
        struct pixlist *oldq, struct pix *p)
{
  struct pixlist *newnode;

  newnode=malloc(sizeof *newnode);
  assert(newnode!=NULL);

  newnode->p=p;
  newnode->ignore=0;
  newnode->next=oldq;
  *newq=newnode;
}





void
popfrompixlist(struct pixlist **queue, struct pix **p)
{
  struct pixlist *tmp;
  tmp=*queue;
  *p=tmp->p;
  *queue=tmp->next;
  free(tmp);
}





void
freepixlist(struct pixlist *queue)
{
  struct pixlist *tmp=queue, *ttmp;;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}





void
printpixlist(struct pixlist *queue)
{
  struct pixlist *tmp;
  for(tmp=queue;tmp!=NULL;tmp=tmp->next)
    printf("%lu, ", tmp->p->i);
  printf("\n");
}




















/****************************************************************
 *****************    pix array functions    ********************
 ****************************************************************/
void 
printpixarray(struct pix *D, size_t size)
{
  size_t i, j;
  for(i=0;i<size;i++)
    {
      printf("%lu: %g\n", D[i].i, D[i].v);
      for(j=0;D[i].ngb[j].p!=NULL;j++)
	printf("%lu, ", D[i].ngb[j].p->i);
      printf("\n\n");
    }
}





/* In the pix structure, there are three pointers to neighbors. Two
   to the list of 4 (ngb4) and 8 (ngb8) connected neighbors. One
   called ngb which can be set to any one of the two above.
   The job of this function is to set this third pointer to 
   either of the 4 or 8 connected neighbors. */
void
setngb(struct pix *D, size_t size, int con_type)
{
  size_t i;
 
  assert(con_type==4 || con_type==8);
    
  if(con_type==4)
    for(i=0;i<size;i++)
      D[i].ngb=D[i].ngb4;
  else
    for(i=0;i<size;i++)
      D[i].ngb=D[i].ngb8;
}





void
addpixval(struct pix *D, float *img, size_t size)
{
  size_t i;
  for(i=0;i<size;i++)
    D[i].v=img[i];
}




void
resetignores(struct pix *D, size_t size)
{
  size_t i;
  struct ngbarrpix *ngbarr;

  ngbarr=D[0].ngb8;
  for(i=0;i<size*NGBSCOLS;i++)
    ngbarr[i].ignore=0;
}




/* Note that the first element of the neighbors array is pointed to
     by the ngb8 pointer of D[0]. */
void 
freepixarray(struct pix *D)
{
  free(D[0].ngb8);
  free(D);
}




















/****************************************************************
 *****************      Sortig pix array     ********************
 ****************************************************************/
int
pix_v_increasing(const void *a, const void *b)
{
  float ta=(*(struct pix * const *)a)->v; 
  float tb=(*(struct pix * const *)b)->v; 
  return (ta > tb) - (ta < tb);
}





int
pix_v_decreasing(const void *a, const void *b)
{
  float ta=(*(struct pix * const *)a)->v; 
  float tb=(*(struct pix * const *)b)->v; 
  return (tb > ta) - (tb < ta);
}





int
pix_i_increasing(const void *a, const void *b)
{
  size_t ta=(*(struct pix * const *)a)->i; 
  size_t tb=(*(struct pix * const *)b)->i; 
  return (ta > tb) - (ta < tb);
}




















/****************************************************************
 *****************     Read image to pix     ********************
 ****************************************************************/
/* Used in imgtopix(), similar for all corner pixels. */
#define CORNERFILL   d[ind].ngb8=&ngbs[ind*NGBSCOLS];       \
  d[ind].ngb4=&ngbs[ind*NGBSCOLS+1];                        \
  for(i=3;i<NGBSCOLS-1;i++) ngbs[ind*NGBSCOLS+i].p=NULL;    \

/* Used in imgtopix(), similar for all side pixels. */
#define SIDEFILL   d[ind].ngb8=&ngbs[ind*NGBSCOLS];	    \
  d[ind].ngb4=&ngbs[ind*NGBSCOLS+2];                        \
  for(k=5;k<NGBSCOLS-1;k++) ngbs[ind*NGBSCOLS+k].p=NULL;    \


/* Convert an array into a graph. 

   Input: img: the array.  s0: its zeroth axis size.  s1: its first
          axis size.  Pointer to an array of "struct pix" that will be
          allocated here.

   The pixel graphs are organized like this. Assuming the image has
   s0*s1=N pixels then after this function finished we have:

   struct pix d[N]: Array of N "struct pix" elements.  

   sturct ngbpixarr ngbs[N*9]: Array of N*9 "struct ngbpixarr"
   elements. This array will show the neighboring relations.

   By default the neighbors will be set to 8 connected.
   
   The reason the neighbors array has 9 columns, while at most we
   have 8 columns is that we want to search blindly (without worring
   about how many neighbors a pixel has for all pixels) in them and
   the only way to do that is to have a NULL element as the last for
   those pixels that do actually have 8 neighbors. 

   Note on the ignore section of struct ngbarrpix: This will be used
   in objects, when a pixel belongs to an object, but one or some of
   its neighbors don't. The NULL value of the pointer in the array
   will only be used for the actual image (which is most of our
   applicaation. So initially, all ignore values are set to zero. 

   pix.ngb8 points to the first element in the corresponding
   row. After the only 8 connected neighbors are finished, the 4
   connected neighbors will begin. pix.ngb4 points to the first of
   those. */
void
imgtopix(size_t s0, size_t s1, struct pix **D)
{
  struct pix *d; 
  size_t ind, i, j, k; 		/* k used in SIDEFILL */
  struct ngbarrpix *ngbs;

  d=malloc(s0*s1*sizeof *d);
  assert(d!=NULL);
  ngbs=malloc(NGBSCOLS*s0*s1*sizeof *ngbs);
  assert(ngbs!=NULL);
  
  for(i=0;i<s0*s1;i++)
    d[i].i=i;    

  for(i=0;i<NGBSCOLS*s0*s1;i++)
    ngbs[i].ignore=0;
  for(i=0;i<s0*s1;i++)
    ngbs[i*NGBSCOLS+NGBSCOLS-1].p=NULL;
  
  ind=0;			/* Bottom Right */
  ngbs[ind*NGBSCOLS].p  =&d[ind+s1+1];
  ngbs[ind*NGBSCOLS+1].p=&d[ind+1];
  ngbs[ind*NGBSCOLS+2].p=&d[ind+s1];
  CORNERFILL

  ind=s1-1;			/* Bottom Left */
  ngbs[ind*NGBSCOLS].p  =&d[ind+s1-1];
  ngbs[ind*NGBSCOLS+1].p=&d[ind-1];
  ngbs[ind*NGBSCOLS+2].p=&d[ind+s1];
  CORNERFILL

  ind=(s0-1)*s1;		/* Top Right */
  ngbs[ind*NGBSCOLS].p  =&d[ind-s1+1];
  ngbs[ind*NGBSCOLS+1].p=&d[ind+1];
  ngbs[ind*NGBSCOLS+2].p=&d[ind-s1];
  CORNERFILL
  
  ind=s0*s1-1;			/* Top Left */
  ngbs[ind*NGBSCOLS].p  =&d[ind-s1-1];
  ngbs[ind*NGBSCOLS+1].p=&d[ind-1];
  ngbs[ind*NGBSCOLS+2].p=&d[ind-s1];
  CORNERFILL
  
  for(j=1;j<s1-1;j++)
    {
      ind=j;			/* Bottom */
      ngbs[ind*NGBSCOLS].p  =&d[ind+s1+1];
      ngbs[ind*NGBSCOLS+1].p=&d[ind+s1-1];
      ngbs[ind*NGBSCOLS+2].p=&d[ind+1];
      ngbs[ind*NGBSCOLS+3].p=&d[ind-1];
      ngbs[ind*NGBSCOLS+4].p=&d[ind+s1];
      SIDEFILL
    }

  for(j=1;j<s1-1;j++)
    {
      ind=(s0-1)*s1+j;		/* Top */
      ngbs[ind*NGBSCOLS].p  =&d[ind-s1+1];
      ngbs[ind*NGBSCOLS+1].p=&d[ind-s1-1];
      ngbs[ind*NGBSCOLS+2].p=&d[ind+1];
      ngbs[ind*NGBSCOLS+3].p=&d[ind-1];
      ngbs[ind*NGBSCOLS+4].p=&d[ind-s1];
      SIDEFILL
    }

  for(i=1;i<s0-1;i++)
    {
      ind=i*s1;			/* Left */
      ngbs[ind*NGBSCOLS].p  =&d[ind-s1+1];
      ngbs[ind*NGBSCOLS+1].p=&d[ind+s1+1];
      ngbs[ind*NGBSCOLS+2].p=&d[ind+1];
      ngbs[ind*NGBSCOLS+3].p=&d[ind+s1];
      ngbs[ind*NGBSCOLS+4].p=&d[ind-s1];
      SIDEFILL
    }

  for(i=1;i<s0-1;i++)
    {
      ind=(i+1)*s1-1;		/* Right */
      ngbs[ind*NGBSCOLS].p  =&d[ind-s1-1];
      ngbs[ind*NGBSCOLS+1].p=&d[ind+s1-1];
      ngbs[ind*NGBSCOLS+2].p=&d[ind-1];
      ngbs[ind*NGBSCOLS+3].p=&d[ind+s1];
      ngbs[ind*NGBSCOLS+4].p=&d[ind-s1];
      SIDEFILL
    }

  for(i=1;i<s0-1;i++)
    for(j=1;j<s1-1;j++)
      {
        ind=i*s1+j;		/* Body */
	ngbs[ind*NGBSCOLS].p  =&d[ind-s1-1];
	ngbs[ind*NGBSCOLS+1].p=&d[ind-s1+1];
	ngbs[ind*NGBSCOLS+2].p=&d[ind+s1-1];
	ngbs[ind*NGBSCOLS+3].p=&d[ind+s1+1];
	ngbs[ind*NGBSCOLS+4].p=&d[ind+1];
	ngbs[ind*NGBSCOLS+5].p=&d[ind-1];
	ngbs[ind*NGBSCOLS+6].p=&d[ind+s1];
	ngbs[ind*NGBSCOLS+7].p=&d[ind-s1];
	d[ind].ngb8=&ngbs[ind*NGBSCOLS];
	d[ind].ngb4=&ngbs[ind*NGBSCOLS+4];
      }

  for(i=0;i<s0*s1;i++)
    d[i].ngb=d[i].ngb8;

  *D=d;
}




















/****************************************************************
 *****************     Separate objects      ********************
 ****************************************************************/
void
freeobjpixarray(struct pix ***objs, size_t size)
{
  size_t i;
  for(i=0;i<size;i++)
    if(objs[i]!=NULL)
      free(objs[i]);
  free(objs);
}





/* A labeled image, with all connected objects labeled, and its
    corresponding graph are input. An array of output graphs are
    generated, with pointers to the pixels in the image that belong to
    each object. The background (with label zero) will also be
    considered here as the first object.

    Note that numlabs (as output from BF_concmp()) is one more than
    the largest label in the image.

    Input: D: Array of "struct pix" with "size" elements.

    The final element of "obj" is not NULL, so if it is necessary to
    remove some objects, they can safely be replaced by a NULL.

    Output: |: array element border (horizontal) or 
               directional flash (vertical)
            *: pointer. 
            p: struct pixel.

        obj: *
             |
            |*|*|*|...     (numlabs elements).
             | |
             | ->|*|*|*|...|NULL| ((num pix in this label)+1 elements)
	     --->|*|*|*|...|NULL| ((num pix in this label)+1 elements)
	          | |
	          | |
	     ------ |
             | ------     
             | |
        D:  |p|p|p|p|p|...("size" elements.)

    We allocate space for (3) here with malloc. But the caller has no
    idea what this address is, so the caller has to be given (4).*/
void
separateobjpixarray(struct pix *D, long *lab, size_t size,
		    size_t numlabs, struct pix ****out)
{
  long lii;
  struct pix ***objs;
  size_t  i, j, ii, *areas, *counters;

  /* Note that all the elements will be initialized with NULL so if an
     ID doesn't exist in the image, it doesn't cause any problems.  */
  objs=calloc(numlabs, sizeof *objs);
  assert(objs!=NULL);

  areas=calloc(numlabs, sizeof *areas);
  assert(areas!=NULL);

  /* counters[] will hold the position in each object's array that we
     have reached so far. */
  counters=calloc(numlabs, sizeof *areas);
  assert(counters!=NULL);

  for(i=0;i<size;i++)
    areas[ lab[D[i].i] ]++;

  /* The labels might not be continues, some objects might have been
     removed prior to feeding the labeled array into this function. If
     this has happened for some of the labeles, just use the initial
     NULL pointer for them so other functions know not to bother with
     this pointer.. */
  for(i=0;i<numlabs;i++)
    if(areas[i]>0)
      {
	objs[i]=malloc( (areas[i]+1) * sizeof **objs);
	assert(objs[i]!=NULL);
	objs[i][areas[i]]=NULL;
      }    

  for(i=0;i<size;i++)
    {
      ii=lii=lab[D[i].i];

      /* Make sure all the neighbors are part of the objects.  */
      for(j=0;D[i].ngb[j].p!=NULL;j++)
	if(lab[D[i].ngb[j].p->i]!=lii) 
	  D[i].ngb[j].ignore=1;

      objs[ii][counters[ii]]=&D[i];
      counters[ii]++;
    }
   
  free(counters);
  free(areas);   
  *out=objs;
}
