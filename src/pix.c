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

#include "pix.h"



/****************************************************************
 *********         Read image to pix graph         **************
 ****************************************************************/
/* Used in imgngbs(), similar for all corner pixels. */
#define CORNERFILL  n[ind*NGBSCOLS  ]=2;                \
  n[ind*NGBSCOLS+4]=-1;

/* Used in imgngbs(), similar for all side pixels. */
#define SIDEFILL  n[ind*NGBSCOLS  ]=3;                \
  n[ind*NGBSCOLS+6]=-1;


/* Convert an array into a graph. 

   Input: 
     s0: its zeroth axis size.  
     s1: its first axis size.  
   Output:
     ngbs: pointer to a size_t array.

   Basic idea: An array of size_t is created that will keep the index
   of the neighbors of all the objects, it can function on 4 or 8
   connectivity.

   Assuming the total number of pixels are 'size=s0*s1', then the ngbs
   array is going to have 'size' rows and 10 columns, the columns
   are. The total number of columns is kept in the NGBSCOLS macro
   defined in pix.h.

   Column 0:      Where the 4 connected neighbors begin.
   Column [1-8]:  The index of neighbors, or -1.
   Column 9:      -1.

   Note that the 8 connected neighbors will always begin from the
   first column.
  
   The number of neighbors is not a priori known, that is why we use a
   concept similar to stings in C, where the last element is marked
   with an impossible character. In this case, the terminating
   character is -1 (=18446744073709551615 in size_t on a 64bit
   system!).
 */
#define IMGTOPIX_SEENGBTABLE 0
void
imgngbs(size_t s0, size_t s1, size_t **ngbs)
{
  size_t ind, i, j;
  size_t *n, numngb; 

  assert( (n=malloc(s0*s1*NGBSCOLS*sizeof *n))!=NULL );
  
  ind=0;			                 /* Bottom Right */
  n[ind*NGBSCOLS+1]=ind+s1+1;		         /* 8-connected. */
  n[ind*NGBSCOLS+2]=ind+1;			 /* 4-connected. */
  n[ind*NGBSCOLS+3]=ind+s1;			 /* 4-connected. */
  CORNERFILL;

  ind=s1-1;			                 /* Bottom Left */
  n[ind*NGBSCOLS+1]=ind+s1-1;      
  n[ind*NGBSCOLS+2]=ind-1;
  n[ind*NGBSCOLS+3]=ind+s1;
  CORNERFILL;

  ind=(s0-1)*s1;		                 /* Top Right */
  n[ind*NGBSCOLS+1]=ind-s1+1;
  n[ind*NGBSCOLS+2]=ind+1;
  n[ind*NGBSCOLS+3]=ind-s1;
  CORNERFILL;
  
  ind=s0*s1-1;			                 /* Top Left */
  n[ind*NGBSCOLS+1]=ind-s1-1;
  n[ind*NGBSCOLS+2]=ind-1;
  n[ind*NGBSCOLS+3]=ind-s1;
  CORNERFILL;
  
  for(j=1;j<s1-1;j++)
    {
      ind=j;			                 /* Bottom */
      n[ind*NGBSCOLS+1]=ind+s1+1;		 /* 8-connected. */
      n[ind*NGBSCOLS+2]=ind+s1-1;		 /* 8-connected. */
      n[ind*NGBSCOLS+3]=ind+1;		         /* 4-connected. */
      n[ind*NGBSCOLS+4]=ind-1;		         /* 4-connected. */
      n[ind*NGBSCOLS+5]=ind+s1;		         /* 4-connected. */
      SIDEFILL;
    }

  for(j=1;j<s1-1;j++)
    {
      ind=(s0-1)*s1+j;	 	                 /* Top */
      n[ind*NGBSCOLS+1]=ind-s1+1;
      n[ind*NGBSCOLS+2]=ind-s1-1;
      n[ind*NGBSCOLS+3]=ind+1;
      n[ind*NGBSCOLS+4]=ind-1;
      n[ind*NGBSCOLS+5]=ind-s1;
      SIDEFILL;
    }

  for(i=1;i<s0-1;i++)
    {
      ind=i*s1;			                 /* Left */
      n[ind*NGBSCOLS+1]=ind-s1+1;
      n[ind*NGBSCOLS+2]=ind+s1+1;
      n[ind*NGBSCOLS+3]=ind+1;
      n[ind*NGBSCOLS+4]=ind+s1;
      n[ind*NGBSCOLS+5]=ind-s1;
      SIDEFILL;
    }

  for(i=1;i<s0-1;i++)
    {
      ind=(i+1)*s1-1;	 	                 /* Right */
      n[ind*NGBSCOLS+1]=ind-s1-1;
      n[ind*NGBSCOLS+2]=ind+s1-1;
      n[ind*NGBSCOLS+3]=ind-1;
      n[ind*NGBSCOLS+4]=ind+s1;
      n[ind*NGBSCOLS+5]=ind-s1;
      SIDEFILL;
    }

  for(i=1;i<s0-1;i++)
    for(j=1;j<s1-1;j++)
      {
        ind=i*s1+j;		                 /* Body */
	n[ind*NGBSCOLS+1]=ind-s1-1;		 /* 8-connected. */
	n[ind*NGBSCOLS+2]=ind-s1+1;		 /* 8-connected. */
	n[ind*NGBSCOLS+3]=ind+s1-1;		 /* 8-connected. */
	n[ind*NGBSCOLS+4]=ind+s1+1;		 /* 8-connected. */
	n[ind*NGBSCOLS+5]=ind+1;		 /* 4-connected. */
	n[ind*NGBSCOLS+6]=ind-1;		 /* 4-connected. */
	n[ind*NGBSCOLS+7]=ind+s1;		 /* 4-connected. */
	n[ind*NGBSCOLS+8]=ind-s1;		 /* 4-connected. */

	n[ind*NGBSCOLS]=5;
	n[ind*NGBSCOLS+9]=-1;
      }

  if(IMGTOPIX_SEENGBTABLE)
    {
      for(i=0;i<s0;i++)
	for(j=0;j<s1;j++)
	  {
	    printf("%lu: ", i*s1+j);
	    numngb=n[(i*s1+j)*NGBSCOLS];;
	    do
	      printf("%lu, ", n[(i*s1+j)*NGBSCOLS+numngb]);
	    while(n[(i*s1+j)*NGBSCOLS+ ++numngb]!=NONINDEX);
	    printf("\b\b.\n");
	  }
    }

  *ngbs=n;
}
