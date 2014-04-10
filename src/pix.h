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

#ifndef PIX_H
#define PIX_H

#define NGBSCOLS 9

struct pix;

/* In some cases (for example when separating objects) we want some of
   the neighbors (in the example, those that are not part of the
   object) to be ignored. */
struct pixlist
{
  int ignore;
  struct pix *p;
  struct pixlist *next;
};

/* Similar to pixlist but will be placed in an array. So no next is
   necessary.  */
struct ngbarrpix
{
  unsigned char ignore;
  struct pix *p;
};

/* We have 4 and 8 connectivity of the pixels.  The ngb4 linked list
   is a subset of the ngb8 linked list. The reason I have added ngb is
   that the desired connectivity might not be known before hand. So
   the user can easily set ngb to ngb4 or ngb8 based on the
   connectivity they want ontop of their function and the rest of the
   job is really easy. */
struct pix
{
  size_t i;
  float v;
  struct ngbarrpix *ngb4;
  struct ngbarrpix *ngb8;
  struct ngbarrpix *ngb;
};

void
addtopixlist(struct pixlist **queue, struct pix *p);

void
addtoanotherpixlist(struct pixlist **newq, 
        struct pixlist *oldq, struct pix *p);

void
popfrompixlist(struct pixlist **queue, struct pix **p);

void
freepixlist(struct pixlist *queue);

void
printpixlist(struct pixlist *queue);



void 
printpixarray(struct pix *D, size_t size);

void
setngb(struct pix *D, size_t size, int con_type);

void
addpixval(struct pix *D, float *img, size_t size);

void
resetignores(struct pix *D, size_t size);

void 
freepixarray(struct pix *D);



int
pix_v_increasing(const void *a, const void *b);

int
pix_v_decreasing(const void *a, const void *b);

int
pix_i_increasing(const void *a, const void *b);









void
imgtopix(size_t s0, size_t s1, struct pix **D);





void
freeobjpixarray(struct pix ***objs, size_t size);

void
separateobjpixarray(struct pix *D, long *lab, size_t size,
		    size_t numlabs, struct pix ****out);

#endif
