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
#ifndef SLL_H
#define SLL_H


/******************* float: */
struct fsll
{
    float v;
    struct fsll *next;
};

void
printarrayoffsll(struct fsll **afsll, size_t num);

void
add_to_fsll(struct fsll **list, float value);

void
pop_from_fsll(struct fsll **list, float *value);

size_t
numinfsll(struct fsll *list);

void
fslltoarray(struct fsll *list, float **f, size_t *num);

void
freefsll(struct fsll *list);

void
freearrayoffsll(struct fsll **afsll, size_t num);



/******************* size_t: */
struct ssll
{
    size_t v;
    struct ssll *next;
};

void
add_to_ssll(struct ssll **list, size_t value);

void
pop_from_ssll(struct ssll **list, size_t *value);

size_t
numinssll(struct ssll *list);

void
sslltoarray(struct ssll *list, size_t **f, size_t *num);

void
freessll(struct ssll *list);


/******************* Two way size_t: */
struct tssll
{
  size_t v;
  struct tssll *next;
  struct tssll *prev;
};

void
add_to_tssll_end(struct tssll **last, size_t value);

void
pop_from_tssll_start(struct tssll **first,  size_t *value);


/******************* Ordered size_t: */
struct ossll
{
  size_t v;			/* The actual value. */
  float s;			/* The parameter to sort by. */
  struct ossll *next;
};

void
add_to_ossll(struct ossll **list, size_t value, float tosort);

void
pop_from_ossll(struct ossll **list,  size_t *value, float *sortvalue);

void
ossll_into_ssll(struct ossll *in, struct ssll **out);


/******************* Two way ordered size_t: */
struct tossll
{
  size_t v;
  float s;
  struct tossll *prev;
  struct tossll *next;
};

void
print_tossll(struct tossll *l, struct tossll *s);

void
add_to_tossll_end(struct tossll **largest, struct tossll **smallest, 
		  size_t value, float tosort);

void
pop_from_tossll_start(struct tossll **lartest, struct tossll **smallest,  
		      size_t *value, float *tosort);

void
smallest_tossll(struct tossll *largest, struct tossll **smallest);

void
tossll_into_ssll(struct tossll *in, struct ssll **out);

#endif
