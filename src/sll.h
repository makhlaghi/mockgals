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

#endif
