/*********************************************************************
forqsort - Functions to be fed into qsort for various types.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

forqsort is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

forqsort is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with forqsort. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/

#ifndef FORQSORT_H
#define FORQSORT_H

int 
intdecreasing(const void * a, const void * b);

int 
intincreasing(const void * a, const void * b);

int
floatdecreasing(const void * a, const void * b);

int
floatincreasing(const void * a, const void * b);

int
doubledecreasing(const void * a, const void * b);

int
doubleincreasing(const void * a, const void * b);

#endif
