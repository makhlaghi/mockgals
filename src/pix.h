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

#define NGBSCOLS 10
#define NONINDEX (size_t)(-1)

void
imgngbs(size_t s0, size_t s1, size_t **ngbs);


#endif
