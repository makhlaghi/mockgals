/*********************************************************************
mockgals - Make mock astronomical profiles (galaxy, star, ...) 
           in a FITS file

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

mockgals is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mockgals is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mockgals. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#ifndef CONVOLVE_H
#define CONVOLVE_H

#define PADDEDCHECK 0		/* View the padded image. */
#define DFTCHECK    0		/* View all dft steps. */

/* Function declarations: */
void
convolve(float *f, size_t fs0, size_t fs1, 
         float *h, size_t hs0, size_t hs1,
	 float **o);

#endif
