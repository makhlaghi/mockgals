/*********************************************************************
noisechisel - Detect and deblend objects.

Detect and deblend (segment) objects in a noisy image.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

noisechisel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

noisechisel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with noisechisel. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef CONVOLVE_H
#define CONVOLVE_H

#define DFTCHECK 0    /* 1: Check all steps of the 
                      descrete fourier transform in convolve(). */

/* Function declarations: */
void
convolve(float *f, size_t fs0, size_t fs1, 
         float *h, size_t hs0, size_t hs1,
	 float **o);

#endif
