/*********************************************************************
freqdomain - Frequency domain applications on C 2d arrays.

Copyright (C) 2013 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

freqdomain is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

freqdomain is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with freqdomain. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/

#ifndef FREQDOMAIN_H
#define FREQDOMAIN_H

#define CONVOLVECHECK 0    /* 1: Check all steps of convolution. */

/* Function declarations: */
void
convolve(float *f, size_t fsize1, size_t fsize2, 
         float *h, size_t hsize1, size_t hsize2);

void
convolve_function(float *f, size_t fsize1, size_t fsize2,
        float p1, float p2, float pa_d, float q, float trunc, 
        float integaccu, float s0_m1_g2, float **conv);

#endif
