/*********************************************************************
arraymanip - Simple functions on arrays.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

arraymanip is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

arraymanip is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with arraymanip. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#ifndef ARRAYMANIP_H
#define ARRAYMANIP_H

void
floatarrcpy(float *in, size_t size, float **out);

void
uchararrcpy(unsigned char *in, size_t size, unsigned char *out);

float 
floatarrsum(float *in, size_t size);

float 
floatarrsumsquared(float *in, size_t size);

float 
floatarrsummask(float *in, long *mask, size_t size);

void
floatarrmwith(float *in, size_t size, float a);

void
floatarrswith(float *in, size_t size, float a);

void
floatsetbelowtozero(float *in, size_t size, float min);

void
floatsetabovetomax(float *in, size_t size, float max);

void
floatshrinkarray(float **in, int size1, int size2,
        int x1, int y1, int x2, int y2);

void
floatvmerge(float *a, float *b, size_t size, float **out);

#endif 
