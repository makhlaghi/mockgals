/*********************************************************************
arraymanip - Simple functions on arrays.

Copyright (C) 2014 Mohammad Akhlaghi 
Tohoku University Astronomical Institute, Sendai, Japan.  
http://astr.tohoku.ac.jp/~akhlaghi/

arraymanip is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

arraymanip is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with arraymanip. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#ifndef ARRAYMANIP_H
#define ARRAYMANIP_H

void
longinit(long *in, size_t size, const long v);

void
ucharinit(unsigned char *in, size_t size, 
        const unsigned char v);

void
floatinit(float *in, size_t size, const float v);




void 
changefvalueinarray(float *in, size_t size, 
        float from, float to);

void
changelvalueinarray(long *in, size_t size, 
        long from, long to);

void
floatinvert(float *in, size_t size, float base);



void 
truncfarray(float *in, size_t size, float thresh);

void 
setabovetozerofarray(float *in, size_t size, float thresh);

void
floatsetbelowtozero(float *in, size_t size, float min);

void
floatsetbelowtomin(float *in, size_t size, float min);

void
floatsetabovetomax(float *in, size_t size, float max);




void
convertanytofloat(void *in, size_t size, int bitpix, float **out,
		  int byte_code, int short_code, int long_code,
		  int float_code, int double_code);

void
convertftd(float *f, size_t size, double **d);





void
removemasked(float **in, size_t size, 
        unsigned char *mask, size_t *nsize);




void
floatarrcpynomalloc(float *in, size_t size, float **out);

void
floatarrcpy(float *in, size_t size, float **out);

void
floatarrcpymask(float *in, size_t size, 
        unsigned char *mask, size_t *nsize, 
        float **out);

void
uchararrcpy(unsigned char *in, size_t size, unsigned char *out);

void
longarrcpy(long *in, size_t size, long *out);





void
maskfarray(float *in, unsigned char *mask, size_t size, 
        unsigned char f1_b0);

void
masklfarray(float *in, long *mask, size_t size, long f1_b0);







void
floatarrmwith(float *in, size_t size, float a);

void
floatarrswith(float *in, size_t size, float a);




void
checkifinarray(int *x1, int *y1, int *x2, int *y2, 
	       int s0, int s1);

void
floatshrinkarray(float **in, int size1, int size2,
		 int x1, int y1, int x2, int y2);

void
floatshrinkarraytonew(float *in, int size1, int size2,
		      int x1, int y1, int x2, int y2, 
		      float **out);

void
longshrinkarraytonew(long *in, int size1, int size2,
		      int x1, int y1, int x2, int y2, 
		     long **out);




void
floatvmerge(float *a, float *b, size_t numrows, float **out);




void
printfarray(float *array, size_t s0, size_t s1, 
	    char *comment, char *filename, int decimals, 
	    char fmttype);

#endif 
