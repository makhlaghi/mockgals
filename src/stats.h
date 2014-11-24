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
#ifndef STATS_H
#define STATS_H

#define MINFD -1e30
#define MAXFD 1e30


void
floatmin(float *in, size_t size, float *min);

void
floatmax(float *in, size_t size, float *max);

void
fminmax(float *in, size_t size, float *min, 
        float *max);

void
dmax_withindex(double *in, size_t size, 
        double *max, size_t *index);

void
fmax_withindex(float *in, size_t size, 
        float *max, size_t *index);

void
dmin_withindex(double *in, size_t size, 
        double *min, size_t *index);

void
fmin_withindex(float *in, size_t size, 
        float *min, size_t *index);






float 
floatsum(float *in, size_t size);

float 
floatsumsquared(float *in, size_t size);

float 
floatsummask(float *in, unsigned char *mask, 
        size_t size, size_t *nsize);

float 
floatsummaskl(float *in, long *mask, 
        size_t size, size_t *nsize);

float 
floatsumsquaredmask(float *in, unsigned char *mask, 
        size_t size, size_t *nsize);

float 
floatsumsquaredmaskl(float *in, long *mask, 
        size_t size, size_t *nsize);






void
fave(float *in, size_t size, float *ave, unsigned char *mask);

void
favel(float *in, size_t size, float *ave, long *mask);

void
favestd(float *in, size_t size, float *ave, float *std, 
    unsigned char *mask);

void
favestdl(float *in, size_t size, float *ave, float *std, 
    long *mask);




void
printhists(float *in, char *filename, 
	   size_t numrows, size_t numcols);

void
histogram(float *in, size_t size, size_t numbins, 
	  float *omin, float *omax, float **outhist, 
	  int abinonzero, float minq, int n01);

void
histogram_save(float *arr, unsigned char *mask, size_t size, 
		 size_t numbins, char *outname, float minq, int n01);

float
savetwohists(float *in1, float *in2, size_t size,
	     size_t numbins, char *filename, 
	     int abinonzero, float minq, int n01);


void
cumulativefp(float *in, size_t size, size_t numbins, 
	     float *omin, float *omax, float **outcum,
	     float minq, int n01);

void
cumulativefp_save(float *arr, unsigned char *mask, size_t size,
		  size_t numbins, char *outname, float minq, 
		  int n01, float max);

void
savetwocfp(float *in1, float *in2, size_t size, size_t numbins, 
	   char *filename, float minq, int n01);








void
sigmaclip_converge(float *array, int o1_n0, size_t num_elem, 
		   float sigma_multiple, float accuracy,
		   float *ave, float *med, float *std);

void
sigmaclip_certainnum(float *array, int o1_n0, size_t num_elem, 
		     float sigma_multiple, size_t numtimes, 
		     float *ave, float *med, float *std);





void
valuefromquantile(float *data, size_t size, 
        float quant, float *quantflux, unsigned char *mask);

void
multivaluefromquantile(float *data, size_t size, float *quants, 
		       float *quantfluxs, size_t numquants, 
		       unsigned char *mask);

void
valuefromquantile_nocopy(float *data, size_t size, 
        float quant, float *quantflux, unsigned char *mask);

void
quantilefromvalue(float *data, size_t size, 
        float *quant, float quantflux, unsigned char *mask);


void
quantilefromvalue_nocopy(float *data, size_t size, 
        float *quant, float quantflux, unsigned char *mask);

#endif
