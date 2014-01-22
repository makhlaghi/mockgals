/*********************************************************************
statistics - Library of statistical functions.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

statistics is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

statistics is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with statistics. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#ifndef STATS_H
#define STATS_H

void
fave(float *in, size_t size, float *ave);

void
fstd(float *in, size_t size, float *std);

void
histogram(float *in, size_t size, size_t numbins, 
        float min, float max, float **outhist);

void
printhist(float *in, char *filename, 
        size_t numrows, size_t numcols);

void
sigmaclip_converge(double *orderedarray, size_t num_elem, 
        double sigma_multiple, double accuracy,
        double *ave, double *med, double *std);

void
sigmaclip_certainnum(double *orderedarray, size_t num_elem, 
        double sigma_multiple, size_t numtimes,
        double *ave, double *med, double *std);

void
valuefromquantile(float *data, size_t size, 
        float quant, float *quantflux);

void
quantilefromvalue(float *data, size_t size, 
        float *quant, float quantflux);

#endif
