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
#ifndef UI_H
#define UI_H

#include <time.h>
#include <sys/time.h>




/* Functions called by argpparser.h*/
void
floatl0(char *optarg, float *var, char *lo, char so);

void
floatel0(char *optarg, float *var, char *lo, char so);

void
sizetlzero(char *optarg, size_t *var, char *lo, char so);

void
anyfloat(char *optarg, float *var, char *lo, char so);






/* Functions called by main.c */
void
setparams(struct mockparams *p, int argc, char *argv[],
	  time_t *rawtime);

void
freeandreporttime(struct mockparams *p, struct timeval *t0);


#endif
