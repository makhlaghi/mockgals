/*********************************************************************
fitsarrayvv - FITS image to C array and vice versa.

Copyright (C) 2013, 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

fitsarrayvv is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fitsarrayvv is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fitsarrayvv. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/

#ifndef FITSMATRIX_H
#define FITSMATRIX_H

#include "fitsio.h"

/* This structure is used to keep any additional keywords
   to add to the final output FITS image:                      */
struct keywords
{
    size_t    num_f;   /* The number of float keywords         */
    char  **names_f;   /* The names of each keyword            */
    float *values_f;   /* The float values of float keyword    */
    char **comments_f; /* Comments for the float keywords      */
    size_t    num_s;   /* The number of string keywords        */ 
    char  **names_s;   /* The names of each keyword            */
    char **values_s;   /* The string values of string keywords */
    char **comments_s; /* Comments for the string keywords     */
};



void
numextinfits(char *inname, int *numext);

void
fits_to_array(char *fits_name, int exten, int *bitpix, 
		void **array, size_t *size1, size_t *size2);

void
array_to_fits(char *fits_name, struct keywords *keys, char *EXTname, 
		int bitpix, void *array, size_t s0, size_t s1);




#endif
