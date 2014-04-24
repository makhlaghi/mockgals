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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "fitsarrayvv.h"



/*************************************************************
 ***********        Find number of extentions      ***********
 *************************************************************/
void
numextinfits(char *inname, int *numext)
{
  int status=0;
  fitsfile *fptr;

  fits_open_file(&fptr, inname, READONLY, &status);  
  fits_report_error(stderr, status);

  fits_get_num_hdus(fptr, numext, &status);
  fits_report_error(stderr, status);

  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);
}




















/*************************************************************
 ***********         FITS to array functions:      ***********
 *************************************************************/
/* Allocate space for the various types */
void
f2a_readalloc(void **array, size_t size, int bitpix, 
	      int *datatype, char *ffname)
{
  short **s;  long **l; unsigned char **u;
  float **f; double **d; /* long long **ll; */

  switch(bitpix)
    {
    case BYTE_IMG:
      *datatype=TBYTE;
      u=(unsigned char **)array;
      assert( ( *u=malloc(size* sizeof**u) ) !=NULL);
      break;

    case SHORT_IMG:
      *datatype=TSHORT;
      s=(short **)array;
      assert( ( *s=malloc(size* sizeof**s) ) !=NULL);
      break;

    case LONG_IMG:
      *datatype=TLONG;
      l=(long **)array;
      assert( ( *l=malloc(size* sizeof**l) ) !=NULL);
      break;

    case FLOAT_IMG:
      *datatype=TFLOAT;
      f=(float **)array;
      assert( ( *f=malloc(size* sizeof**f) )!=NULL);
      break;

    case DOUBLE_IMG:
      *datatype=TDOUBLE;
      d=(double **)array;
      assert( ( *d=malloc(size* sizeof**d) ) !=NULL);
      break;

    case LONGLONG_IMG:
      /*
       *datatype=TLONGLONG;
      ll=(long long **)array;
      assert( ( *ll=malloc(size* sizeof**ll) ) !=NULL);
      */
      printf("\n\n%s. BITPIX=%d (long long) Not supported!\n\n",
	     ffname, bitpix);
      exit(EXIT_FAILURE);

    default:
      printf("\n\n%s. BITPIX=%d, Not recognized!\n\n",
	     ffname, bitpix);
      exit(EXIT_FAILURE);
    }
}





/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array.  */
void 
fits_to_array(char *fits_name, int exten, int *bitpix, 
        void **array, size_t *s0, size_t *s1)
{
  long b;
  int datatype;
  fitsfile *fptr;
  float *nulval=NULL;
  long fpixel[2]={1,1}, naxes[2];
  int status=0, anynul=0, maxdim=2;
  char *ffname, comment[FLEN_COMMENT];
  
  /* Add extension to the fits_name: */
  ffname=malloc(FLEN_FILENAME*sizeof(char));   
  sprintf(ffname, "%s[%d]", fits_name, exten);

  fits_open_file(&fptr, ffname, READONLY, &status);  
  fits_report_error(stderr, status);

  fits_read_key(fptr, TLONG, "BITPIX", &b, comment, &status);
  fits_report_error(stderr, status);
  *bitpix=b;

  fits_get_img_size(fptr, maxdim, naxes, &status);
  fits_report_error(stderr, status);
  *s0=naxes[1];
  *s1=naxes[0];

  f2a_readalloc(array, *s0 * *s1, *bitpix, &datatype, ffname);

  fits_read_pix(fptr, datatype, fpixel,
                naxes[0]*naxes[1], nulval,
                *array, &anynul, &status);
  fits_report_error(stderr, status);

  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);

  free(ffname);
}




















/*************************************************************
 ***********         Array to FITS functions:      ***********
 *************************************************************/
void
a2f_setdatatype(int bitpix, int *datatype)
{
  switch(bitpix)
    {
    case BYTE_IMG:
      *datatype=TBYTE;
      break;
    case SHORT_IMG:
      *datatype=TSHORT;
      break;
    case LONG_IMG:
      *datatype=TLONG;
      break;
    case FLOAT_IMG:
      *datatype=TFLOAT;
      break;
    case DOUBLE_IMG:
      *datatype=TDOUBLE;
      break;
    case LONGLONG_IMG:
      *datatype=TLONGLONG;
      break;
    default:
      printf("\n\nInput bixpix=%d, Not recognized!\n\n", bitpix);
      exit(EXIT_FAILURE);
    }
}





void
array_to_fits(char *fits_name, struct keywords *keys, char *EXTname, 
	      int bitpix, void *array, size_t s0, size_t s1)
{
  size_t i;
  fitsfile *fptr;
  char comm[100];
  int status=0, datatype;
  long fpixel=1, naxis=2, nelements, naxes[2];

  a2f_setdatatype(bitpix, &datatype);

  naxes[1]=s0;
  naxes[0]=s1;
  nelements=naxes[0]*naxes[1];

  if(access(fits_name,F_OK) != -1 )
    fits_open_file(&fptr,fits_name, READWRITE, &status);
  else 
    fits_create_file(&fptr, fits_name, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  fits_write_key(fptr, TSTRING, "EXTNAME", EXTname, "", &status);
  if(keys!=NULL)
    {
      for(i=0;i<keys->num_f;i++)
	fits_write_key(fptr, TFLOAT, keys->names_f[i], 
		       &keys->values_f[i], keys->comments_f[i], 
		       &status);
      for(i=0;i<keys->num_s;i++)
	fits_write_key(fptr, TSTRING, keys->names_s[i], 
		       keys->values_s[i], keys->comments_s[i], 
		       &status);
    }
  sprintf(comm, "Created with CFITSIO v%.2f", CFITSIO_VERSION);
  fits_write_comment(fptr, comm, &status);

  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);
}
