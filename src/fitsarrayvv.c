/*********************************************************************
fitsarrayvv - FITS image to C array and vice versa.

Copyright (C) 2014 Mohammad Akhlaghi
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

#include "fitsio.h"
#include "fitsarrayvv.h"




















/*************************************************************
 ***********         FITS to array functions:      ***********
 *************************************************************/
#define FITSTOARRAYTOP   int datatype;                \
  fitsfile *fptr;                                     \
  float *nulval=NULL;                                 \
  long fpixel[2]={1,1}, naxes[2];                     \
  int status=0, anynul=0, maxdim=2;                   \





#define FITSTOARRAYBODY fits_open_file(&fptr,         \
              fits_name, READONLY, &status);          \
  fits_report_error(stderr, status);                  \
  fits_get_img_size(fptr, maxdim, naxes, &status);    \
  fits_report_error(stderr, status);                  \
  *size1=naxes[1];                                    \
  *size2=naxes[0];                                    \
  *array=malloc(sizeof **array * *size1 * *size2);    \
  assert(*array!=NULL);                               \
  fits_read_pix(fptr, datatype, fpixel,               \
                naxes[0]*naxes[1], nulval,            \
                *array, &anynul, &status);            \
  fits_report_error(stderr, status);                  \
  fits_close_file(fptr, &status);                     \
  fits_report_error(stderr, status);                  \





void 
fits_to_array_float(char *fits_name, void **array0, 
        size_t *size1, size_t *size2)
{
  FITSTOARRAYTOP
  float **array=(float **)array0;
  datatype=TFLOAT;
  FITSTOARRAYBODY
}





void 
fits_to_array_long(char *fits_name, void **array0, 
        size_t *size1, size_t *size2)
{
  FITSTOARRAYTOP
  long **array=(long **)array0;
  datatype=TLONG;
  FITSTOARRAYBODY
}





void 
fits_to_array_short(char *fits_name, void **array0, 
        size_t *size1, size_t *size2)
{
  FITSTOARRAYTOP
  short **array=(short **)array0;
  datatype=TSHORT;
  FITSTOARRAYBODY
}






void 
fits_to_array_uchar(char *fits_name, void **array0, 
        size_t *size1, size_t *size2)
{
  FITSTOARRAYTOP
  unsigned char **array=(unsigned char **)array0;
  datatype=TBYTE;
  FITSTOARRAYBODY
}





/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array.  */
void 
fits_to_array(char *fits_name, int exten, char *fitstype, 
        void **array, size_t *size1, size_t *size2)
{
  char *finalfitsname;

  /* Add extension to the fits_name: */
  finalfitsname=malloc(strlen(fits_name)*200
		       *sizeof(char));   
  sprintf(finalfitsname, "%s[%d]", fits_name, exten);

  /* Set the output type of the final FITS image: */
  if (strcmp(fitstype, "FLOAT")==0)
    fits_to_array_float(finalfitsname, array, size1, size2);
  else if (strcmp(fitstype, "SHORT")==0)
    fits_to_array_short(finalfitsname, array, size1, size2);
  else if (strcmp(fitstype, "UCHAR")==0)
    fits_to_array_uchar(finalfitsname, array, size1, size2);
  else 
    {
      printf("\n\nfitstype should be FLOAT, SHORT, UCHAR\n\n");
      exit(EXIT_FAILURE);
    }

  free(finalfitsname);
}




















/*************************************************************
 ***********         Array to FITS functions:      ***********
 *************************************************************/
/* Convert a float array to a float FITS image. 
   Note that for a float:
   bitpix=FLOAT_IMG;        datatype=TFLOAT;*/
void
array_to_fits_float(char *fits_name, void *array0, 
        struct keywords *keys, char *EXTname, 
        size_t size1, size_t size2)
{
  size_t i;
  fitsfile *fptr;
  char comm[100];
  int status=0, bitpix, datatype;
  long fpixel=1, naxis=2, nelements, naxes[2];

  float *array=(float *)array0;
  bitpix=FLOAT_IMG;        datatype=TFLOAT;

  /* Set the naxes values: */
  naxes[1]=size1;
  naxes[0]=size2;
  nelements=naxes[0]*naxes[1];

  /* Check if the file exists, if so, open it, if not create one: */
  if(access(fits_name,F_OK) != -1 )
    fits_open_file(&fptr,fits_name, READWRITE, &status);
  else 
    fits_create_file(&fptr, fits_name, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  /* Write the keywords:*/
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

  /* Close the FITS file: */
  fits_close_file(fptr, &status);

  /* Report any possible errors: */
  fits_report_error(stderr, status);
}





/* Convert a long array to a long FITS image.
   Note that for a long:
   bitpix=LONG_IMG;        datatype=TLONG;
   to array_to_fits_float() that is why there is no 
   comments in this function.*/
void
array_to_fits_long(char *fits_name, void *array0, 
        struct keywords *keys, char *EXTname, 
        size_t size1, size_t size2)
{
  size_t i;
  fitsfile *fptr;
  char comm[100];
  int status=0, bitpix, datatype;
  long fpixel=1, naxis=2, nelements, naxes[2];

  long *array=(long *)array0;
  bitpix=LONG_IMG;        datatype=TLONG;

  naxes[1]=size1;
  naxes[0]=size2;
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




/* Convert a short array to a short FITS image.
   Note that for a short:
   bitpix=SHORT_IMG;        datatype=TSHORT;
   to array_to_fits_float() that is why there is no 
   comments in this function.  */
void
array_to_fits_short(char *fits_name, void *array0, 
        struct keywords *keys, char *EXTname, 
        size_t size1, size_t size2)
{
  size_t i;
  fitsfile *fptr;
  char comm[100];
  int status=0, bitpix, datatype;
  long fpixel=1, naxis=2, nelements, naxes[2];

  short *array=(short *)array0;
  bitpix=SHORT_IMG;        datatype=TSHORT;

  naxes[1]=size1;
  naxes[0]=size2;
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





/* Convert an unsigned char array to a byte FITS image.
   Note that for an unsigned char:
   bitpix=BYTE_IMG;        datatype=TBYTE;
   All the steps in the function are identical 
   to array_to_fits_float() that is why there is no 
   comments in this function.  */
void
array_to_fits_uchar(char *fits_name, void *array0, 
        struct keywords *keys, char *EXTname, 
        size_t size1, size_t size2)
{
  size_t i;
  fitsfile *fptr;
  char comm[100];
  int status=0, bitpix, datatype;
  long fpixel=1, naxis=2, nelements, naxes[2];

  unsigned char *array=(unsigned char *)array0;
  bitpix=BYTE_IMG;        datatype=TBYTE;

  naxes[1]=size1;
  naxes[0]=size2;
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





/* Save an array as a FITS image based on the array type.*/
void
array_to_fits(char *fits_name, struct keywords *keys, 
        char *EXTname, char *fitstype, void *array, 
        size_t size1, size_t size2)
{
  /* Set the output type of the final FITS image: */
  if (strcmp(fitstype, "FLOAT")==0)
    array_to_fits_float(fits_name, array, 
			keys, EXTname, size1, size2);
  else if (strcmp(fitstype, "UCHAR")==0)
    array_to_fits_uchar(fits_name, array, 
			keys, EXTname, size1, size2);
  else if (strcmp(fitstype, "SHORT")==0)
    array_to_fits_short(fits_name, array, 
			keys, EXTname, size1, size2);
  else if (strcmp(fitstype, "LONG")==0)
    array_to_fits_long(fits_name, array, 
		       keys, EXTname, size1, size2);
  else 
    {
      printf("\n\nfitstype should be FLOAT, SHORT, LONG, UCHAR\n\n");
      exit(EXIT_FAILURE);
    }
}
