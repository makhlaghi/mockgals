/*********************************************************************
mockgals - Detect and deblend objects.

Make any number of mock profiles in an array.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

mockgals is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mockgals is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mockgals. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#include "attaavv.h"
#include "mock.h"

/* Check to two size paramters */
void
checksize(char *optarg, size_t *var, int opt)
{
  long tmp;
  char *tailptr;
  tmp=strtol(optarg, &tailptr, 0);
  if(tmp<0)
    {
      printf("\n\n Error: argument to -%c ", opt); 
      printf("should be positive\n\n");
      exit(1);
    }
  *var=tmp;  
}





void
printmockgalshelp()
{
  printf("\n\nThis is mockgals version 0.01:\n\n\n");
  printf("The options to run mockgals:\n\n");
  printf(" No options:\n\tMake 50 random profiles"); 
  printf("with default parametrs\n\n");
  printf(" -h:\n\tPrint this help message.\n\n");
  printf(" -x INTIGER:\n\tThe NAXIS0 size of the ");
  printf("output FITS image.\n");
  printf("\tdefault: 201 pixels\n\n");
  printf(" -y INTIGER:\n\tThe NAXIS1 size of the ");
  printf("output FITS image.\n");
  printf("\tdefault: 201 pixels\n\n");
  printf(" -i FILENAME:\n\tInput ASCII table with 9 columns.\n");
  printf("\tdefault: mockinfo.txt\n\n");
  printf(" -o FILENAME:\n\tOutput FITS image name\n");
  printf("\tdefault: mock.fits\n\n");
  printf(" -s FLOAT:\n\tBackground value of the image.\n");
  printf("\tdefault: 10000.0\n\n");
  printf(" -a FLOAT:\n\tPSF (Moffat function) FWHM.\n");
  printf("\tdefault: 3.0\n\n");
  printf(" -b FLOAT:\n\tPSF (Moffat function) beta.\n");
  printf("\tdefault: 3.0\n\n");
  printf(" -v:\n\tView the PSF used, no argument necessary.\n");
  printf("\tBy default it is off (PSF will not be saved).\n\n");
  exit(0);
}




int
main(int argc, char *argv[])
{
  float *img;
  FILE *tmpfile;
  double *params=NULL;
  struct ArrayInfo intable;
  int c, am=0, im=0, om=0, vpsf=0;
  float psf_p1=3, psf_p2=3, sky=10000;
  size_t size1=201, size2=201, nummock=50;
  char *tailptr, *infoname=NULL, *outname=NULL;

  while( (c=getopt(argc, argv, "vhx:y:i:o:a:b:s:")) != -1 )
    switch(c)
      {
      case 'x':			/* NAXIS1 value of output image. */
	checksize(optarg, &size1, c);
	break;
      case 'y':			/* NAXIS2 value of output image. */
	checksize(optarg, &size2, c);
	break;
      case 'i':			/* Input table. */
	infoname=optarg;
	break;
      case 'o': 		/* Output fits name. */
	outname=optarg;
	break;
      case 's':			/* Sky value */
	sky=strtof(optarg, &tailptr);
	break;
      case 'a':			/* First PSF parameter. */
	psf_p1=strtof(optarg, &tailptr);
	break;
      case 'b':			/* Second PSF parameter. */
	psf_p2=strtof(optarg, &tailptr);
	break;
      case 'v':
	vpsf=1;
      case 'h':
	printmockgalshelp();
      case '?':
	return 1;
      default:
	abort();
      }

  if(infoname==NULL)		/* No input table specified: */
    {
      im=1;
      assert( (infoname=malloc(20*sizeof *infoname)) != NULL);
      strcpy(infoname, "mockinfo.txt");
    }
  
  if(outname==NULL)
    {
      om=1;
      assert( (outname=malloc(20*sizeof *outname)) != NULL);
      strcpy(outname, "mock.fits");
    }

  /* If the input table exists, read it, if not, just make a file with
     that input name. */
  if ((tmpfile = fopen(infoname, "r")) != NULL) 
    {
      am=1;
      fclose(tmpfile);
      readasciitable(infoname, &intable);
      nummock=intable.s0;
      params=intable.d;
      if(intable.s1!=10)
	{
	  printf("\n\n\tERROR: %s doesn't have 9 columns.\n", argv[1]);
	  printf("\t\tAborted\n\n");
	  exit(EXIT_FAILURE);
	}
    }

  mockimg(size1, size2, sky, nummock, params, psf_p1, 
	  psf_p2, vpsf, &img, outname, infoname);
  
  free(img);
  if(am) freeasciitable(&intable);
  if(im) free(infoname);
  if(om) free(outname);
  return 0;
}
