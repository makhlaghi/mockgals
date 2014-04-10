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
#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sys/time.h>

#include "mock.h"
#include "attaavv.h"

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
printmockgalshelp(size_t nummock, size_t s0, size_t s1, 
		  char *dinfoname, char *doutname, float sky,
		  float psf_p1, float psf_p2, int vhist, 
		  float histmin, float histmax)
{
  printf("\n\nThis is mockgals version 0.01:\n\n\n");
  printf(" No options:\n\tMake %lu random profiles\n\n", nummock); 
  printf("\n\n###### Options with no argument:\n");
  printf("###### In default, they are all inactive.\n");
  printf(" -h:\n\tPrint this help message.\n\n");
  printf(" -v:\n\tVerbose: report steps.\n\n");
  printf(" -p:\n\tView the PSF used, no argument necessary.\n\n");
  printf(" -m:\n\tView unconvolved mock image.\n");
  printf("\tAs a prior extension to main output. \n\n");
  printf(" -n:\n\tView convolved (before adding noise) image.\n");
  printf("\tAs a prior extension to main output. \n\n");
  printf("\n\n###### Options with an argument:\n");
  printf(" -x INTEGER:\n\tThe NAXIS0 size of the ");
  printf("output FITS image.\n");
  printf("\tdefault: %lu pixels\n\n", s1);
  printf(" -y INTEGER:\n\tThe NAXIS1 size of the ");
  printf("output FITS image.\n");
  printf("\tdefault: %lu pixels\n\n", s0);
  printf(" -i FILENAME:\n\tInput ASCII table with 9 columns.\n");
  printf("\tdefault: '%s'\n\n", dinfoname);
  printf(" -o FILENAME:\n\tOutput FITS image name\n");
  printf("\tdefault: '%s'\n\n", doutname);
  printf(" -s FLOAT:\n\tBackground value of the image.\n");
  printf("\tdefault: %.2f\n\n", sky);
  printf(" -a FLOAT:\n\tPSF (Moffat function) FWHM.\n");
  printf("\tdefault: %.2f\n\n", psf_p1);
  printf(" -b FLOAT:\n\tPSF (Moffat function) beta.\n");
  printf("\tdefault: %.2f\n\n", psf_p2);
  printf(" -t INTEGER:\n\tIf positive, print histogram.\n");
  printf("\tdefault: %d. Integer (argument) number of bins.\n\n",
	 vhist);
  printf(" -c FLOAT:\n\tHistogram minimum value \n");
  printf("\tdefault: %.2f.\n\n", histmin);
  printf(" -d FLOAT:\n\tHistogram maximum value \n");
  printf("\tdefault: %.2f.\n\n", histmax);
  exit(0);
}




int
main(int argc, char *argv[])
{
  int counter=0;
  FILE *tmpfile;
  time_t rawtime;
  double *params=NULL;
  struct timeval t0, t1;
  struct ArrayInfo intable;
  char doutname[]="mock.fits";
  int vconv=0, vnoconv=0, verb=0;
  char dinfoname[]="mockinfo.txt";
  size_t s1=201, s0=201, nummock=50;
  float psf_p1=3, psf_p2=3, sky=10000;
  int c, am=0, im=0, om=0, vpsf=0, vhist=0;
  float *img, histmin=-250.0f, histmax=700.0f;
  char *tailptr, *infoname=NULL, *outname=NULL;

  time(&rawtime);
  gettimeofday(&t0, NULL);

  while( (c=getopt(argc, argv, "pmnvhx:y:i:o:a:b:s:t:c:d:")) != -1 )
    switch(c)
      {
      case 'x':			/* NAXIS1 value of output image. */
	checksize(optarg, &s1, c);
	break;
      case 'y':			/* NAXIS2 value of output image. */
	checksize(optarg, &s0, c);
	break;
      case 'i':			/* Input table. */
	infoname=optarg;
	break;
      case 'o': 		/* Output fits name. */
	outname=optarg;
	break;
      case 's':			/* Sky value. */
	sky=strtof(optarg, &tailptr);
	break;
      case 'a':			/* First PSF parameter. */
	psf_p1=strtof(optarg, &tailptr);
	break;
      case 'b':			/* Second PSF parameter. */
	psf_p2=strtof(optarg, &tailptr);
	break;
      case 't':			/* View histogram. */
	vhist=strtol(optarg, &tailptr, 0);
	if(vhist<=0) vhist=0;
	break;
      case 'c':			/* Histogram minimum. */
	histmin=strtof(optarg, &tailptr);
	break;
      case 'd':			/* Histogram maximum. */
	histmax=strtof(optarg, &tailptr);
	break;
      case 'v':			/* View PSF. */
	verb=1;
	break;
      case 'p':			/* View PSF. */
	vpsf=1;
	break;
      case 'm':			/* View not convolved. */
	vnoconv=1;
	break;
      case 'n':			/* View no noised image. */
	vconv=1;
	break;
      case 'h':			/* Print help. */
	printmockgalshelp(nummock, s0, s1, dinfoname, doutname, 
			  sky, psf_p1, psf_p2, vhist, histmin, 
			  histmax);
	return 1;
      case '?':
	fprintf(stderr, "Unknown option: '-%c'.\n\n", optopt);
	return 1;
      default:
	abort();
      }

  if(infoname==NULL)		/* No input table specified: */
    {
      im=1;
      assert( (infoname=malloc(20*sizeof *infoname)) != NULL);
      strcpy(infoname, dinfoname);
    }
  
  if(outname==NULL)
    {
      om=1;
      assert( (outname=malloc(20*sizeof *outname)) != NULL);
      strcpy(outname, doutname);
    }

  if(verb)
    {
      printf("\n\n--------------------------\n");
      printf("mockgals started on %s\n", ctime(&rawtime));
    }

  /* Check if the output image exists: */
  if ((tmpfile = fopen(outname, "r")) != NULL) 
    {
      fclose(tmpfile);
      if(unlink(outname)==-1)
	{
	  fprintf(stderr, "'%s' already exists and could", outname); 
	  fprintf(stderr, " not be removed");
	  exit(EXIT_FAILURE);
	}
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
	  fprintf(stderr, "\n\n\tERROR: %s doesn't have 9 columns.\n", 
		  argv[1]);
	  fprintf(stderr, "\t\tAborted\n\n");
	  exit(EXIT_FAILURE);
	}
      if(verb)
	printf("- Information for %lu profile%sread from '%s'.\n\n",
	     nummock, nummock>1 ? "s " : " ", infoname);
    }
  else if(verb)
    printf("- %lu random profiles will be made.\n\n", nummock);

  mockimg(s0, s1, sky, nummock, params, psf_p1, psf_p2, 
	  vpsf, vnoconv, vconv, vhist, histmin, histmax, 
	  &img, outname, infoname);

  if(verb)
    {
      if(nummock>1)
	printf("- All %lu profiles made, convolved and noised.\n\n",
	       nummock);
      else
	printf("- Profile is made, convolved and noised.\n\n");
      printf("- Profile info saved in '%s'\n\n", infoname);

      if(vnoconv)
	printf("- Not convolved profiles saved in '%s' (ext %d)\n",
	       outname, counter++);
      if(vconv)
	printf("- Convolved profiles saved in '%s' (ext %d)\n",
	       outname, counter++);
      if(counter==0)
	printf("- Final image saved in '%s'\n\n", outname);
      else
	printf("- Final image saved in '%s' (ext %d)\n\n", 
	       outname, counter);
    }
  
  free(img);
  if(am) freeasciitable(&intable);
  if(im) free(infoname);
  if(om) free(outname);

  if(verb)
    {
      gettimeofday(&t1, NULL);
      printf("mockgals finished in %.4f (seconds)\n",
	     ((double)t1.tv_sec+(double)t1.tv_usec/1e6) - 
	     ((double)t0.tv_sec+(double)t0.tv_usec/1e6));
      printf("--------------------------\n\n");
    }

  return 0;
}
