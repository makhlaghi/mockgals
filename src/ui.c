/*********************************************************************
mockgals - Make mock astronomical profiles (galaxy, star, ...) 
           in a FITS file

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
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include "attaavv.h"
#include "mock.h"
#include "ui.h"





/****************************************************************
 *****************     Small functions used    ******************
 *********************      in main()      **********************
 ****************************************************************/
/* Set the default values for the inputs into mockimg() (mock.h). */
void
setdefaultoptions(struct mockparams *p)
{
  /* On or off options: */
  p->vhist     =0;
  p->verb      =0;
  p->vpsf      =0;
  p->vnoconv   =0;
  p->vconv     =0;

  /* Options with arguments: */
  p->s0        =201;
  p->s1        =201;
  p->sky       =10000.0f;
  p->trunc     =5;
  p->psfname   ="";
  p->psf_mg    =1;
  p->psf_p1    =5;
  p->psf_p2    =3;
  p->histmin   =-250;
  p->histmax   =700;
  p->nummock   =50;
  p->infoname  ="mockinfo.txt";
  p->outname   ="mock.fits";

  /* Internal settings: */
  p->profileparams=NULL;
  p->numppcols    =10;
}





/* If the input table exists, read it and put its pointer into
   mockparams. If it doesn't exist, then let prflprms(), which is
   is mock.c, make a random set of parameters. */
void
readinputinfo(struct mockparams *p)
{
  FILE *tmpfile;
  struct ArrayInfo intable;

  if ((tmpfile = fopen(p->infoname, "r")) != NULL) 
    {
      fclose(tmpfile);
      readasciitable(p->infoname, &intable);
      if(intable.s1!=10)
	{
	  fprintf(stderr, "\n\n\tERROR: %s not %lu columns.\n", 
		  p->infoname, p->numppcols);
	  fprintf(stderr, "\t\tAborted\n\n");
	  exit(EXIT_FAILURE);
	}
      p->nummock=intable.s0;
      p->profileparams=intable.d;
      intable.d=malloc(sizeof *(intable.d));/*freeasciitable() has */
      assert(intable.d!=NULL);		    /*something to free! */
      freeasciitable(&intable);
      if(p->verb)
	printf("- Information for %lu profile%sread from '%s'.\n\n",
	     p->nummock, p->nummock>1 ? "s " : " ", p->infoname);
    }
  else 
    {
      setprflprms(&p->profileparams, p->numppcols, p->nummock, 
		  p->s1, p->s0);
      if(p->verb)
	{
	  printf("- %lu random profiles will be made.\n", 
		 p->nummock);
	  printf("  Their info will be saved in '%s'.\n\n", 
		 p->infoname);
	}
    }
}





/* Check if the output image exists. If so, remove it. */
void
checkremoveoutimage(char *outname)
{
  FILE *tmpfile;
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
}




















/****************************************************************
 *****************        Read options:      ********************
 ****************************************************************/
/* Check if the two size paramters are positive. */
void
checksize(char *optarg, size_t *var, int opt)
{
  long tmp;
  char *tailptr;
  tmp=strtol(optarg, &tailptr, 0);
  if(tmp<=0)
    {
      printf("\n\n Error: argument to -%c ", opt); 
      printf("should be positive\n\n");
      exit(EXIT_FAILURE);
    }
  *var=tmp;  
}





/* Print the help menu. */
void
printmockgalshelp(struct mockparams *p)
{
  printf("\n\nThis is mockgals version 0.01:\n\n\n");
  printf(" No options:\n\tMake %lu random profiles\n\n", p->nummock); 
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
  printf("\tdefault: %lu pixels\n\n", p->s1);
  printf(" -y INTEGER:\n\tThe NAXIS1 size of the ");
  printf("output FITS image.\n");
  printf("\tdefault: %lu pixels\n\n", p->s0);
  printf(" -i FILENAME:\n\tInput ASCII table with 9 columns.\n");
  printf("\tdefault: '%s'\n\n", p->infoname);
  printf(" -o FILENAME:\n\tOutput FITS image name\n");
  printf("\tdefault: '%s'\n\n", p->outname);
  printf(" -s FLOAT:\n\tBackground value of the image.\n");
  printf("\tdefault: %.2f\n\n", p->sky);
  printf(" -t FLOAT:\n\tProfile truncation, a multiple of radius.\n");
  printf("\tdefault: %.2f\n\n", p->trunc);
  printf(" -u INTEGER:\n\tPSF radial function.\n");
  printf("\t1: Moffat function.\n");
  printf("\t2: Gaussian function.\n");
  printf("\tAny other input value will be changed to default.\n");
  printf("\tdefault: %d.\n", p->psf_mg);
  printf("\tdefault: %.2f\n\n", p->psf_p1);
  printf(" -a FLOAT:\n\tPSF (Moffat function) FWHM.\n");
  printf("\tdefault: %.2f\n\n", p->psf_p1);
  printf(" -b FLOAT:\n\tPSF (Moffat function) beta.\n");
  printf("\tdefault: %.2f\n\n", p->psf_p2);
  printf(" -g INTEGER:\n\tIf positive, print histogram.\n");
  printf("\tdefault: %d. Integer (argument) number of bins.\n\n",
	 p->vhist);
  printf(" -c FLOAT:\n\tHistogram minimum value \n");
  printf("\tdefault: %.2f.\n\n", p->histmin);
  printf(" -d FLOAT:\n\tHistogram maximum value \n");
  printf("\tdefault: %.2f.\n\n", p->histmax);
  exit(0);
}





/* Read all the options into the program */
void
getsaveoptions(struct mockparams *p, 
	       int argc, char *argv[])
{
  int c;
  char *tailptr; 

  while( (c=getopt(argc, argv, "pmnvhx:y:i:o:a:b:s:g:c:d:f:t:u:")) 
	 != -1 )
    switch(c)
      {
      case 'x':			/* NAXIS1 value of output image. */
	checksize(optarg, &p->s1, c);
	break;
      case 'y':			/* NAXIS2 value of output image. */
	checksize(optarg, &p->s0, c);
	break;
      case 'i':			/* Input table. */
	p->infoname=optarg;
	break;
      case 'o': 		/* Output fits name. */
	p->outname=optarg;
	break;
      case 's':			/* Sky value. */
	p->sky=strtof(optarg, &tailptr);
	break;
      case 't':			/* PSF FWHM.4 */
	p->trunc=strtof(optarg, &tailptr);
	break;
      case 'f':			/* Input PSF name */
	p->psfname=optarg;
	break;
      case 'u':			/* PSF function type */
	p->psf_mg=strtol(optarg, &tailptr, 0);
	if(p->psf_mg<1 || p->psf_mg>2) 
	  p->psf_mg=1;
	break;
      case 'a':			/* PSF FWHM.4 */
	p->psf_p1=strtof(optarg, &tailptr);
	break;
      case 'b':			/* If Moffat, the beta. */
	p->psf_p2=strtof(optarg, &tailptr);
	break;
      case 'g':			/* View histogram. */
	p->vhist=strtol(optarg, &tailptr, 0);
	if(p->vhist<=0) p->vhist=0;
	break;
      case 'c':			/* Histogram minimum. */
	p->histmin=strtof(optarg, &tailptr);
	break;
      case 'd':			/* Histogram maximum. */
	p->histmax=strtof(optarg, &tailptr);
	break;
      case 'v':			/* Verbatim mode. */
	p->verb=1;
	break;
      case 'p':			/* View PSF. */
	p->vpsf=1;
	break;
      case 'm':			/* View not convolved. */
	p->vnoconv=1;
	break;
      case 'n':			/* View no noised image. */
	p->vconv=1;
	break;
      case 'h':			/* Print help. */
	printmockgalshelp(p);
	exit(EXIT_FAILURE);
      case '?':
	fprintf(stderr, "Unknown option: '-%c'.\n\n", optopt);
	exit(EXIT_FAILURE);
      default:
	abort();
      }
}
