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
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include "main.h"
#include "mock.h"
#include "stats.h"
#include "attaavv.h"
#include "arraymanip.h"

#include "ui.h"			/* Needs mock.h */





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
  p->initcomments=NULL;
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
      assert(intable.d!=NULL);		    /*something to free!   */

      p->initcomments=intable.c;
      intable.c=malloc(sizeof *(intable.c));/*freeasciitable() has */
      assert(intable.c!=NULL);		    /*something to free!   */ 

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





void
printversioninfo()
{
  printf("\n\nMockGals %.1f\n", MOCKGALSVERSION);
  printf("============\n");
  printf("Make mock stars and galaxies in a FITS image.\n");
  printf("\nCopyright (C) 2014  Mohammad Akhlaghi\n");
  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to\n");
  printf("modify and redistribute it under the\n");
  printf("GNU Public License v3 or later.\n\n\n");
}




/* Print the help menu. */
void
printmockgalshelp(struct mockparams *p)
{
  printversioninfo();
  printf("No options:\n\tMake %lu random profiles\n\n", p->nummock); 

  printf("\n\n###### Options that won't run Mockgals\n");
  printf(" -h:\n\tPrint this help message.\n\n");
  printf(" -v:\n\tPrint version and copyright information.\n\n\n");

  printf("\n\n###### Options with no argument:\n");
  printf("###### In default, they are all inactive.\n");
  printf(" -e:\n\tVerbose: report all steps.\n\n");
  printf(" -p:\n\tView the PSF used, no argument necessary.\n\n");
  printf(" -m:\n\tView unconvolved mock image.\n");
  printf("\tAs a prior extension to main output. \n\n");
  printf(" -n:\n\tView convolved (before adding noise) image.\n");
  printf("\tAs a prior extension to main output. \n\n\n");

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

  printf(" -f FILENAME:\n\tInput PSF fits file.\n");
  printf("\tIf a file is specified, the other\n");
  printf("\tPSF parameters are ignored.\n\n");

  printf(" -u INTEGER:\n\tPSF radial function.\n");
  printf("\t1: Moffat function.\n");
  printf("\t2: Gaussian function.\n");
  printf("\tAny other input value will be changed to default.\n");
  printf("\tdefault: %d.\n\n", p->psf_mg);

  printf(" -a FLOAT:\n\tPSF FWHM.\n");
  printf("\tdefault: %.2f\n\n", p->psf_p1);

  printf(" -b FLOAT:\n\tPSF (Moffat function) beta.\n");
  printf("\tNot used for a Gaussian.\n");
  printf("\tdefault: %.2f\n\n", p->psf_p2);

  printf(" -g INTEGER:\n\tIf positive, print histogram.\n");
  printf("\tdefault: %d. Integer (argument) number of bins.\n\n",
	 p->vhist);

  printf(" -c FLOAT:\n\tHistogram minimum value \n");
  printf("\tdefault: %.2f.\n\n", p->histmin);

  printf(" -d FLOAT:\n\tHistogram maximum value \n");
  printf("\tdefault: %.2f.\n\n\n", p->histmax);
  exit(0);
}





/* Read all the options into the program */
void
getsaveoptions(struct mockparams *p, 
	       int argc, char *argv[])
{
  int c;
  char *tailptr; 

  while( (c=getopt(argc, argv, "pmnevhx:y:i:o:a:b:s:g:c:d:f:t:u:")) 
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
      case 'e':			/* Verbose mode. */
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
      case 'v':			/* Print version and copyright. */
	printversioninfo();
	exit(EXIT_FAILURE);
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




















/****************************************************************
 *****************  Save mock image and info ********************
 ****************************************************************/
/* If the mock image is to be saved, save the information
   of the galaxies and the actual mock image. */
void
savemockinfo(struct mockparams *p)
{
  char temp[1000];
  struct ArrayInfo ai;
  int int_cols[]={0, 1, 6,-1}, accu_cols[]={2,3,8,9,-1};
  int space[]={6,8,11}, prec[]={2,4};

  ai.c=malloc(MAXALLCOMMENTSLENGTH*sizeof(char));
  assert(ai.c!=NULL);
  ai.s0=p->nummock;
  ai.s1=p->numppcols;
  ai.d=p->profileparams;

  if(p->initcomments==NULL)
    {
      sprintf(temp, "# Properties of %lu mock profiles.\n", 
	      p->nummock);
      strcpy(ai.c, temp);
      sprintf(temp, "# The sky valued is assumed to be: %.2f\n", 
	      p->sky);
      strcat(ai.c, temp);
      sprintf(temp, "# Truncation at %.2f * radial parameter\n# \n", 
	      p->trunc);
      strcat(ai.c, temp);
      strcat(ai.c, "# 0: ID.\n");
      strcat(ai.c, "# 1: 0: Sersic, 1: Moffat, 2: Gaussian, ");
      strcat(ai.c, "3: Point.\n");
      strcat(ai.c, "# 2: X position (FITS definition).\n");
      strcat(ai.c, "# 3: Y position (FITS definition).\n");
      strcat(ai.c, "# 4: Sersic re or Moffat FWHM.\n");    
      strcat(ai.c, "# 5: Sersic n or Moffat beta.\n");    
      strcat(ai.c, "# 6: Position angle, degrees.\n");    
      strcat(ai.c, "# 7: Axis ratio.\n");    
      strcat(ai.c, "# 8: Signal to noise: ");
      strcat(ai.c, "(average profile flux-sky)/sqrt(sky).\n");    
      strcat(ai.c, "# 9: Total flux (Sky subtracted).\n\n"); 
    }
  else ai.c=p->initcomments;
  writeasciitable (p->infoname, &ai, int_cols, 
		   accu_cols, space, prec);
}





void
printmockhist(float *img, size_t size, int numbins, float histmin,
	      float histmax, float *nonoisehist)
{
  char temp[1000];
  double *dallhist;
  struct ArrayInfo ai;
  float *noisedhist, *allhist;

  int int_cols[]={1, 2, -1}, accu_cols[]={-1};
  int space[]={6,8,15}, prec[]={3,4};

  histogram(img, size, numbins, &histmin, &histmax, 
	    &noisedhist, 1, 0, 0);

  floatvmerge(nonoisehist, noisedhist, numbins+1, &allhist);
  convertftd(allhist, (numbins+1)*3, &dallhist);

  ai.c=malloc(MAXALLCOMMENTSLENGTH*sizeof(char));
  assert(ai.c!=NULL);
  ai.s0=numbins+1;
  ai.s1=3;
  ai.d=dallhist;

  sprintf(temp, "# Histogram of noised and no noised mock image.\n");
  strcpy(ai.c, temp);
  sprintf(temp, "# Range: %.3f-%.3f\n", histmin, histmax);
  strcat(ai.c, temp);
  strcat(ai.c, "# NOTES:\n");  
  strcat(ai.c, "# \t-One lower bin flux is set to zero.\n");  
  strcat(ai.c, "# \t because of that, min and max of the\n");  
  strcat(ai.c, "# \t whole histogram are slightly shifted.\n");  
  strcat(ai.c, "# \t-There is one extra row (last) with zero\n");  
  strcat(ai.c, "# \t values. This is to plot with pgfplots.\n#\n");  
  strcat(ai.c, "# Columns:\n");  
  strcat(ai.c, "# 0: Left (lower) bin value\n");
  strcat(ai.c, "# 1: Number of pixels in no-noised image\n");
  strcat(ai.c, "# 2: Number of pixels in noised image.\n");

  writeasciitable ("mockhist.txt", &ai, int_cols, 
		   accu_cols, space, prec);  
  
  free(allhist);
  free(dallhist);
  free(noisedhist); 
  free(nonoisehist); 
}
