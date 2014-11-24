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
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include "mock.h"
#include "stats.h"
#include "config.h"
#include "attaavv.h"
#include "arraymanip.h"
#include "argpparser.h"
#include "fitsarrayvv.h"

#include "ui.h"			/* Needs mock.h */





/****************************************************************
 *****************     Check input values     *******************
 ****************************************************************/
void
anyfloat(char *optarg, float *var, char *lo, char so)
{
  float tmp;
  char *tailptr;
  tmp=strtof(optarg, &tailptr);
  if(strlen(tailptr))
    {
      fprintf(stderr, PACKAGE": the argument to option `-%c`: `%s` was not "
	     "readable as a number!\n", so, optarg);
      exit(EXIT_FAILURE);
    }
  *var=tmp;
}





void
floatl0(char *optarg, float *var, char *lo, char so)
{
  float tmp;
  char *tailptr;
  tmp=strtof(optarg, &tailptr);
  if(strlen(tailptr))
    {
      fprintf(stderr, PACKAGE": the argument to option `-%c`: `%s` was not "
	     "readable as a number!\n", so, optarg);
      exit(EXIT_FAILURE);
    }
  if(tmp<=0)
    {
      fprintf(stderr, PACKAGE": `--%s (-%c)` > 0, but it is: %.3f\n", 
	     lo, so, tmp);
      exit(EXIT_FAILURE);      
    }
  *var=tmp;
}





void
floatel0(char *optarg, float *var, char *lo, char so)
{
  float tmp;
  char *tailptr;
  tmp=strtof(optarg, &tailptr);
  if(strlen(tailptr))
    {
      fprintf(stderr, PACKAGE": the argument to option `-%c`: `%s` "
	     "was not readable as a number!\n", so, optarg);
      exit(EXIT_FAILURE);
    }
  if(tmp<0)
    {
      fprintf(stderr, PACKAGE": `--%s (-%c)` => 0, but it is: %.3f\n", 
	     lo, so, tmp);
      exit(EXIT_FAILURE);      
    }
  *var=tmp;
}





void
sizetlzero(char *optarg, size_t *var, char *lo, char so)
{
  long tmp;
  char *tailptr;
  tmp=strtol(optarg, &tailptr, 0);
  if(strlen(tailptr))
    {
      fprintf(stderr, PACKAGE": the argument to option `-%c`: `%s` was not "
	     "readable as a number!\n", so, optarg);
      exit(EXIT_FAILURE);
    }
  if(tmp<=0)
    {
      fprintf(stderr, PACKAGE": argument to `--%s (-%c)` should be >0, it "
	     "is: %ld\n", lo, so, tmp);
      exit(EXIT_FAILURE);
    }
  *var=tmp;
}




















/****************************************************************
 *****************    Output names and file    ******************
 ****************************************************************/
void
changenameedning(char *in, char *append, char **outname, int removenamedir)
{
  char *out;
  size_t i, l, al, offset=0;

  l=strlen(in);
  al=strlen(append);
  assert( (*outname=out=malloc((l+al+5)*sizeof *out))!=NULL );
  strcpy(out, in);
  for(i=l;i!=0;--i)
    if(out[i]=='.')
      {
	out[i]='\0';
	strcat(out, append);
	break;
      }
  if(i==0)
    {
      fprintf(stderr, "%s: %s does not have a '.' in the name\n", 
	     PACKAGE, in);
      exit(EXIT_FAILURE);
    }

  /* If it is desired to remove the directory information from the
     name, do it here: */
  if(removenamedir)
    {
      l=strlen(out);
      for(i=l;i!=0;--i)	 	  /* Find the last forward slash.      */
	if(out[i]=='/')
	  {offset=i+1; break;}
      if(offset)
	for(i=offset;i<=l;++i)	  /* <= because we want to shift the   */
	  out[i-offset]=out[i]; /* '\0' character in the string too. */
    }
}



















/****************************************************************
 *****************     Input/Output reading    ******************
 ****************************************************************/
/* If the input table exists, read it and put its pointer into
   mockparams. If it doesn't exist, then let prflprms(), which is
   is mock.c, make a random set of parameters. */
void
readinputinfo(struct mockparams *p)
{
  double *pp;
  size_t i, j;
  FILE *tmpfile;
  struct ArrayInfo intable;

  /* Set the number of columns in the output catalog: */
  p->numppcols=10;

  if ((tmpfile = fopen(p->up.infoname, "r")) != NULL) 
    {
      fclose(tmpfile);
      readasciitable(p->up.infoname, &intable);
      if(intable.s1!=p->numppcols-1)
	{
	  fprintf(stderr, PACKAGE": %s not %lu columns.", 
		  p->up.infoname, p->numppcols-1);
	  exit(EXIT_FAILURE);
	}
      p->nummock=intable.s0;

      /* Put the information of the input table into an output table: */
      assert( (pp=malloc(p->nummock*p->numppcols*sizeof *pp))!=NULL );
      for(i=0;i<p->nummock;++i)
	for(j=0;j<p->numppcols-1;++j)
	  pp[i*p->numppcols+j]=intable.d[i*intable.s1+j];
      p->profileparams=pp;
      
      /* Free the input parameters: */
      freeasciitable(&intable);
    }
  else 
    {
      fprintf(stderr, PACKAGE": %s could not be read.\n", p->up.infoname);
      exit(EXIT_FAILURE);
    }
}





/* Read the point spread function FITS file: */
void
readpsf(struct mockparams *p)
{
  void *tmp;
  int bitpix;
  float psfsum;
  FILE *tmpfile;

  /* If a name is not specified, go on to build the PSF later. */
  if(p->up.psfname==NULL) return;

  /* The name is give, so check if the file exists. */
  if ((tmpfile = fopen(p->up.psfname, "r")) != NULL) 
    {			/* The file exists! */
      fclose(tmpfile);
      fits_to_array(p->up.psfname, 0, &bitpix, &tmp,
		    &p->psf_s0, &p->psf_s1);

      /* Check if the sides are odd */
      if(p->psf_s0%2==0)
	{
	  fprintf(stderr, PACKAGE": NAXIS2 of %s (PSF) must be odd! But "
		  "it is %lu.\n", p->up.psfname, p->psf_s0);
	  exit(EXIT_FAILURE);
	}
      if(p->psf_s1%2==0)
	{
	  fprintf(stderr, PACKAGE": NAXIS1 of %s (PSF) must be odd! But "
		  "it is %lu.\n", p->up.psfname, p->psf_s1);
	  exit(EXIT_FAILURE);
	}

      /* Everything is fine, make sure the output is in float. */
      convertanytofloat(tmp, p->psf_s0 * p->psf_s1, bitpix, &p->psf,
			BYTE_IMG, SHORT_IMG, LONG_IMG, FLOAT_IMG, 
			DOUBLE_IMG);

      /* Check if the sum of the PSF is unity. */
      psfsum=floatsum(p->psf, p->psf_s0 * p->psf_s1);
      if(psfsum!=1)
	floatarrmwith(p->psf, p->psf_s0 * p->psf_s1, 1/psfsum);
    }
  else
    {
      fprintf(stderr, PACKAGE": %s could not be opened.\n",
	      p->up.psfname);
      exit(EXIT_FAILURE);
    }
}





void
changenameending(char *in, char *append, char **outname, int removenamedir)
{
  char *out;
  size_t i, l, al, offset=0;

  l=strlen(in);
  al=strlen(append);
  assert( (*outname=out=malloc((l+al+5)*sizeof *out))!=NULL );
  strcpy(out, in);
  for(i=l;i!=0;--i)
    if(out[i]=='.')
      {
	out[i]='\0';
	strcat(out, append);
	break;
      }

  /* If it is desired to remove the directory information from the
     name, do it here: */
  if(removenamedir)
    {
      l=strlen(out);
      for(i=l;i!=0;--i)	 	  /* Find the last forward slash.      */
	if(out[i]=='/')
	  {offset=i+1; break;}
      if(offset)
	for(i=offset;i<=l;++i)	  /* <= because we want to shift the   */
	  out[i-offset]=out[i]; /* '\0' character in the string too. */
    }
}





/* Check if a file exists. If so, remove it. */
void
checkremovefile(char *filename, int dontdelete)
{
  FILE *tmpfile;
  if ((tmpfile = fopen(filename, "r")) != NULL) 
    {
      fclose(tmpfile);
      if(dontdelete)
	{
	  fprintf(stderr, PACKAGE": '%s' already exists and you have "
		  "asked to not remove it with the `--dontdelete` "
		  "(`-D`) option.\n", filename);
	  exit(EXIT_FAILURE);
	}
      else
	if(unlink(filename)==-1)
	  {
	    fprintf(stderr, PACKAGE": '%s' already exists and "
		    "could not be removed\n", filename); 
	    exit(EXIT_FAILURE);
	  }
    }
}





void
makecheckoutputnames(struct mockparams *p)
{
  /* If no output name is specified, then use the input name. */
  if(p->up.outname==NULL)
    p->up.outname=p->up.infoname;

  /* Set the name endings: */
  changenameending(p->up.outname, ".fits", &p->fitsname,
		   p->up.removenamedir);
  changenameending(p->up.outname, "_cat.txt", &p->catname,
		   p->up.removenamedir);
  checkremovefile(p->fitsname, p->up.dontdelete);
  checkremovefile(p->catname, p->up.dontdelete);

  if(p->outpsfname || p->onlypsf)
    {
      changenameending(p->up.outname, "_psf.fits", &p->outpsfname,
		       p->up.removenamedir);
      checkremovefile(p->outpsfname, p->up.dontdelete);
    }
  if(p->sepsfname)
    {
      changenameending(p->up.outname, "_psf.conv", &p->sepsfname,
		       p->up.removenamedir);
      checkremovefile(p->sepsfname, p->up.dontdelete);
    }
}















/****************************************************************
 *****************        Read options:      ********************
 ****************************************************************/
void
setparams(struct mockparams *p, int argc, char *argv[],
	  time_t *rawtime)
{
  /* Default values: */
  p->up.infoname      = NULL;
  p->up.outname       = NULL;
  p->up.psfname       = strlen(DP_PSF) ? DP_PSF : NULL;
  p->up.dontdelete    = 0;
  p->up.removenamedir = 1;
  
  p->truncation       = DP_TRUNCATION_V;
  p->tolerance        = DP_TOLERANCE_V;
  
  p->psf              = NULL;
  p->outpsfname       = NULL;
  p->sepsfname        = NULL;
  p->psffunction      = DP_PSFFUNCTION_V;
  p->psf_p1           = DP_PSFFWHM_V;
  p->psf_p2           = DP_MOFFATBETA_V;
  p->psf_t            = DP_PSFTRUNC_V;
  p->onlypsf          = 0;

  p->s0               = DP_NAXIS2_V; /* Recall that the C standard and   */
  p->s1               = DP_NAXIS1_V; /* FITS have different side lengths.*/
  p->background       = DP_BACKGROUND_V;
  p->zeropoint        = DP_ZEROPOINT_V;
  p->viewnoconv       = 0;
  p->viewconv         = 0; 

  p->verb             = 1;
  p->numthreads       = 1;

  /* Read the arguments: */
  argp_parse(&argp, argc, argv, 0, 0, p);

  /* Check output names: */
  makecheckoutputnames(p);
  
  /* Check and read input catalog. */
  readinputinfo(p);
  
  /* Check and read PSF. */
  readpsf(p);

  /* There were no inconsistencies and input values are read, report
     the starting of MockGals. */
  if(p->verb)
    {
      printf("\n\n--------------------------\n");
      printf("mockgals started on %s\n", ctime(rawtime));
      printf(" - Information of %lu profile%sread from '%s'.\n",
	     p->nummock, p->nummock>1 ? "s " : " ", p->up.infoname);
    }
}




















/****************************************************************
 **********  Free allocated space and report timing    **********
 ****************************************************************/
void
freeandreporttime(struct mockparams *p, struct timeval *t0)
{
  struct timeval t1;
  
  free(p->catname);
  free(p->fitsname);
  free(p->profileparams);

  gettimeofday(&t1, NULL);
  printf("mockgals finished in %.4f (seconds)\n",
	 ((double)t1.tv_sec+(double)t1.tv_usec/1e6) - 
	 ((double)t0->tv_sec+(double)t0->tv_usec/1e6));
  if(p->verb)
    {
      printf("--------------------------\n\n");
    }
}
