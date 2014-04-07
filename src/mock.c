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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h> 

#include <gsl/gsl_sort.h>          /* gsl_sort2_float */
#include <gsl/gsl_integration.h>   /* gsl_integration_qng*/
#include <gsl/gsl_rng.h>           /* used in setrandoms*/
#include <gsl/gsl_randist.h>       /* To make noise.*/
#include <sys/time.h>              /* generate random seed*/

#include "mock.h"
#include "stats.h"
#include "attaavv.h"
#include "raddist.h"
#include "freqdomain.h"
#include "arraymanip.h"
#include "fitsarrayvv.h"
#include "macrofunctions.h"

/* Integration relative error, look in GSL manual: */
#define EPSREL_FOR_INTEG 2 





/****************************************************************
 *****************         Profiles:         ********************
 ****************************************************************/
/* The Gaussian function at a point: junk is only used here to make it
   conform to the general format of functions in this file: like
   Sersic(). */
double
Gaussian(double r, double junk, double a)
{
  junk=1;
  return exp( a*r*r );
}





/* This function will find the moffat function alpha value based on
   the explantions here:
 
   http://labs.adsabs.harvard.edu/adsabs/abs/2001MNRAS.328..977T/

   When alpha=(fwhm/2)/(2^(1.0/beta)-1)^(0.5). Then the moffat
   function at r is: (1.0 + (r/alpha)^2.0)^(-1.0*beta)*/
double
moffat_alpha(double fwhm, double beta)
{
    return (fwhm/2)/pow((pow(2, 1/beta)-1), 0.5f);
}





/* Find the Moffat profile for a certain radius. 
 
   rda=r/alpha     and nb=-1*b.

   This is done before hand to speed up the process. Junk is only used
   here to make it conform to the general format of functions in this
   file: like Sersic(). */
double
Moffat(double rda, double nb, double junk)
{
  junk=1;
  return pow(junk+rda*rda, nb);
}





/* This approximation of b(n) for n>0.35 is taken from McArthur,
   Courteau and Holtzman 2003:
   http://adsabs.harvard.edu/abs/2003ApJ...582..689 */
double
sersic_b(double n)
{
  if(n<=0.35f) 
    {
      printf("\n\n   ##ERROR in mock.c's sersic_b()\n");
      printf("\t Sersic n (=%f) must be smaller than 0.35\n\n", n);
      exit(EXIT_FAILURE);
    }
  return 2*n-(1/3)+(4/(405*n))+(46/(25515*n*n))+
    (131/(1148175*n*n*n)-(2194697/(30690717750*n*n*n*n)));
}





/* Find the Sersic profile for a certain radius. rdre=r/re, inv_n=1/n,
   nb= -1*b.  */
double
Sersic(double rdre, double inv_n, double nb)
{
    return exp( nb*( pow(rdre,inv_n)-1 ) );
}




















/****************************************************************
 *****************      2D integration:      ********************
 ****************************************************************/
/* Calculate a profile on position (i,j) in a cartesian
   grid. While the ellipse of the profile is rotated and has axis
   ratio q. i,j: positions on the grid. */
double
rot_ell(double x, struct mockparams *p)
{
  double c=p->c, s=p->s;
  double r;
  r=sqrt((x*c+p->y*s)*(x*c+p->y*s)+
	 ((p->y*c-x*s)*(p->y*c-x*s)/p->q/p->q));
  return p->profile(r/p->p2, p->p1, p->co);
}





/* This function is used in the integration of a profile. It
   assumes a fixed y and integrates over a range of x values.  */
double
twod_over_x(double x, void *params)
{
  struct mockparams *p;
  p=(struct mockparams *)params;

  return rot_ell(x, params);
}





/* Find the 2d integration over the region. */
double
twod_over_xy(double y, void *params)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  struct mockparams *p;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result;

  F.function = &twod_over_x;
  F.params = params;

  p=(struct mockparams *)params;
  p->y=y;
  gsl_integration_qng(&F, p->xl, p->xh, epsabs, epsrel, 
		      &result, &abserr, &neval);
  return result;
}




/* 2D integration of a profile.*/
double
integ2d(struct mockparams *params)
{
    gsl_function F;
    static double abserr;
    static size_t neval=0;
    double epsabs=0, epsrel=EPSREL_FOR_INTEG, result; 

    F.function = &twod_over_xy;
    F.params = params;   
    gsl_integration_qng(&F, params->yl, params->yh, epsabs, 
            epsrel, &result, &abserr, &neval);
    return result;
}




















/****************************************************************
 *****************       One profile:        ********************
 ****************************************************************/
void
fillmock(float **mock, float x_c, float y_c, struct mockparams *p, 
        float integaccu, float trunc_r, size_t *x_w, size_t *y_w, 
        double *totalflux, size_t *numpixs)
{
  float *pmock;
  size_t i, size, stride=1;
  float *actual_index;  /* Used to sort the mock image twice.        */
  double integ;         /* To temporarily keep the integrated value. */
  double point;         /* To compare the integrated value with.     */
  float t_i, t_j;       /* To read the position from the 1D index.   */
  float m_i_c, m_j_c;   /* The profile's center relative to its grid.*/  
  double p1=p->p1, p2=p->p2, co=p->co;

  makecanvas(trunc_r, p->q*trunc_r, p->pa_r, x_c, 
	     y_c, &m_i_c, &m_j_c, x_w, y_w, mock);

  *numpixs=0;
  *totalflux=0;
  pmock=*mock; 

  size=(*x_w)*(*y_w);
  actual_index=malloc(size*sizeof(float));
  assert(actual_index!=NULL);
  for(i=0;i<size;i++) actual_index[i]=i;
  gsl_sort2_float(pmock, stride, 
		  actual_index, stride, size);

  /* For the points that need integration: */
  for(i=0;i<size;i++)
    {
      t_j=(size_t)actual_index[i]%(*y_w);
      t_i=(actual_index[i]-t_j)/(*y_w);
      t_i=t_i-m_i_c; t_j=t_j-m_j_c;
      p->xl=t_i-0.5; p->xh=t_i+0.5;
      p->yl=t_j-0.5; p->yh=t_j+0.5;
      integ=integ2d(p);
      point=p->profile(pmock[i]/p2, p1, co);
      *totalflux+=pmock[i]=integ;
      (*numpixs)++;
      if (fabs(integ-point)/integ<integaccu) 
        {
	  i++;
	  break;
        }
    }
 
  /* For the points inner to the truncation radius
     where the elliptical radius is already calculated. */
  for(;i<size;i++)
    {
      if(pmock[i]>trunc_r) break;
      else 
        { 
          *totalflux+=pmock[i]=p->profile(pmock[i]/p2, p1, co);
          (*numpixs)++;
        }
    }


  /* For those points out of the truncation radius: */
  for(;i<size;i++) pmock[i]=0;

  gsl_sort2_float(actual_index, stride, pmock, stride, size);

  free(actual_index);
}





/* Make a profile based on the input parameters. */
#define SHOWONEPROFILE 0
void
oneprofile(float x_c, float y_c, float p1, float p2, float pa_d, 
        float q, float trunc, float integaccu, int av0_tot1, 
        float flux, float s0_m1_g2, float **mock, size_t *x_w, 
        size_t *y_w, size_t *numpixs)
{
  float trunc_r;
  char profletter;
  double totalflux=0;
  struct mockparams p;
  int is0_m1_g2=s0_m1_g2;
  static size_t fitscounter=1;
  char fitsname[100];

  switch(is0_m1_g2)
    {
    case 0: /* Sersic: p1: n, p2: re */
      trunc_r=p2*trunc;
      p.p1=1/p1; p.p2=p2;
      p.co=-1*sersic_b(p1);
      p.profile=&Sersic;
      profletter='s';
      break;
    case 1: /* Moffat: p1: beta. p2: fwhm, later alpha*/
      trunc_r=(p2/2)*trunc;
      p.p1=-1*p1; p.co=0;
      p.p2=moffat_alpha(p2, p1);
      p.profile=&Moffat;
      profletter='c';
      break;
    case 2: /* Gaussian: p1: sigma*/
      trunc_r=p1*trunc;
      p.p1=0; p.p2=1;
      p.co=-0.5f/p1/p1;
      p.profile=&Gaussian;
      profletter='g';
      break;
    default:
      printf("\n\n\tError: s0_m1_g2 must be 0, 1 or 2\n");
      printf("\t\tIt is: %d\n\n", is0_m1_g2);
      exit(EXIT_FAILURE);
      break;
    }
  p.pa_r=pa_d*3.141592f/180;
  p.c=cos(p.pa_r);
  p.s=sin(p.pa_r);
  p.q=q;

  fillmock(mock, x_c, y_c, &p, integaccu, trunc_r, 
	   x_w, y_w, &totalflux, numpixs);

  /* if av0_tot1==0, flux corresponds to the average flux,
     so the total flux will be flux*numpixs.  */
  if(av0_tot1==0)
    flux *= *numpixs;   
  /* Set the total flux of the final image. */
  totalflux=floatsum(*mock, *x_w * *y_w);
  floatarrmwith(*mock, *x_w * *y_w, flux/totalflux);  

  if(SHOWONEPROFILE)
    {
      sprintf(fitsname, "%c%lu.fits", profletter, fitscounter++);
      array_to_fits(fitsname, NULL, "MOCK", "FLOAT", 
		    *mock, *x_w, *y_w);
    }
}




















/****************************************************************
 *****************      Random variables     ********************
 ****************************************************************/
/* Simple function to generate a seed for a random number generator.*/
unsigned long int 
random_seed()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return(tv.tv_sec + tv.tv_usec);
}





/* Fill the array, with random numbers from a uniform
   distribution. Since the array might be two dimentional, the "array"
   pointer need not point to the first element.  It just needs to
   point to the first element in a column.  In case you want to fill a
   whole column of a 2D array, you point "array" to the first element
   in that column of the array, "num_rand" is the number of rows and
   "stride" is the number of clumns.

   The output of gsl_rng_uniform gives results in the range [0,1). So
   by multiplying the result by (range=higher_limit-lower_limit) we
   scale it to[0,range). Finally by adding the lower limit to the
   random variable, we ensure that the result is in the range:
   [lower_limit,higher_limit).*/
void
setunifrandvalues(double *array, size_t num_rand, size_t stride, 
        double lower_limit, double higher_limit)
{
  size_t i;
  gsl_rng * r;
  double range;
  const gsl_rng_type * T;

  range=higher_limit-lower_limit;

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,random_seed());

  for (i=0; i<num_rand; i++)
      array[i*stride] =lower_limit + gsl_rng_uniform(r)*range;

  gsl_rng_free (r);
}





/* An already allocated array is input into this function and is
   filled with random Gaussians, with mode of zero and standard
   deviation of sigma. */
void addnoise(float *array, size_t size, double sky)
{
  size_t i;
  gsl_rng * r;
  const gsl_rng_type * T;

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);   
  gsl_rng_set(r,random_seed());

  for(i=0;i<size;i++)
    array[i]+=gsl_ran_gaussian(r,sqrt(sky+array[i]));

  gsl_rng_free (r);
}




















/****************************************************************
 *****************    Set random parameters  ********************
 ****************************************************************/
/* Set the profile parameters needed for the mock image. */
void
setprflprms(double **prflprms, size_t numprflprms, size_t *nummock,
	    int size1, int size2)
{
  size_t i;  
  double *pprflprms;
  float randparamranges[14]={-10,+10,-10,+10,0.5,8,1,30,
			     0,360,0.2,1,0.005,0.1};
  randparamranges[1]+=size1;
  randparamranges[3]+=size2;
    
  *prflprms=malloc(numprflprms * *nummock*sizeof(double));
  assert(*prflprms!=NULL);
  pprflprms=*prflprms;
    
  /* We will assume all desired profiles are Sersic: */
  for(i=0;i<*nummock;i++)
    {
      pprflprms[i*numprflprms]=i+1;
      pprflprms[i*numprflprms+1]=0;
      pprflprms[i*numprflprms+numprflprms-1]=0;
    }

  /* Note that the last column is set for the total flux
     which is filled after the profile has been made. */
  for(i=2;i<numprflprms-1;i++)
    setunifrandvalues(&pprflprms[i], *nummock, numprflprms, 
		      randparamranges[(i-2)*2], 
                      randparamranges[(i-2)*2+1]);
}





/* This function will make a mock profile and place it in the larger image.*/
void
profileinimg(float *img, size_t size1, size_t size2, float trunc, 
        float integaccu, float s0_m1_g2, float x_c, float y_c, 
        float p1, float p2, float pa_d, float q, float avflux, 
        double *totflux)
{
  float *mock;
  int av0_tot1=0;
  int i, j, ilimit, jlimit;
  size_t x_w, y_w, xc, yc, numpixs;
  int x1l, y1l, x2l, y2l, x1s, y1s, x2s, y2s;

  if(x_c-(int)x_c>=0.5) xc=x_c+1;
  else xc=x_c;
  if(y_c-(int)y_c>=0.5) yc=y_c+1;
  else yc=y_c;

  oneprofile(x_c, y_c, p1, p2, pa_d, q, trunc, integaccu, 
	     av0_tot1, avflux, s0_m1_g2, &mock, &x_w, &y_w, 
	     &numpixs);
  *totflux=numpixs*avflux;

  smallinlarge((int)size1, (int)size2, (int)x_w, (int)y_w, 
	       (int)xc, (int)yc, &x1l, &y1l, &x2l, &y2l, &x1s, 
	       &y1s, &x2s, &y2s);

  ilimit=x2s-x1s; jlimit=y2s-y1s;
  for(i=0;i<ilimit;i++)
    for(j=0;j<jlimit;j++)
      img[(x1l+i)*size2+y1l+j]+=mock[(x1s+i)*y_w+(y1s+j)];

  free(mock);    
}




















/****************************************************************
 *****************  Save mock image and info ********************
 ****************************************************************/
/* If the mock image is to be saved, save the information
   of the galaxies and the actual mock image. */
void
saveimgandinfo(double *prflprms, size_t nummock, size_t numprflprms, 
	       float sky, float trunc, size_t hsize1, size_t hsize2, 
	       char *infoname)
{
  size_t i;
  float ftemp;
  char temp[1000];
  struct ArrayInfo ai;
  int int_cols[]={0, 1, 6,-1}, accu_cols[]={2,3,9,-1};
  int space[]={6,8,11}, prec[]={2,4};

  ai.c=malloc(MAXALLCOMMENTSLENGTH*sizeof(char));
  assert(ai.c!=NULL);
  ai.s0=nummock;
  ai.s1=numprflprms;
  ai.d=prflprms;
  /* Change to FITS standard. */
  for(i=0;i<nummock;i++){/* Remove extra spacing. */  
    ftemp=prflprms[i*numprflprms+2]-hsize1+1;
    prflprms[i*numprflprms+2]=prflprms[i*numprflprms+3]-hsize2+1;
    prflprms[i*numprflprms+3]=ftemp;
  }

  sprintf(temp, "# Properties of %lu mock galaxies.\n", nummock);
  strcpy(ai.c, temp);
  sprintf(temp, "# The sky valued is assumed to be: %.2f\n", sky);
  strcat(ai.c, temp);
  sprintf(temp, "# Truncation at %.2f * radial parameter\n# \n", trunc);
  strcat(ai.c, temp);
  strcat(ai.c, "# 0: ID.\n");
  strcat(ai.c, "# 1: 0: Sersic, 1: Moffat, 2: Gaussian.\n");
  strcat(ai.c, "# 2: X position (FITS definition).\n");
  strcat(ai.c, "# 3: Y position (FITS definition).\n");
  strcat(ai.c, "# 4: Sersic n or Moffat beta.\n");    
  strcat(ai.c, "# 5: Sersic re or Moffat FWHM.\n");    
  strcat(ai.c, "# 6: Position angle, degrees.\n");    
  strcat(ai.c, "# 7: Axis ratio.\n");    
  strcat(ai.c, "# 8: Signal to noise: ");
  strcat(ai.c, "(average profile flux-sky)/sqrt(sky).\n");    
  strcat(ai.c, "# 9: Total flux (Sky subtracted).\n\n"); 

  writeasciitable (infoname, &ai, int_cols, accu_cols, space, prec);
}




















/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************     Main output program   ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/* Put one or more mock profiles into and image, convolve it and add
   noise to the final result.  The convolution is going to make the
   sides darker.  So the actual image where the galaxies will be
   placed is going to be larger than the desired image.  After
   convolution those sides will be trimed and the pixels on the sides
   of the result will be equally convolved as the central pixels.

   Note on MINFLOAT: We are using float values here, due to roundoff
   errors, after convolution we might have pixel values less than
   NUMMOCK, which are pure error, so they will all be set to zero. */
#define REPORTMOCK 1
#define MINFLOAT 1e-6f
void
mockimg(size_t size1, size_t size2, float sky, size_t nummock, 
        double *prflprms, float psf_fwhm, float psf_beta, int vpsf, 
        float **mock, char *outname, char *infoname)
{
  int free_1_0=0;
  size_t hsize1, hsize2;
  float omin=-250, omax=700;
  int psf_smg=1, psf_av0_tot1=1;
  float integaccu=0.01, psfsum=1;
  size_t i, psf_size1, psf_size2, junk;
  size_t nsize1, nsize2, numprflprms=10;
  float *psf, psf_q=1, psf_pa=0, trunc=5;
  float *nonoisehist, *noisedhist, *allhist;

  oneprofile(0.0f, 0.0f, psf_beta, psf_fwhm, psf_pa, psf_q, 
	     trunc, integaccu, psf_av0_tot1, psfsum, psf_smg, 
	     &psf, &psf_size1, &psf_size2, &junk);

  if(vpsf)
    array_to_fits("PSF.fits", NULL, "PSF", "FLOAT", 
		  psf, psf_size1, psf_size2);

  hsize1=psf_size1/2;       hsize2=psf_size2/2;
  nsize1=size1+2*hsize1;    nsize2=size2+2*hsize2;

  if(prflprms==NULL)
    {
      free_1_0=1;
      setprflprms(&prflprms, numprflprms, &nummock, 
		  nsize1, nsize2);
    }

  *mock=calloc(nsize1*nsize2,sizeof(float));
  assert(*mock!=NULL);
  
  for(i=0;i<nummock;i++)
    profileinimg(*mock, nsize1, nsize2, trunc, integaccu, 
		 prflprms[i*numprflprms+1], prflprms[i*numprflprms+2], 
		 prflprms[i*numprflprms+3], prflprms[i*numprflprms+4], 
		 prflprms[i*numprflprms+5], prflprms[i*numprflprms+6], 
		 prflprms[i*numprflprms+7], 
		 sqrt(sky)*prflprms[i*numprflprms+8], 
		 &prflprms[i*numprflprms+9]);

  convolve(*mock, nsize1, nsize2, psf, psf_size1, psf_size2);

  floatshrinkarray(mock, nsize1, nsize2, hsize1, hsize2, 
		   size1+hsize1, size2+hsize2);

  floatsetbelowtozero(*mock, size1*size2, MINFLOAT);

  array_to_fits(outname, NULL, "NONOISE", "FLOAT", 
		*mock, size1, size2);
   

  if(MOCKHIST)
    histogram(*mock, size1*size2, MOCKHISTNUMBINS, 
	      &omin, &omax, &nonoisehist, 1, 0, 0);

  addnoise(*mock, size1*size2, sky);

  array_to_fits(outname, NULL, "WITHNOISE", "FLOAT", 
		*mock, size1, size2);
   

  if(MOCKHIST)
    {
      histogram(*mock, size1*size2, MOCKHISTNUMBINS, 
		&omin, &omax, &noisedhist, 1, 0, 0);
      floatvmerge(nonoisehist, noisedhist, MOCKHISTNUMBINS, &allhist);
      printfarray(allhist, MOCKHISTNUMBINS, 3, "", "mockhists.txt", 6);
      free(nonoisehist); free(noisedhist); free(allhist);
    }

  if(REPORTMOCK)
    if(free_1_0==1)
      saveimgandinfo(prflprms, nummock, numprflprms, sky, 
		     trunc, hsize1, hsize2, infoname);
   
  if(free_1_0==1) free(prflprms);    
  free(psf);
}
