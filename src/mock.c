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
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <sys/time.h>		 /* generate random seed*/
#include <gsl/gsl_rng.h>	 /* used in setrandoms*/
#include <gsl/gsl_sort.h>	 /* gsl_sort2_float */
#include <gsl/gsl_randist.h>	 /* To make noise.*/
#include <gsl/gsl_sf_gamma.h>	 /* For total Sersic flux. */
#include <gsl/gsl_integration.h> /* gsl_integration_qng*/

#include "pix.h"
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

   alpha=(FWHM/2)/(2^(1.0/beta)-1)^(0.5). Then the moffat
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
  junk=1; /* So it can keep the same functional form as the others. */
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




/* Find the total flux in a Sersic profile. From equation 4 in Peng
   2010. This assumes the surface brightness at the effective radius
   is 1.*/
double
totsersic(double n, double re, double b, double q)
{
  return (2*M_PI*re*re*exp(b)*n*pow(b, -2*n)*q*
	  gsl_sf_gamma(2*n));
}





/* For a point source. */
double
Point(double j1, double j2, double j3)
{
  j1=j2=j3;
  return 1;
}




















/****************************************************************
 *****************      2D integration:      ********************
 ****************************************************************/
/* Calculate a profile on position (i,j) in a cartesian
   grid. While the ellipse of the profile is rotated and has axis
   ratio q. i,j: positions on the grid. */
double
rot_ell(double x, struct integparams *p)
{
  double c=p->c, s=p->s;
  double r;
  r=sqrt((x*c+p->y*s)*(x*c+p->y*s)+
	 ((p->y*c-x*s)*(p->y*c-x*s)/p->q/p->q));
  return p->profile(r/p->p1, p->p2, p->co);
}





/* This function is used in the integration of a profile. It
   assumes a fixed y and integrates over a range of x values.  */
double
twod_over_x(double x, void *params)
{
  struct integparams *p;
  p=(struct integparams *)params;

  return rot_ell(x, params);
}





/* Find the 2d integration over the region. */
double
twod_over_xy(double y, void *params)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  struct integparams *p;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result;

  F.function = &twod_over_x;
  F.params = params;

  p=(struct integparams *)params;
  p->y=y;
  gsl_integration_qng(&F, p->xl, p->xh, epsabs, epsrel, 
		      &result, &abserr, &neval);
  return result;
}




/* 2D integration of a profile.*/
double
integ2d(struct integparams *params)
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





/* Set the general parameters of struct integparams. In all profiles,
   p1 is the parameter that the radius is divided by. */
void
setintegparams(int s0_m1_g2_p3, float p1, float p2, float pa_d, 
	       float q, float trunc, float *trunc_r, char *profletter, 
	       struct integparams *p)
{
  switch(s0_m1_g2_p3)
    {
    case 0: /* Sersic: p1: re, p2: n */
      *trunc_r=p1*trunc;	/* Trunctation in units of re. */
      p->p1=p1; 		/* re. */
      p->p2=1/p2;		/* n is used as 1/n in Sersic. */
      p->co=-1*sersic_b(p2);	/* Constant in the power. */
      p->profile=&Sersic;
      *profletter='s';
      break;
    case 1: /* Moffat: p1: fwhm(input) -> alpha. p2: beta*/
      *trunc_r=(p1/2)*trunc;	/* Truncation in units of FWHM/2. */
      p->p1=moffat_alpha(p1, p2);
      p->p2=-1*p2;              /* Beta is always negative!  */
      p->co=0;	                /* No constant terms needed */
      p->profile=&Moffat;
      *profletter='m';
      break;
    case 2: /* Gaussian: p1: FWHM(input) -> sigma */
      p1/=2.35482;              /* Convert FWHM to sigma. */
      *trunc_r=p1*trunc;	/* Truncation in units of sigma */
      p->co=-1.0f/(2.0f*p1*p1);	/* Constant to multiply */
      p->p1=1;			/* r gets divided by p->p1! */
      p->p2=0;			/* Not needed! */
      p->profile=&Gaussian;   
      *profletter='g';
      break;
    case 3: /* Point:*/
      *trunc_r=2;
      p->p1=p->p2=0;
      p->co=1;
      p->profile=&Point;
      *profletter='p';
      break;
    default:
      printf("\n\n\tError: s0_m1_g2 must be 0, 1 or 2\n");
      printf("\t\tIt is: %d\n\n", s0_m1_g2_p3);
      exit(EXIT_FAILURE);
      break;
    }
  p->pa_r=pa_d*M_PI/180;
  p->c=cos(p->pa_r);
  p->s=sin(p->pa_r);
  p->q=q;
}


















/****************************************************************
 *****************       One profile:        ********************
 ****************************************************************/
/* Sometimes (for example in making a PSF) you just want one profile
   that is fully enclosed in an array. The two functions here, namely
   oneprofile() and fillmock(). Were made for this purpose. These were
   the precursors to the currently used makeprofile(), they were
   called by makeprofile_old(), and the intersection of the profile
   made here and the image (depending on where the actual image was
   asked to be) would be added to the main output image. 

   makeprofile_old() was removed from this source code on April 12th,
   2014. You can see it in the last commit before that.*/
void
fillmock(float **mock, float x_c, float y_c, struct integparams *p, 
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
        float flux, float s0_m1_g2_p3, float **mock, size_t *x_w, 
        size_t *y_w, size_t *numpixs)
{
  float trunc_r;
  char profletter;
  double totalflux=0;
  struct integparams p;
  int is0_m1_g2_p3=s0_m1_g2_p3;
  static size_t fitscounter=1;
  char fitsname[100];

  setintegparams(is0_m1_g2_p3, p1, p2, pa_d, q, trunc, &trunc_r, 
		 &profletter, &p);

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
      array_to_fits(fitsname, NULL, "ONEPROFILE", FLOAT_IMG, 
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





/* Set the profile parameters needed for the mock image. */
void
setprflprms(double **prflprms, size_t numprflprms, 
	    size_t nummock, int size1, int size2)
{
  size_t i;  
  double *pprflprms;
  double starSN[5]={400, 150, 1000, 300, 110};
  float randparamranges[14]={-20,+20,     /* x axis range. */
			     -20,+20,     /* y axis range. */
			     8,12,        /* re range. */
			     0.5,8,       /* n range */
			     0,360,       /* position angle range */
			     0.2,1,       /* Axis ration range. */
			     0.0005,0.001}; /* S/N range. */
  randparamranges[1]+=size1;
  randparamranges[3]+=size2;
    
  *prflprms=malloc(numprflprms * nummock*sizeof(double));
  assert(*prflprms!=NULL);
  pprflprms=*prflprms;
    
  /* We will assume all desired profiles are Sersic: */
  for(i=0;i<nummock;i++)
    {
      pprflprms[i*numprflprms]=i+1;
      pprflprms[i*numprflprms+1]=0;
      pprflprms[i*numprflprms+numprflprms-1]=0;
    }

  /* Note that the last column is set for the total flux
     which is filled after the profile has been made. */
  for(i=2;i<numprflprms-1;i++)
    setunifrandvalues(&pprflprms[i], nummock, numprflprms, 
		      randparamranges[(i-2)*2], 
                      randparamranges[(i-2)*2+1]);

  /* Put 5 stars in the image: */
  for(i=0;i<5;i++)
    {
      pprflprms[i*numprflprms+1]=3;
      pprflprms[i*numprflprms+8]=starSN[i];
    }
}





/* An already allocated array is input into this function and is
   filled with random Gaussians, with mode of zero and standard
   deviation of sigma. */
void 
addnoise(float *array, size_t size, double sky)
{
  size_t i;
  gsl_rng * r;
  const gsl_rng_type * T;

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);   
  gsl_rng_set(r,random_seed());

  for(i=0;i<size;i++)
    array[i]+=gsl_ran_gaussian(r,sqrt(sky+array[i]));
  
  gsl_rng_free(r);
}




















/****************************************************************
 *****************  Fill in the mock values  ********************
 ****************************************************************/





/* Find the first pixel in the image to begin building the profile.

   The input sizes and positions are based on the FITS standard,
   But in main(), we reversed the sizes to fits the C standard and
   when calling this function, we reversed the positions to fit
   the C standard. So by the time they get here, the inputs are
   all in the C standard.*/
void
findstartingpixel(struct pix *D, size_t s0, size_t s1, float truncr, 
		  struct elraddistp *e, struct pix **p)
{
  float rmin, x_c, y_c, fs0, fs1;
  size_t x_w, y_w, xmin=-1, ymin=-1;
  int is0, is1, i, j, x, y, x1, y1, x2, y2;

  is0=s0; is1=s1;
  x_c=e->xc; y_c=e->yc;

  /* Find the central pixel, this will be needed if it is inside the
     image or outside it. */
  if(x_c-(int)x_c>0.5) x=x_c+1;
  else x=x_c;
  if(y_c-(int)y_c>0.5) y=y_c+1;
  else y=y_c;

  /* The central pixel is in the image, set the pixel and return. */
  fs0=s0;			/* Just to make things easier! */
  fs1=s1;
  if(x_c>=0 && x_c<fs0 && y_c>=0 && y_c<fs1)
    {
      *p=D+x*s1+y;
      return;
    }

  /* The center is out of the image. Use encloseellipse() from
   raddist.c to see if any of the pixels within the truncation
   radius fit into the image.*/
  encloseellipse(truncr, e->q*truncr, e->t, &x_w, &y_w);
  x1=x-x_w/2;
  y1=y-y_w/2;
  x2=x+x_w/2;
  y2=y+y_w/2;
  
  /* Check if any of the four corners of the box inclosing the
     profile are in the mock image. If they are not, x1=x2 or
     y1=y2*/
  checkifinarray(&x1, &y1, &x2, &y2, s0, s1);
  if(x1==x2 || y1==y2)	    
    {			     /* The profile's region is */
      *p=NULL;		     /* Completely out of the image. */
      return;		     /* Return NULL. */
    }
  else			     /* The profile and the image overlap */
    {			     /* Find the point on the side of the */
      rmin=1e10;	     /* image with the smallest radius. */
      if(x1==0)		     /* This is important, because the  */
	for(j=y1;j<y2;j++)   /* first check later will be  */
	  if(elraddist(e, 0, j)<rmin)
	    {		     /* integration and we want to be sure  */
              xmin=0; ymin=j; /* we start with the smallest radius. */
	      rmin=elraddist(e, 0, j); /* in the image. */
	    }
      if(x2==is0)		       
	for(j=y1;j<y2;j++)
	  if(elraddist(e, s0-1, j)<rmin)
	    {
              xmin=s0-1; ymin=j;
	      rmin=elraddist(e, s0-1, j);
	    }
      if(y1==0)		       
	for(i=x1;i<x2;i++)
	  if(elraddist(e, i, 0)<rmin)
	    {
              xmin=i; ymin=0;
	      rmin=elraddist(e, i, 0);
	    }
      if(y2==is1)		       
	for(i=x1;i<x2;i++)
	  if(elraddist(e, i, s1-1)<rmin)
	    {
              xmin=i; ymin=s1-1;
	      rmin=elraddist(e, i, s1-1);
	    }
      if(rmin<truncr && xmin!=(size_t)-1 && ymin!=(size_t)-1)
	*p=D+xmin*s1+ymin;
      else *p=NULL;
    }
}





/* Make a profile in an array.

   The logic: byt is an array the same size as the image, that 
   should be completely zero upon making each profile. If it is
   zero at first, it will be zero once this function is finished
   with it. It is used to mark which pixels have been checked in 
   the image. indexs will keep the index (1D) of those pixels that
   have been marked so after the job is finished, it can set them
   all back to zero and not bother with resetting the whole array!

   We begin with the pixel in the image that is closest to the desired
   center of the profile. A pixel queue (a FIFO) will keep all the
   neighbors of pixels in order to check them all. Since some profiles
   will fall on the sides of the image, there is no way we can
   calculate the whole flux by summing of the pixels and then setting
   them to the desired value. So we have to use integration, the
   profiles were all made with their constant set to 1. We integrate
   over the profile to infinity over a surface and consider that as
   the total flux. Note that when the truncation radius is small, this
   theoretical total flux will be much larger (>10%) than the actual
   total flux of what is actually put in the image if it all fits in.
   Since that total flux will be set to the desired value, then if the
   truncation radius is too small, the real total flux in the image is
   slightly lower than the desired value. This problem did not exist
   in the makeprofile_old() function which actually allocated an array
   for each profile, summed over it to find the total flux and then
   only used the intersection with the main image to put the profile's
   image into the main image. But that was too slow, especially if a
   large number of profiles were needed. It was removed on April 12th,
   2014. So if you are interested, you can see it in the commits before
   this date.

   For all pixels, first it is checked if the pixel is within the
   truncation radius. If it isn't, its neighbors will not be added to
   the queue, but it's byt value will be set to one and a value of
   zero will be put into its D[i]->v value.
 */
int
makeprofile(float *img, unsigned char *byt, size_t *bytind, 
	    struct pix *D, size_t s0, size_t s1, float trunc, 
	    float integaccu, float s0_m1_g2_p3, float x_c, float y_c, 
	    double p1, double p2, float pa_d, float q, float avflux, 
	    double *totflux)
{
  struct pix *p;
  float t_i, t_j;		/* 2D position from 1D. */
  char profletter;
  struct elraddistp e;
  size_t j, counter=0;
  struct integparams ip;
  struct pixlist *Q=NULL;
  double sum=0, area=0, co;
  double (*func)(double, double, double);
  float r, truncr, maxir, integ, tmp, multiple=0;
  
  e.q=q;
  e.xc=x_c;
  e.yc=y_c;
  e.t=M_PI*pa_d/180;
  e.cos=cos(e.t);
  e.sin=sin(e.t);

  setintegparams(s0_m1_g2_p3, p1, p2, pa_d, q, trunc, 
		 &truncr, &profletter, &ip);

  findstartingpixel(D, s0, s1, truncr, &e, &p);
  if(p==NULL)
    return 0;			/* Profile is completely out. */

  if(s0_m1_g2_p3==0)
    {
      sum=totsersic(p1, p2, sersic_b(p1), q);
      area=M_PI*truncr*truncr*q;
      *totflux=avflux*area;
      multiple=*totflux/sum;
    }
  else if(s0_m1_g2_p3==3)
    {
      img[p->i]=avflux;
      *totflux=avflux;
      return 1;			/* Successful. */
    }

  addtopixlist(&Q, p);

  maxir=truncr;
  co=ip.co;
  func=ip.profile;
  p1=ip.p1;
  p2=ip.p2;

  while(Q!=NULL)
    {
      popfrompixlist(&Q, &p);

      /* A pixel might be added to this list more than once.
         check if it has already been checked:*/
      if(byt[p->i]==1) continue;

      byt[p->i]=1;		/* Mark as checked. */
      bytind[counter++]=p->i;

      t_i=p->i/s1;
      t_j=p->i%s1;
      if( (r=elraddist(&e, t_i, t_j)) > truncr) 
	{
	  p->v=0;		/* Read above. */
	  continue;
	}
      
      /* Find the value for this pixel: */
      tmp=func(r/p1, p2, co); 
      if(r<maxir)		
	{
	  t_i-=x_c;             t_j-=y_c;
	  ip.xl=t_i-0.5;       ip.xh=t_i+0.5;
	  ip.yl=t_j-0.5;       ip.yh=t_j+0.5;
	  integ=integ2d(&ip);	  
	  if (fabs(integ-tmp)/integ>integaccu) 
	    tmp=integ;
	  else maxir=r;
	}
      img[p->i]+=tmp*multiple;

      for(j=0;p->ngb8[j].p!=NULL;j++)
	if(byt[p->ngb8[j].p->i]==0)
	  addtopixlist(&Q, p->ngb8[j].p);
    }

  /* Clean the byt image for the next profile. */
  for(j=0;j<counter;j++)
    byt[bytind[j]]=0;
  return 1;
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

  sprintf(temp, "# Properties of %lu mock galaxies.\n", p->nummock);
  strcpy(ai.c, temp);
  sprintf(temp, "# The sky valued is assumed to be: %.2f\n", p->sky);
  strcat(ai.c, temp);
  sprintf(temp, "# Truncation at %.2f * radial parameter\n# \n", 
	  p->trunc);
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

















/****************************************************************
 *****************     Main output program   ********************
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
#define MINFLOAT 1e-8f
void
mockimg(struct mockparams *p)
{
  struct pix *D;
  double *pp, ss;
  unsigned char *byt;
  int extcounter=0, suc;
  int psf_smg=1, psf_av0_tot1=1;
  float *img, integaccu=0.01, psfsum=1;
  float *nonoisehist, *preconv, trunc=10;
  size_t i, psf_s0, psf_s1, junk, ns0, ns1;
  size_t nc, nsize, size, hs0, hs1, *bytind;
  float *psf, psf_q=1, psf_pa=0, psf_trunc=10;

  oneprofile(0.0f, 0.0f, p->psf_p1, p->psf_p2, psf_pa, psf_q, 
	     psf_trunc, integaccu, psf_av0_tot1, psfsum, psf_smg, 
	     &psf, &psf_s0, &psf_s1, &junk);
  if(p->vpsf)
    array_to_fits("PSF.fits", NULL, "PSF", FLOAT_IMG, psf, 
		  psf_s0, psf_s1);

  hs0=psf_s0/2;          hs1=psf_s1/2;
  ns0=p->s0+2*hs0;       ns1=p->s1+2*hs1;
  size=p->s0*p->s1;		/* Shorter name ;-). */
  nsize=ns0*ns1;		/* Shorter name ;-). */
  nc=p->numppcols;		/* Shorter name ;-). */
  pp=p->profileparams;		/* Shorter name ;-). */
  ss=sqrt(p->sky);		/* Shorter name ;-). */

  assert( (img=calloc(nsize,sizeof *img))!=NULL );
  assert( (byt=calloc(nsize,sizeof *byt))!=NULL );
  assert( (bytind=malloc(nsize*sizeof *bytind))!=NULL );

  imgtopix(ns0, ns1, &D);

  if(p->verb)
    {
      printf("(x,y): profile position.\n");
      printf("\t Y : At least part of it was in the image.\n");
      printf("\t-*-: Profile was not in the image.\n");
    }

  for(i=0;i<p->nummock;i++)
    {
      suc=makeprofile(img, byt, bytind, D, ns0, 
		      ns1, trunc, integaccu, 
		      pp[i*nc+1],	  /* Profile function. */
		      pp[i*nc+3]+hs0-1,   /* x_c (C format) */
		      pp[i*nc+2]+hs1-1,   /* y_c (C format) */
		      pp[i*nc+4],	  /* p1: sersic re. */
		      pp[i*nc+5],	  /* p2: sersic n. */
		      pp[i*nc+6],	  /* position angle. */
		      pp[i*nc+7],	  /* axis ratio. */
		      ss*pp[i*nc+8],	  /* average flux.*/
		      &pp[i*nc+9]);	  /* Total flux of profile*/ 
      if(p->verb)
	printf(" - (%-.2f, %-.2f)\t%s\n", pp[i*nc+2], pp[i*nc+3], 
	       suc ? " Y" : "-*-");
    }
  if(p->verb)
    printf("\n\n");
  
  if(p->vnoconv)
    {
      floatshrinkarraytonew(img, ns0, ns1, hs0, hs1, 
			    p->s0+hs0, p->s1+hs1, &preconv);
      array_to_fits(p->outname, NULL, "NOCONV", FLOAT_IMG, 
		    preconv, p->s0, p->s1);
      free(preconv);

      if(p->verb)
	printf("- Pre-convolved profiles saved in '%s' (ext %d)\n",
	       p->outname, extcounter++);
    }

  convolve(img, ns0, ns1, psf, psf_s0, psf_s1);

  floatshrinkarray(&img, ns0, ns1, hs0, hs1, p->s0+hs0, p->s1+hs1);

  floatsetbelowtozero(img, size, MINFLOAT);

  if(p->vconv)
    {
      array_to_fits(p->outname, NULL, "NONOISE", FLOAT_IMG, img, 
		    p->s0, p->s1);
      if(p->verb)
	printf("- Convolved profiles saved in '%s' (ext %d)\n",
	       p->outname, extcounter++);
    }
   
  if(p->vhist)
    histogram(img, size, p->vhist, &p->histmin, 
	      &p->histmax, &nonoisehist, 1, 0, 0);

  addnoise(img, size, p->sky);

  array_to_fits(p->outname, NULL, "WITHNOISE", FLOAT_IMG, img, 
		p->s0, p->s1);

  if(p->verb)
    printf("- Noised image saved in '%s' (ext %d)\n",
	   p->outname, extcounter++);
   
  if(p->vhist)
    printmockhist(img, size, p->vhist, p->histmin, 
		  p->histmax, nonoisehist);

  savemockinfo(p);
  if(p->verb)
    printf("- Profile info saved in '%s'\n\n", p->infoname);

  free(psf);
  free(img);
  free(byt);
  free(bytind);
  freepixarray(D);
}
