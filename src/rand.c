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

#include <sys/time.h>		 /* generate random seed*/
#include <gsl/gsl_rng.h>	 /* used in setrandoms*/
#include <gsl/gsl_randist.h>	 /* To make noise.*/

#include "integtwod.h"
#include "rand.h"



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
			     -90,90,      /* position angle range */
			     0.2,1,       /* Axis ration range. */
			     0.001,0.01}; /* S/N range. */
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





/* Fill pixel with random values */
float
randompoints(struct integparams *p)
{
  size_t i;
  double (*prof)(double, double, double);
  double r, q, c, s, p1, p2, co, x, y, sum=0;

  gsl_rng * rng;
  double xrange, yrange;
  const gsl_rng_type * T;

  prof=p->profile;
  xrange=p->xh-p->xl;
  yrange=p->yh-p->yl;
  p1=p->p1;        p2=p->p2;         co=p->co;
  c=p->c;          s=p->s;           q=p->q;

  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng,random_seed());
  
  for(i=0;i<NUMRANDPOINTS;i++)
    {
      x= p->xl + gsl_rng_uniform(rng)*xrange;
      y= p->yl + gsl_rng_uniform(rng)*yrange;
      r=sqrt( (x*c+y*s)*(x*c+y*s)+((y*c-x*s)*(y*c-x*s)/q/q) );
      sum+=prof(r/p1, p2, co);
    }

  gsl_rng_free(rng);

  return sum/NUMRANDPOINTS;
}
