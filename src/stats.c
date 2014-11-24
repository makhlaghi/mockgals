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
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_statistics_double.h>

#include "stats.h"
#include "attaavv.h"
#include "forqsort.h"
#include "arraymanip.h"

/****************************************************************
 *****************    Mininum and Maximum    ********************
 ****************************************************************/
void
floatmin(float *in, size_t size, float *min)
{
  float tmin=FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(*in<tmin) tmin=*in;
  while (++in<fpt);
  *min=tmin;
}





void
floatmax(float *in, size_t size, float *max)
{
  float tmax=FLT_MIN, *fpt;
  fpt=in+size;
  do
    if(*in>tmax) tmax=*in;
  while(++in<fpt);
  *max=tmax;
}





void
fminmax(float *in, size_t size, float *min, float *max)
{
  float tmin=FLT_MAX, tmax=FLT_MIN, *fpt;
  fpt=in+size;
  do
    {
      if(*in>tmax) tmax=*in;
      else if(*in<tmin) tmin=*in;
    }
  while(++in<fpt);
  *max=tmax;    
  *min=tmin;    
}





void
dmax_withindex(double *in, size_t size, 
        double *max, size_t *index)
{
  size_t tindex=0;
  double *fpt, *pt=in, tmax=MINFD;

  fpt=pt+size;
  do
    if(*pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *max=tmax;
}





void
fmax_withindex(float *in, size_t size, 
	       float *max, size_t *index)
{
  size_t tindex=0;
  float *pt=in, *fpt, tmax=MINFD;

  fpt=pt+size;
  do
    if(*pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *max=tmax;
}





void
dmin_withindex(double *in, size_t size, 
	       double *min, size_t *index)
{
  size_t tindex=0;
  double *pt=in, *fpt, tmin=MAXFD;

  fpt=pt+size;
  do
    if(*pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *min=tmin;
}





void
fmin_withindex(float *in, size_t size, 
	       float *min, size_t *index)
{
  size_t tindex=0;
  float *pt=in, *fpt, tmin=MAXFD;

  fpt=pt+size;
  do
    if(*pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *min=tmin;
}




















/****************************************************************
 *****************            Sum            ********************
 ****************************************************************/
float 
floatsum(float *in, size_t size)
{
  float *fpt;
  double sum=0;
  fpt=in+size;
  while(in<fpt) 
    sum+=*in++;
  return sum;
}





float 
floatsumsquared(float *in, size_t size)
{
  float *fpt;
  double sum=0;
  fpt=in+size;
  do
    sum+=*in * *in;
  while(++in<fpt);
  return sum;
}





/* Sum over all elements of the array that are not covered by a
   mask. Any non-zero masked pixel is considered to be a masked
   pixel. */
float 
floatsummask(float *in, unsigned char *mask, 
	     size_t size, size_t *nsize)
{
  double sum=0;
  float *pt, *fpt;
  size_t counter=0;
  fpt=in+size;
    
  for(pt=in;pt<fpt;pt++) 
    if(mask[pt-in]==0)
      {
	sum+=*pt;
	counter++;
      }

  *nsize=counter;
  return sum;
}





float 
floatsummaskl(float *in, long *mask, 
              size_t size, size_t *nsize)
{
  double sum=0;
  float *pt, *fpt;
  size_t counter=0;
  fpt=in+size;
    
  for(pt=in;pt<fpt;pt++) 
    if(mask[pt-in]==0)
      {
	sum+=*pt;
	counter++;
      }
  *nsize=counter;
  return sum;
}





float 
floatsumsquaredmask(float *in, unsigned char *mask, 
                    size_t size, size_t *nsize)
{
  size_t counter=0;
  float *pt, *fpt;
  double sum=0;
  fpt=in+size;
  pt=in;
  for(pt=in;pt<fpt;pt++)
    if(mask[pt-in]==0)
      { 
	sum+=*pt * *pt;
	counter++;
      }
  *nsize=counter;
  return sum;
}





float 
floatsumsquaredmaskl(float *in, long *mask, 
                     size_t size, size_t *nsize)
{
  size_t counter=0;
  float *pt, *fpt;
  double sum=0;
  fpt=in+size;
  pt=in;
  for(pt=in;pt<fpt;pt++)
    if(mask[pt-in]==0)
      {
	sum+=*pt * *pt;
	counter++;
      }
  *nsize=counter;
  return sum;
}




















/****************************************************************
 *****************      Average and          ********************
 ****************    Standard deviation      ********************
 ****************************************************************/
void
fave(float *in, size_t size, float *ave, unsigned char *mask)
{
  float sum;
  size_t nsize;
  if(mask==NULL)
    sum=floatsum(in, size);
  else 
    {
      sum=floatsummask(in, mask, size, &nsize);
      size=nsize;
    }
  *ave=sum/size;
}





void
favel(float *in, size_t size, float *ave, long *mask)
{
  float sum;
  size_t nsize;
  if(mask==NULL)
    sum=floatsum(in, size);
  else 
    {
      sum=floatsummaskl(in, mask, size, &nsize);
      size=nsize;
    }
  *ave=sum/size;
}

/* Find the average and standard deviation of an array, assuming that
   there is a mask array. Any mask array pixel that is not zero will
   not be included in the average and standard deviations.  Here the
   mask is assumed to be unsigned char.  */
void
favestd(float *in, size_t size, float *ave, float *std, 
    unsigned char *mask)
{
  size_t nsize1, nsize2;
  float sum, sum2;
  if(mask==NULL)
    {
      sum=floatsum(in, size);
      sum2=floatsumsquared(in, size);
    }
  else 
    {
      sum=floatsummask(in, mask, size, &nsize1);
      sum2=floatsumsquaredmask(in, mask, size, &nsize2);
      assert(nsize1==nsize2);
      size=nsize1;
    }
  *ave=sum/size;
  *std=sqrt( (sum2-sum*sum/size)/size );
}





/* Similar to favestd, but when the mask is assumed to be a long
   array.  */
void
favestdl(float *in, size_t size, float *ave, float *std, 
    long *mask)
{
  size_t nsize1, nsize2;
  float sum, sum2;
  if(mask==NULL)
    {
      sum=floatsum(in, size);
      sum2=floatsumsquared(in, size);
    }
  else
    {
      sum=floatsummaskl(in, mask, size, &nsize1);
      sum2=floatsumsquaredmaskl(in, mask, size, &nsize2);
      assert(nsize1==nsize2);
      size=nsize1;
    }
  *ave=sum/size;
  *std=sqrt( (sum2-sum*sum/size)/size );
}
