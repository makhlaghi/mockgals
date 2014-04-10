/*********************************************************************
statistics - Library of statistical functions.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

statistics is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

statistics is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with statistics. If not, see <http://www.gnu.org/licenses/>.

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
  float tmin=FLT_MAX, *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;pt++)
    if(*pt<tmin) tmin=*pt;
  *min=tmin;
}





void
floatmax(float *in, size_t size, float *max)
{
  float tmax=FLT_MIN, *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;pt++)
    if(*pt>tmax) tmax=*pt;
  *max=tmax;
}





void
fminmax(float *in, size_t size, float *min, float *max)
{
  float tmin=FLT_MAX, tmax=FLT_MIN, *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;pt++) 
    {
      if(*pt>tmax) tmax=*pt;
      else if(*pt<tmin) tmin=*pt;
    }
  *max=tmax;    
  *min=tmin;    
}





void
dmax_withindex(double *in, size_t size, 
        double *max, size_t *index)
{
  size_t tindex=0;
  double *pt=in, *fpt, tmax=MINFD;

  fpt=pt+size;
  for(;pt<fpt;pt++)
    if(*pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
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
  for(;pt<fpt;pt++)
    if(*pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
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
  for(;pt<fpt;pt++)
    if(*pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
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
  for(;pt<fpt;pt++)
    if(*pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
  *index=tindex;
  *min=tmin;
}




















/****************************************************************
 *****************            Sum            ********************
 ****************************************************************/
float 
floatsum(float *in, size_t size)
{
  float *pt, *fpt;
  double sum=0;
  fpt=in+size;
  pt=in;
  while(pt<fpt) 
    sum+=*pt++;
  return sum;
}





float 
floatsumsquared(float *in, size_t size)
{
  float *pt, *fpt;
  double sum=0;
  fpt=in+size;
  for(pt=in;pt<fpt;pt++) 
    sum+=*pt * *pt;
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




















/****************************************************************
 *****************        Histogram          ********************
 ****************************************************************/
void
printhists(float *in, char *filename, size_t numrows, size_t numcols)
{
  size_t i, j;
  FILE *fp;

  fp=fopen(filename, "w");
  for(i=0;i<numrows;i++)
    {
      fprintf(fp, "%-10.3f", in[i*numcols]);
      for(j=1;j<numcols;j++)
	fprintf(fp, "%-10.0f", in[i*numcols+j]);
      fprintf(fp, "\n");
    }
  fclose(fp);   
}





/* Find the histogram of an array, the output is a two column array
   with the left bin width as the first column and the number of
   pixels in the right.  min and max are pointers so that if they are
   chosen to be set here, the calling function can use their
   value.  

   abinonzero==1: one of the bins will start at zero.

   minq: minimum and maximum quantiles only for when the minimum and
         maximum are not already defined. For example, if minq=0.1,
         the histogram range will be betwen 0.1 quantile to the 0.9
         quantile. Set it to zero so the full range is displayed.

   n01==1: Normalized.
*/
void
histogram(float *in, size_t size, size_t numbins, 
	  float *omin, float *omax, float **outhist,
	  int abinonzero, float minq, int n01)
{
  float min, max;
  size_t numcols=2, i, histrow;
  float  tosubtract, *tmp, *hist, binwidth;
    
  floatarrcpy(in, size, &tmp);
  hist=calloc((numbins+1)*numcols,sizeof(float));
  assert(hist!=NULL);

  qsort(tmp, size, sizeof(float), floatincreasing);

  min=*omin; max=*omax;
  if(min==max)
    {
      min=tmp[(size_t)(minq*size)];
      max=tmp[(size_t)((1-minq)*(size-1))]; 
    }
  binwidth=(max-min)/numbins;

  for(i=0;i<numbins+1;i++) 
    hist[i*numcols]=min+i*binwidth;

  if(abinonzero)
    {
      for(i=0;i<numbins+1;i++) 
	if(hist[i*numcols]>0) break;
      tosubtract=hist[i*numcols];
      for(i=0;i<numbins+1;i++) 
	hist[i*numcols]-=tosubtract;
    }

  histrow=0;
  for(i=0;i<size;i++)
    {
      if(tmp[i]<hist[histrow*numcols]) continue;
      while (tmp[i]>=hist[(histrow+1)*numcols]) 
	if(++histrow>=numbins) break;
      if(histrow>numbins) break;
      hist[histrow*numcols+1]++;        
    }

  /* In case a normalized histogram is desired: */
  if(n01)
    for(i=0;i<numbins+1;i++)
      hist[i*numcols+1]/=size;

  /* Since I will be plotting with pgfplots, it is best that the
     last bin be zero: */
  hist[numbins*numcols+1]=0;
 
  *outhist=hist;
  *omin=min; *omax=max;
  free(tmp);  
}





/* will be ordered here. */
void
histogram_save(float *arr, unsigned char *mask, size_t size,
		 size_t numbins, char *outname, float minq, int n01)
{
  size_t i, nsize;
  float *forhist, omin=0, omax=0, *hist;

  floatarrcpymask(arr, size, mask, &nsize, &forhist);
  histogram(forhist, nsize, numbins, &omin, &omax, &hist, 
	    0, minq, n01);

  if(n01)
    {
      for(i=0;i<numbins+1;i++)
	hist[i*2+1]*=1000;
      printfarray(hist, numbins+1, 2, "# Histogram\n", 
		  outname, 3);
    }
  else
    printhists(hist, outname, numbins+1, 2);

  free(forhist);
  free(hist);
}





/* Take two arrays and save the two histograms in one ascii file. The
   scales will be set based on the minimum and maximum of in2.
   
   This function returns the maximum value of the histogram if it is 
   asked to be normalized. This can be used to scale a cumulative
   frequency plot that is to be over plotted on the histogram.*/
float
savetwohists(float *in1, float *in2, size_t size,
	     size_t numbins, char *filename, int abinonzero,
	     float minq, int n01)
{
  size_t i, numcols=3;
  float *hist, *hist1, *hist2;
  float min=0, max=0, maxhist=1;

  histogram(in2, size, numbins, &min, &max, &hist2, 
	    abinonzero, minq, n01);
  histogram(in1, size, numbins, &min, &max, &hist1, 
	    abinonzero, minq, n01);

  if(n01)
    for(i=0;i<numbins+1;i++)
      {
	hist1[i*2+1]*=1000;
	if(hist1[i*2+1]>maxhist) maxhist=hist1[i*2+1];
	hist2[i*2+1]*=1000;
	if(hist2[i*2+1]>maxhist) maxhist=hist2[i*2+1];
      }

  floatvmerge(hist1, hist2, numbins+1, &hist);

  if(n01)
    printfarray(hist, numbins+1, numcols, "# Histogram\n", 
		filename, 3);
  else
    printhists(hist, filename, numbins+1, numcols);

  free(hist);
  free(hist1);
  free(hist2);
  
  return maxhist;
}




















/****************************************************************
 ***************** Cumulative Frequency Plot ********************
 ****************************************************************/
/* Find the cumulative frequency plot of an array, the output is a two
   column array. The first column is the x axis (value) and the second
   column is its index. min and max are pointers so that if they are
   chosen to be set here, the calling function can use their value.

   The limit is one bin larger than the given value in order to be the
   same as the histogram function. So when you want to over plot them,
   the two have the same limits.

   minq: minimum and maximum quantiles only for when the minimum and
         maximum are not already defined. For example, if minq=0.1,
         the histogram range will be betwen 0.1 quantile to the 0.9
         quantile. Set it to zero so the full range is displayed.
*/
void
cumulativefp(float *in, size_t size, size_t numbins, 
	     float *omin, float *omax, float **outcum,
	     float minq, int n01)
{
  float min, max;
  size_t numcols=2, i, cumrow;
  float  *tmp, *cum, binwidth;
    
  floatarrcpy(in, size, &tmp);
  cum=calloc((numbins+1)*numcols,sizeof(float));
  assert(cum!=NULL);

  qsort(tmp, size, sizeof(float), floatincreasing);

  min=*omin; max=*omax;
  if(min==max)
    {
      min=tmp[(size_t)(minq*size)];
      max=tmp[(size_t)((1-minq)*(size-1))]; 
      max+=0.01*max;
    }
  binwidth=(max-min)/numbins;

  for(i=0;i<numbins+1;i++) 
    cum[i*numcols]=min+i*binwidth;
 
  cumrow=0;
  for(i=0;i<size;i++)
    if(tmp[i]>cum[cumrow*numcols])
      {
	cum[cumrow++*numcols+1]=i;
	if(cumrow>numbins)
	  break;
      }
  
  /* In case a normalized cumulative frequency plot is desired: */
  if(n01)
    for(i=0;i<numbins+1;i++)
      cum[i*numcols+1]/=size;

  /* In case the last bins are not filled (where larger than the data
     range: */
  for(i=numbins;i!=0;i--)
    if(cum[i*numcols+1]==0.0f)
      cum[i*numcols+1]=1.0f;

  *outcum=cum;
  *omin=min; *omax=max;
  free(tmp);  
}





/* will be ordered here. 

  max: the maximum value you want the cumulative distribution to
  reach. Normally it is 1. If you want to overplot it over a
  histogram, the value comes from the maximum value in the histogram.
 */
void
cumulativefp_save(float *arr, unsigned char *mask, size_t size,
		  size_t numbins, char *outname, float minq, 
		  int n01, float max)
{
  size_t i, nsize;
  float *forcum, omin=0, omax=0, *cum;

  floatarrcpymask(arr, size, mask, &nsize, &forcum);

  cumulativefp(forcum, nsize, numbins, &omin, &omax, &cum, minq, n01);

  if(n01)
    {
      for(i=0;i<numbins+1;i++)
	cum[i*2+1]*=max;
      printfarray(cum, numbins+1, 2, "# Cumulative frequency plot\n", 
		  outname, 3);
    }
  else
    printhists(cum, outname, numbins+1, 2);

  free(forcum);
  free(cum);
}





/* Take two arrays and save the two histograms in one ascii file. The
   scales will be set based on the minimum and maximum of in2. */
void
savetwocfp(float *in1, float *in2, size_t size, size_t numbins, 
	   char *filename, float minq, int n01)
{
  size_t numcols=3;
  float min=0, max=0;
  float *cum, *cum1, *cum2;

  cumulativefp(in2, size, numbins, &min, &max, &cum2, minq, n01);
  cumulativefp(in1, size, numbins, &min, &max, &cum1, minq, n01);

  floatvmerge(cum1, cum2, numbins+1, &cum);

  if(n01)
    printfarray(cum, numbins+1, numcols, "# Histogram\n", 
		filename, 6);
  else
    printhists(cum, filename, numbins+1, numcols);

  free(cum);
  free(cum1);
  free(cum2);
}




















/****************************************************************
 *****************        Sigma clip         ********************
 ****************************************************************/
/* This function will repeatedly sigma clip the data and return the
   median. It is assumed that the data have been ordered.  */
#define REPORTSIGMACLIP 0
void
sigmaclip_converge(float *array, int o1_n0, size_t num_elem, 
		   float sigma_multiple, float accuracy,
		   float *ave, float *med, float *std)
{
  size_t counter=0;
  float *start, *oldstart, *dpt;
  float oldparam=0, *orderedarray;

  if(o1_n0==0)
    {
      floatarrcpy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray, 
	    floatincreasing);      
    }
  else orderedarray=array;

  start=orderedarray;
  while(1)
    {
      oldstart=start;

      *med=*(start+num_elem/2);
      favestd(start, num_elem, ave, std, NULL);
       
      if(REPORTSIGMACLIP)
	{
	  printf("%lu: mean: %.2f, med: %.2f, std= %.2f\n", 
		 counter+1, *ave, *med, *std);
	  if(counter>0)
	    {
	      printf("\tDiff with previous: %f\n", 
		     fabs((*med-oldparam)/(*med)));
	      printf("\tIndexes: Start: %lu, Number of elements: %lu\n",
		     start-orderedarray, num_elem);
	    }
	}

      if(counter>0 && fabs((*med-oldparam)/(*med))<accuracy) 
	return;

      for(dpt=start; dpt<start+num_elem; dpt++)
	if (*dpt>*med-sigma_multiple*(*std)) 
	  {
	    start=dpt; 
	    break;
	  }

      for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
	if (*dpt<*med+sigma_multiple*(*std))
	  {
	    num_elem=dpt-start+1; 
	    break;
	  }
           

      oldparam=*med;
      counter++;
    }
  printf("\n\n\tError: sigmaclip_converge() in"); 
  printf("stats.c did not converge\n");
  exit(0);
}





/* This function will do a certain number of sigma clips and return
   the final average, median and std. o1_n0: 1: initially ordered. 2:
   initially not ordered.*/
void
sigmaclip_certainnum(float *array, int o1_n0, size_t num_elem, 
		     float sigma_multiple, size_t numtimes, 
		     float *ave, float *med, float *std)
{
  size_t counter=0;
  float *start, *oldstart, *dpt, *orderedarray;

  if(o1_n0==0)
    {
      floatarrcpy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray, 
	    floatincreasing);      
    }
  else orderedarray=array;

  start=orderedarray;
  for(counter=0;counter<numtimes;counter++)
    {
      oldstart=start;

      *med=*(start+num_elem/2);
      favestd(start, num_elem, ave, std, NULL);
       
      printf("%lu: mean: %.3g, med: %.3g, std= %.3g\n", 
	     counter+1, *ave, *med, *std);

      for(dpt=start; dpt<start+num_elem; dpt++)
	if (*dpt>*ave-sigma_multiple*(*std)) 
	  {
	    start=dpt; 
	    break;
	  }
            

      for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
	if (*dpt<*ave+sigma_multiple*(*std)) 
	  {
	    num_elem=dpt-start+1; 
	    break;
	  }
    }

  if(o1_n0==0)
    free(orderedarray);
}




















/****************************************************************
 *****************         Quantiles         ********************
 ****************************************************************/
/* Find the value corresponding to a certain quantile.  If a mask is
   supplied (mask!=NULL) then any pixel that is masked will not be
   included in the quantile finding procedure.*/ 
void
valuefromquantile(float *data, size_t size, float quant, 
        float *quantflux, unsigned char *mask)
{
  int chvalue=0;
  size_t fsize;
  float *tmp, ave=0;
  
  assert(quant>=0 && quant<=1);

  floatarrcpymask(data, size, mask, &fsize, &tmp);  

  if(quant==1)
    {
      floatmax(tmp, fsize, quantflux);
      free(tmp);
      return;
    }
     
  /* This is to increase the speed for a symmetric
     or positively skewed distribution, which is the case
     I have so far. */
  if(quant<=0.5)
    {
      chvalue=1;
      fave(tmp, fsize, &ave, mask);
      floatsetabovetomax(tmp, fsize, ave);
    }
  else if(quant>0.8)
    {
      chvalue=1;
      fave(tmp, fsize, &ave, mask);
      floatsetbelowtomin(tmp, fsize, ave);
    }

  qsort(tmp, fsize, sizeof*tmp, floatincreasing);

  *quantflux=tmp[(size_t)(quant*fsize)];

  if(chvalue==1 && *quantflux==ave)
    {
      printf("\n\n\tWarning from valuefromquantile()");
      printf("in stats.c.\n\tThe approximation did not work. ");
      printf("Speed will be slower!\n\n");
      free(tmp);
      floatarrcpymask(data, size, mask, &fsize, &tmp);
      qsort(tmp, fsize, sizeof(float), floatincreasing);
      *quantflux=tmp[(size_t)(quant*fsize-1)];
    }

  free(tmp);
}





void
multivaluefromquantile(float *data, size_t size, float *quants, 
		       float *quantfluxs, size_t numquants, 
		       unsigned char *mask)
{
  float *tmp;
  size_t i, fsize;
  
  for(i=0;i<numquants;i++)
    assert(quants[i]>=0 && quants[i]<=1);

  floatarrcpymask(data, size, mask, &fsize, &tmp);  

  qsort(tmp, fsize, sizeof*tmp, floatincreasing);

  for(i=0;i<numquants;i++)
    quantfluxs[i]=tmp[(size_t)(quants[i]*(fsize-1))];

  free(tmp);
}





/* Similar to valuefromquantile(), but the input array will not be
   copied. Since this function will sort the array, the result is that
   the input array will be sorted after this function is run.*/
void
valuefromquantile_nocopy(float *data, size_t size, 
        float quant, float *quantflux, unsigned char *mask)
{
  size_t fsize;

  assert(quant>=0 && quant<=1);

  if(quant==1)
    {
      floatmax(data, size, quantflux);
      return;
    }

  if(mask!=NULL)
    {
      removemasked(&data, size, mask, &fsize);
      size=fsize;
    }

  qsort(data, size, sizeof(float), floatincreasing);
  *quantflux=data[(size_t)(quant*size)];
}





/* Find quantile corresponding to a certain value.  If a mask is
   supplied (mask!=NULL) then any pixel that is masked will not be
   included in the quantile finding procedure.  */
void
quantilefromvalue(float *data, size_t size, 
        float *quant, float quantflux, unsigned char *mask)
{
  float *tmp;
  size_t i, fsize;

  floatarrcpymask(data, size, mask, &fsize, &tmp);

  qsort(tmp, fsize, sizeof(float), floatincreasing);

  for(i=0;i<fsize;i++) 
    if(tmp[i]>quantflux) 
      {
	*quant=(float)i/(float)fsize;
	break;
      }
  if(i==fsize)
    {
      printf("Input flux (%f)>largest flux (%f)\n", 
	     quantflux, tmp[fsize-1]);
      *quant=1;
    }

  free(tmp);
}





/* Similar to quantilefromvalue(), but the input array will not be
   copied. Since this function will sort the array, the result is
   that the input array will be sorted after this function is run.  */
void
quantilefromvalue_nocopy(float *data, size_t size, 
        float *quant, float quantflux, unsigned char *mask)
{
  size_t i, fsize;

  if(mask!=NULL) 
    {
      removemasked(&data, size, mask, &fsize);
      size=fsize;
    }

  qsort(data, size, sizeof(float), floatincreasing);

  for(i=0;i<size;i++) 
    if(data[i]>quantflux) 
      {
	*quant=(float)i/(float)size;
	break;
      }
  if(i==size)
    {
      printf("Input flux (%f)>largest flux (%f)\n", 
	     quantflux, data[size-1]);
      *quant=1;
    }
}
