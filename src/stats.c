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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>

#include "stats.h"
#include "attaavv.h"
#include "forqsort.h"
#include "arraymanip.h"

/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************      Average and          ********************
 ****************    Standard deviation      ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
void
fave(float *in, size_t size, float *ave)
{
    float sum;
    sum=floatarrsum(in, size);
    *ave=sum/size;
}

void
fstd(float *in, size_t size, float *std)
{
    float sum, sum2;
    sum=floatarrsum(in, size);
    sum2=floatarrsumsquared(in, size);
    *std=sqrt( (sum2-sum*sum/size)/size );
}

/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************        Histogram          ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
void
printhist(float *in, char *filename, size_t numrows, size_t numcols)
{
    size_t i;
    FILE *fp;
    fp=fopen(filename, "w");
    for(i=0;i<numrows;i++)
        fprintf(fp, "%-10.3f%-10.0f%-10.0f\n", in[i*numcols], 
                in[i*numcols+1], in[i*numcols+2]);
    fclose(fp);   
}

/* Find the histogram of an array, the output is a two
    column array with the left bin width as the first
    column and the number of pixels in the right. */
void
histogram(float *in, size_t size, size_t numbins, 
        float min, float max, float **outhist)
{
    size_t numcols=2, i, histrow;
    float  tosubtract, *tmp, *hist, binwidth;
    
    floatarrcpy(in, size, &tmp);
    hist=calloc(numbins*numcols,sizeof(float));
    assert(hist!=NULL);

    qsort(tmp, size, sizeof(float), floatincreasing);

    if(min==max){
        min=tmp[0]; max=tmp[size-1];
    }
    binwidth=(max-min)/numbins;

    for(i=0;i<numbins;i++) 
        hist[i*numcols]=min+i*binwidth;
    for(i=0;i<numbins;i++) 
        if(hist[i*numcols]>0) break;
    tosubtract=hist[i*numcols];
    for(i=0;i<numbins;i++) 
        hist[i*numcols]-=tosubtract;

    histrow=0;
    for(i=0;i<size;i++){
        if(tmp[i]<hist[histrow*numcols]) continue;
        while (tmp[i]>=hist[(histrow+1)*numcols]) 
            if(++histrow>=numbins-1) break;
        if(histrow>numbins-1) break;
        hist[histrow*numcols+1]++;        
    }

    *outhist=hist;

    free(tmp);  
}








/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************        Sigma clip         ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/*
 * This function will repeatedly sigma clip the data 
 * and return the median. It is assumed that the data
 * have been ordered.
 */
#define REPORTSIGMACLIP 0
void
sigmaclip_converge(double *orderedarray, size_t num_elem, 
        double sigma_multiple, double accuracy,
        double *ave, double *med, double *std)
{
    double oldparam=0;
    size_t counter=0, stride=1;
    double *start, *oldstart, *dpt;

    start=orderedarray;
    while(1){
        oldstart=start;

        *med=gsl_stats_median_from_sorted_data(start, stride, num_elem);
        *ave=gsl_stats_mean(start,stride,num_elem);
        *std=gsl_stats_sd_m(start, stride, num_elem, *med);
       
        #if REPORTSIGMACLIP
        printf("%lu: mean: %.2f, med: %.2f, std= %.2f\n", 
            counter+1, *ave, *med, *std);
        if(counter>0){
            printf("\tDiff with previous: %f\n", 
                    fabs((*med-oldparam)/(*med)));
            printf("\tIndexes: Start: %lu, Number of elements: %lu\n",
                    start-orderedarray, num_elem);
        }
        #endif

        if(counter>0 && fabs((*med-oldparam)/(*med))<accuracy) 
            return;

        for(dpt=start; dpt<start+num_elem; dpt++)
            if (*dpt>*med-sigma_multiple*(*std)) {
                start=dpt; break;
            }

        for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
            if (*dpt<*med+sigma_multiple*(*std)) {
                num_elem=dpt-start+1; break;
            }

        oldparam=*med;
        counter++;
    }
    printf("\n\n\tError: sigmaclip_converge() in"); 
    printf("stats.c did not converge\n");
    exit(0);
}

/*
 * This function will do a certain number of sigma clips
 * and return the final average, median and std.
 */
void
sigmaclip_certainnum(double *orderedarray, size_t num_elem, 
        double sigma_multiple, size_t numtimes,
        double *ave, double *med, double *std)
{
    double *start, *oldstart, *dpt;
    size_t counter=0, stride=1;

    start=orderedarray;
    for(counter=0;counter<numtimes;counter++)
    {
        oldstart=start;

        *med=gsl_stats_median_from_sorted_data(start, stride, num_elem);
        *ave=gsl_stats_mean(start,stride,num_elem);
        *std=gsl_stats_sd_m(start, stride, num_elem, *med);
       
        printf("%lu: mean: %.3g, med: %.3g, std= %.3g\n", 
               counter+1, *ave, *med, *std);

        for(dpt=start; dpt<start+num_elem; dpt++)
            if (*dpt>*ave-sigma_multiple*(*std)) {
               start=dpt; break;
            }

        for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
            if (*dpt<*ave+sigma_multiple*(*std)) {
               num_elem=dpt-start+1; break;
            }

    }
}





/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************         Quantiles         ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/* 
 * Find the value corresponding to a certain quantile.
 */
void
valuefromquantile(float *data, size_t size, float quant, 
        float *quantflux)
{
    float *tmp, ave;

    assert(quant>=0 && quant<=1);

    floatarrcpy(data, size, &tmp);

    /* This is to increase the speed for a symmetric
    or positively skewed distribution, which is the case
    I have so far. */
    if(quant<=0.5){
        fave(tmp, size, &ave);
        floatsetabovetomax(tmp, size, ave);
    }

    qsort(tmp, size, sizeof(float), floatincreasing);

    *quantflux=tmp[(size_t)(quant*size)];

    if(*quantflux==ave){
         printf("\n\n\tWarning from valuefromquantile()");
         printf("in stats.c.\n\tNegatively skewed distribution");
         printf("speed will be slower!\n\n");
         free(tmp);
         floatarrcpy(data, size, &tmp);
         qsort(tmp, size, sizeof(float), floatincreasing);
         *quantflux=tmp[(size_t)(quant*size)];
    }

    free(tmp);
}

/* 
 * Find quantile corresponding to a certain value.
 */
void
quantilefromvalue(float *data, size_t size, 
        float *quant, float quantflux)
{
    float *tmp;
    size_t i;

    floatarrcpy(data, size, &tmp);

    qsort(tmp, size, sizeof(float), floatincreasing);

    for(i=0;i<size;i++) 
        if(tmp[i]>quantflux) {
            *quant=(float)i/(float)size;
            break;
        }
    if(i==size){
        printf("Input flux (%f)>largest flux (%f)\n", 
               quantflux, tmp[size-1]);
        *quant=1;
    }

    free(tmp);
}
