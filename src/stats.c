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

void
quantilefromvalue(float *data, size_t size, 
        float *quant, float quantflux)
/* 
 * Find quantile corresponding to a certain value.
 */
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


#if 0
/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************           Mode            ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
#define PRINTCM 0
int
check_with_mirrored(float *ordered, size_t NFF_index, 
        size_t size, size_t maskcounter, size_t steps)
/*
 * Based on the NFF value, this function will divide the
 * data to regions above and below the NFF value (not
 * including the NFF value its self). It will then 
 * return the difference of the mirrored negative and   
 * actual positive histogram bins                      
 */
{
    /* Declarations: */
    float NFF;
    int diff, minneg=0;
    size_t i, j, NFFpi, nomaskNFFindex;

    #if PRINTCM
    int junk, ID=174;
    char command[200];
    system("reset");
    printf("\n\n\n\n\n\n");
    printf("%-8s %-8s %-15s %-10s %-10s\n", "p_ind", "mir_ind", 
            "F (NFF subt)", "Difference", "Counting Error");
    #endif

    NFF=ordered[NFF_index];
    minneg=-300;/*-1*sqrt(NFF_index)/10;*/

    /* The masked pixels were set to have the minimum value
    of the array, so in the ordered array, they take up the
    first "maskcounter" elements. So "j" which counts the 
    backwards steps that we take must not become larger than
    nomaskNFFindex. */
    nomaskNFFindex=NFF_index-maskcounter;

    /* Go over the above and bellow regions of the NFF symmetrically. */
    j=0;i=1;
    NFFpi=NFF_index+i;
    for(;NFFpi<size;i+=steps)
    {

        /* This is to increase speed: */
        NFFpi=NFF_index+i;

        /* Lets call: 
                NDiff=NFF-ordered[NFF_index-j];
                PDIFF=ordered[NFFpi]-NFF;
            Since we start with j=0 and i=1, we start with: NDIFF=0 and PDIFF>0. 
            Note: the array is ordered. Therefore at starting point: NDIFF<PDIFF.
            So the moment NDIFF>PDIFF, stop increasing j.   */
        while(NFF-ordered[NFF_index-j]<ordered[NFFpi]-NFF && j<nomaskNFFindex) 
            j++;
        if(j>nomaskNFFindex) break;

        /* Check to see if j should be increased by one pixel. 
        note that: 
         ----> increasing flux and j:
        (1)NFF-[NFF_index-j]  >  (2)[NFFpi]-NFF  >  (3)NFF-[NFF_index-j+1] 

        The question is: which one of (1) or (3) is closer to (2).*/
        if( NFF-ordered[NFF_index-j]   -   ordered[NFFpi]-NFF  > 
            ordered[NFFpi]-NFF   -   NFF-ordered[NFF_index-j+1])
            /* If (3) is closer to (2) than (1) then decrease j: */
            j--;

        /* Set the difference in number between the two: */
        diff=((int)i-(int)j);

        #if PRINTCM
        printf("%-8lu %-8lu %-15f %-+10d\n", i, j, 
                ordered[NFFpi]-NFF, (int)i-(int)j);
        #endif

        if (diff<minneg ) {diff=-1*size; break;}
    }
    if(diff==-1*(int)size) diff=-1;
    else diff=1;

    #if PRINTCM
    sprintf(command, "./src/PythonPlots/MirroredNeg.py ./PS/%d.fits %.10f %d", 
            ID, NFF, ID);
    system(command);    
    printf("Continue? (Type an int)  ");
    scanf("%d", &junk);
    printf("\n\n");
    #endif

    /* Note that minneg is negative. */
    return diff;
}

size_t
findNFFindex(float *ordered, size_t low, size_t mid, 
        size_t high, float tol, size_t size, 
        size_t maskcounter, size_t steps)
/* 
 * This function will recursively try to find the 
 * index of NFF, by dividing the region by 2. 
 * Initially the state is (CWM: check_with_mirrored)
 * The value shown is the sign of CWM output.
 *
 *  Position:   low       mid        high
 *    CWM:       +         ?          -
 * 
 * We want the least possible positive value and we know
 * that the higher position will always be more negative
 * (or less positive).
 * So if mid>0: check the point in the middle mid and high
 *    if mid<0: check the point in the middle low and mid
 */
{
    /*printf("low: %lu, mid=%lu, high=%lu\n", low, mid, high);*/
    if(high-low<tol*mid) return mid;
    
    if(check_with_mirrored(ordered, mid, size, maskcounter, steps)>0)
        return findNFFindex(ordered, mid, mid+(high-mid)/2, 
                high, tol, size, maskcounter, steps);
    else
        return findNFFindex(ordered, low, low+(mid-low)/2, 
                mid, tol, size, maskcounter, steps);
}

float 
sym_mode_by_cdf(float *data, short *mask, size_t size, size_t steps)
/*
 * Find the symmetric mode of a distribution. By symmetric mode,
 * it is meant the point where the distribution is most symmetric about.
 * It assumes the mode is below the 51 percentile of the image.
 * data: 
 *       A one dimentional array, keeping the data+noise.
 *
 * mask: 
 *       A one dimentional array, keeping any data pixels
 *       that must not be processed. If all pixels are to
 *       be considered, set mask to NULL.
 *
 * size: 
 *       The size of both the input arrays.
 *
 * steps:
 *       The pixels that will be checked, if 1, all pixels
 *       around the NFF will be checked, if larger, some 
 *       pixels will be ignored.
 */
{
    /* Declarations: */
    int diff=0;
    size_t maskcounter=0, NFFi=0;
    char funcname[]="sym_mode_by_cdf()";
    size_t j=0, stride=1, high=0, low=0;
    float *actual_index=NULL, *ordered=NULL;
    float i=0, quant=0, tol=0.01, mindata=1e10;

    /* Prepare for the mask by changing the masked values 
    in the data to the global minimum of the image. */
    if(mask!=NULL)
    {
        for(j=0;j<size;j++)if(data[j]<mindata)mindata=data[j];
        for(j=0;j<size;j++)if(mask[j]!=0)
             {maskcounter++;data[j]=mindata;}
    }

    /* Make the actual index array, this will be	
    used to approximate the mirrored negative cdf
    in make_weight. */
    actual_index=malloc(size*sizeof(float));
    ordered=malloc(size*sizeof(float));
    if(actual_index==NULL) 
        array_alloc_err(filename, funcname, "actual_index", size);
    if(ordered==NULL) 
        array_alloc_err(filename, funcname, "ordered", size);


    /* Fill in the two created arrays, i is float,
    j is size_t, that is why they are incremented 
    separately. This is done to not waste computer
    energy on type conversion. */
    i=0;j=0;
    while(j<size)
    {
        ordered[j]=data[j];
        actual_index[j]=i;
        i++;j++;
    }

    /* Order actual_pixels and 
    actual_index based on the flux: */
    gsl_sort2_float(ordered, stride, actual_index, stride, size);

    /* Find the region (with 0.1% quantile accuracy) 
    that contains the lowest. */
    quant=0.51;
    while(quant>=0)
    {
        j=maskcounter+(size_t)((size-maskcounter)*quant);
        diff=check_with_mirrored(ordered, j, size, maskcounter, steps);
        if(diff<=0) {high=j; quant-=0.02; continue;}
        else {low=j;break;}
    }
    if(quant<=0) 
        printerrorquit(filename, funcname, "Lower bound not found.");

    /* This comes from the golden section minimization: */
    NFFi=findNFFindex(ordered, low, low+(size_t)((high-low)/2),  
            high, tol, size, maskcounter, steps);

    /* Put the result in NFF: */
    if (NFFi==0)
        printerrorquit(filename, funcname, 
                "Mode index (NFFi) should not be zero.");
    *NFF=ordered[NFFi];

    free(ordered);
    free(actual_index);
}
#endif
