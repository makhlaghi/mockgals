/*********************************************************************
freqdomain - Frequency domain applications on C 2d arrays.

Copyright (C) 2013 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

freqdomain is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

freqdomain is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with freqdomain. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
/* System libraries: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include <fftw3.h>

/* Other libraries: */
#include "mock.h"
#include "arraymanip.h"
#include "freqdomain.h"
#include "fitsarrayvv.h"



/****************************************************************
 ****************************************************************
 *****************                           ********************
 ***********         To check the result:        ****************
 ***********    (save the arrays as FITS files   ****************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/* This function will receive a 2D complex array of type 
    fftw_complex. It will then save the desired aspect of 
    the complex array into a FITS file.                    
    The desired aspect can be: The spectrum (absorphase=1) 
    and the phase angle (absorphase=0) */
void
savecomplexasfits(char *name, fftw_complex *array, size_t size1, 
        size_t size2, size_t absorphase)
{
    size_t i, size;
    float *farray;
 
    size=size1*size2;
    farray=malloc(sizeof(float)*size);
    assert(farray!=NULL);

    if (absorphase==1)       /* In case you want the spectrum. */
        for(i=0;i<size;i++)
            farray[i]=sqrt(array[i][0]*array[i][0]+array[i][1]*array[i][1]);

    if (absorphase==0)       /* In case you want the phase angle. */
        for(i=0;i<size;i++)
            farray[i]=atan2(array[i][1], array[i][0]);

    array_to_fits(name, NULL, "", "FLOAT", farray, size1, size2);

    free(farray);
}

/* Saving a double array is hardly ever necessary:
    it will just take too much space and that level of
    accuracy is hardly ever needed. */
void
savedoubleasfits(char *name, double *array, 
        size_t size1, size_t size2)
{
    size_t i,j;
    float *farray;

    farray=malloc(sizeof(float)*size1*size2);
    assert(farray!=NULL);

    for(i=0;i<size1;i++)
        for(j=0;j<size2;j++) 
            farray[i*size2+j]=array[i*size2+j];

    array_to_fits(name, NULL, "", "FLOAT", farray, size1, size2);

    free(farray);
}

/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************   Background function(s)  ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/* Prior to convolution, both arrays have to be padded to 
    the same size, that is the job of this function. */
void
makepadded(float *f, size_t fsize1, size_t fsize2, 
           float *h, size_t hsize1, size_t hsize2, 
           double **tempf, double **temph, 
           size_t *paddedsize1, size_t *paddedsize2)
{
    double *ttempf, *ttemph;
    size_t i, j, num_pix, num_pix1, num_pix2;

    *paddedsize1=fsize1+hsize1-1;
    *paddedsize2=fsize2+hsize2-1;
    if (*paddedsize1%2==1) (*paddedsize1)++;
    if (*paddedsize2%2==1) (*paddedsize2)++;
    num_pix = *paddedsize1 * *paddedsize2;
    num_pix1=fsize1*fsize2;
    num_pix2=hsize1*hsize2;
    *tempf=calloc(num_pix, sizeof(double));
    assert(*tempf!=NULL);
    *temph=calloc(num_pix, sizeof(double));
    assert(*temph!=NULL);

    ttempf=*tempf;
    for(i=0;i<fsize1;i++)
        for(j=0;j<fsize2;j++)
            ttempf[ i* *paddedsize2 + j ]=f[i*fsize2+j];

    ttemph=*temph;
    for(i=0;i<hsize1;i++)
        for(j=0;j<hsize2;j++)
            ttemph[ i* *paddedsize2 + j ]=h[i*hsize2+j];
   
    #if CONVOLVECHECK
    savedoubleasfits("padded.fits", *tempf, 
        *paddedsize1, *paddedsize2);
    savedoubleasfits("padded.fits", *temph, 
        *paddedsize1, *paddedsize2);
    #endif
}


/****************************************************************
 ****************************************************************
 *****************                           ********************
 *****************      Main Functions:      ********************
 *****************                           ********************
 ****************************************************************
 ****************************************************************/
/* This function get two arrays: f(img) and h(psf) and 
    convolves them, the result will be saved in c(convolved).
    Space for c will be allocated in this function. */
void
convolve(float *f, size_t fsize1, size_t fsize2, 
         float *h, size_t hsize1, size_t hsize2)
{
    fftw_plan p;
    fftw_complex *ch, *cf;
    size_t i, j, size, psize1, psize2;
    double *tempf, *temph, *tempc, temp;

    makepadded(f, fsize1, fsize2, h, hsize1, hsize2, 
           &tempf, &temph, &psize1, &psize2);

    cf=fftw_malloc(sizeof(fftw_complex)*psize1*(psize2/2+1));
    p=fftw_plan_dft_r2c_2d(psize1, psize2, tempf, cf, FFTW_ESTIMATE);
    fftw_execute(p); 

    ch=fftw_malloc(sizeof(fftw_complex)*psize1*(psize2/2+1));
    p=fftw_plan_dft_r2c_2d(psize1, psize2, temph, ch, FFTW_ESTIMATE);
    fftw_execute(p);

    #if CONVOLVECHECK
    savecomplexasfits("dft.fits", cf, psize1, psize2/2+1, 1);
    savecomplexasfits("dft.fits", ch, psize1, psize2/2+1, 1);
    #endif

    size=psize1*(psize2/2+1);
    for(i=0;i<size;i++){
        temp=cf[i][0]*ch[i][0]-cf[i][1]*ch[i][1];
        cf[i][1]=cf[i][0]*ch[i][1]+cf[i][1]*ch[i][0];
        cf[i][0]=temp;
    }

    #if CONVOLVECHECK
    savecomplexasfits("dft.fits", cf, psize1, psize2/2+1, 1);
    #endif

    tempc=malloc(psize1*psize2*sizeof(double));
    p=fftw_plan_dft_c2r_2d(psize1, psize2, cf, tempc, FFTW_ESTIMATE);
    fftw_execute(p);  

    for(i=0;i<fsize1;i++)
        for(j=0;j<fsize2;j++) 
            f[i*fsize1+j]=
               tempc[(i+hsize1/2)*psize2+(j+hsize2/2)]/(psize1*psize2);

    free(tempf);       free(temph);      free(tempc);
    fftw_free(cf);     fftw_free(ch);
    fftw_destroy_plan(p);
}

/*
 * This function will make the PSF based on the input
 * and use that to convolve the image. For the parameters
 * see mock.h.
 */
#define CONVFUNCCHECK 0
void
convolve_function(float *f, size_t fsize1, size_t fsize2,
        float p1, float p2, float pa_d, float q, float trunc, 
        float integaccu, float s0_m1_g2, float **conv)
{
    float *psf;
    int av0_tot1=1; 
    float sumflux=1;
    size_t numpixs, x_w, y_w;   

    oneprofile(0.0, 0.0, p1, p2, pa_d, q, trunc, integaccu, 
            av0_tot1, sumflux, s0_m1_g2, &psf, &x_w, &y_w, 
            &numpixs);

    #if CONVFUNCCHECK
    array_to_fits("conv.fits", NULL, "CONVFUNC", "FLOAT", 
            psf, x_w, y_w);
    #endif

    floatarrcpy(f, fsize1*fsize2, conv);

    convolve(*conv, fsize1, fsize2, psf, x_w, y_w);
}
