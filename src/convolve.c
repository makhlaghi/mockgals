/*********************************************************************
noisechisel - Detect and deblend objects.

Detect and deblend (segment) objects in a noisy image.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

noisechisel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

noisechisel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with noisechisel. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include <fftw3.h>

#include "convolve.h"
#include "arraymanip.h"
#include "fitsarrayvv.h"



/****************************************************************
 ***********         To check the result:        ****************
 ******   (save the complex arrays as FITS files)    ************
 ****************************************************************/
/* This function will receive a 2D complex array of type
    fftwf_complex. It will then save the desired aspect of the complex
    array into a FITS file.  The desired aspect can be: The spectrum
    (absorphase=1) and the phase angle (absorphase=0) */
void
savecomplexasfits(char *name, fftwf_complex *array, 
		  size_t s0, size_t s1, size_t absorphase)
{
  float *farray;
  char *extname;
  size_t i, size;
 
  size=s0*s1;
  farray=malloc(sizeof(float)*size);
  assert(farray!=NULL);

  if (absorphase==1)       /* In case you want the spectrum. */
    {
      for(i=0;i<size;i++)
	farray[i]=sqrt(array[i][0]*array[i][0]
		       +array[i][1]*array[i][1]);
      extname="Spectrum";
    }

  if (absorphase==0)       /* In case you want the phase angle. */
    {
      for(i=0;i<size;i++)
	farray[i]=atan2(array[i][1], array[i][0]);
      extname="Phase angle";
    }

  array_to_fits(name, NULL, extname, FLOAT_IMG, farray, s0, s1);
  free(farray);
}




















/****************************************************************
 *****************   Background function(s)  ********************
 ****************************************************************/
/* Prior to convolution, both arrays have to be padded to 
    the same size, that is the job of this function. */
#define PADDEDCHECK 0
void
makepadded(float *f, size_t fs0, size_t fs1, 
           float *h, size_t hs0, size_t hs1, 
           float **tf, float **th, size_t *ps0, size_t *ps1)
{
  float *ttf, *tth;
  size_t i, j, num_pix, num_pix1, num_pix2;

  *ps0=fs0+hs0-1;
  *ps1=fs1+hs1-1;

  if (*ps0%2==1) (*ps0)++;
  if (*ps1%2==1) (*ps1)++;

  num_pix = *ps0 * *ps1;
  num_pix1=fs0*fs1;
  num_pix2=hs0*hs1;

  assert( (*tf=calloc(num_pix, sizeof **tf))!=NULL );
  assert( (*th=calloc(num_pix, sizeof **th))!=NULL );

  ttf=*tf;
  for(i=0;i<fs0;i++)
    for(j=0;j<fs1;j++)
      ttf[ i* *ps1 + j ]=f[i*fs1+j];

  tth=*th;
  for(i=0;i<hs0;i++)
    for(j=0;j<hs1;j++)
      tth[ i* *ps1 + j ]=h[i*hs1+j];
   
  if(PADDEDCHECK)
    {
      array_to_fits("padded.fits", NULL, "Padded Image", 
		    FLOAT_IMG, *tf, *ps0, *ps1);
      array_to_fits("padded.fits", NULL, "Padded Kernel", 
		    FLOAT_IMG, *th, *ps0, *ps1);
    }
}




















/****************************************************************
 *****************      Main Function:      ********************
 ****************************************************************/
/* This function get two arrays: f(img) and h(psf) and convolves them,
    the result will be saved in c(convolved).  Space for output (o)
    will be allocated in this function. */
void
convolve(float *f, size_t fs0, size_t fs1, 
         float *h, size_t hs0, size_t hs1,
	 float **o)
{
  fftwf_plan p;
  fftwf_complex *ch, *cf;
  size_t i, size, ps0, ps1;
  float *tf, *th, *tout, temp;

  if(hs0%2==0 || hs1%2==0)
    {
      printf("\n\nThe convolution kernel has to have odd sides.\n");
      printf("the given kernel has size: %lux%lu\n\n\n", hs0, hs1);
      exit(EXIT_FAILURE);
    }

  makepadded(f, fs0, fs1, h, hs0, hs1, &tf, &th, &ps0, &ps1);

  cf=fftwf_malloc(sizeof(fftwf_complex)*ps0*(ps1/2+1));
  p=fftwf_plan_dft_r2c_2d(ps0, ps1, tf, cf, FFTW_ESTIMATE);
  fftwf_execute(p); 

  ch=fftwf_malloc(sizeof(fftwf_complex)*ps0*(ps1/2+1));
  p=fftwf_plan_dft_r2c_2d(ps0, ps1, th, ch, FFTW_ESTIMATE);
  fftwf_execute(p);

  if(DFTCHECK)
    {
      savecomplexasfits("dft.fits", cf, ps0, ps1/2+1, 1);
      savecomplexasfits("dft.fits", ch, ps0, ps1/2+1, 1);
    }

  size=ps0*(ps1/2+1);
  for(i=0;i<size;i++)
    {
      temp=cf[i][0]*ch[i][0]-cf[i][1]*ch[i][1]; 
      cf[i][1]=cf[i][0]*ch[i][1]+cf[i][1]*ch[i][0]; 
      cf[i][0]=temp;
    }

  if(DFTCHECK)
    savecomplexasfits("dft.fits", cf, ps0, ps1/2+1, 1);

  tout=malloc(ps0*ps1*sizeof *tout);
  p=fftwf_plan_dft_c2r_2d(ps0, ps1, cf, tout, FFTW_ESTIMATE);
  fftwf_execute(p);  

  floatshrinkarraytonew(tout, ps0, ps1, hs0/2, hs1/2,
			fs0+hs0/2, fs1+hs1/2, o);

  free(tf);       
  free(th);      
  free(tout);
  fftwf_free(cf);         
  fftwf_free(ch);
  fftwf_destroy_plan(p);                
  fftwf_cleanup();
}
