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
#ifndef MOCK_H
#define MOCK_H





/* Inputs into mockimg(). */
struct mockparams
{
  char *infoname;		/* Name of file with galaxy info. */
  char *outname;		/* Name of output FITS file. */
  size_t s0;			/* C standard axis 0 size. */
  size_t s1;			/* C standard axis 1 size. */
  float sky;			/* Sky value in the image. */
  float zeropoint;		/* Magnitude zero point. */
  float trunc;			/* Truncation radius of the profiles. */
  char *psfname;		/* Name of PSF FITS name. */
  int psf_mg;			/* PSF Moffat or Gaussian. */
  float psf_p1;			/* First parameter of PSF. */
  float psf_p2;			/* Second parameter of PSF. */
  int vhist;			/* View histogram (>0) or not(0)? */
  float histmin;		/* Minimum of histogram. */
  float histmax;		/* Maximum of histogram. */
  int verb;			/* Verbatim mode (if ==1). */
  int vpsf;			/* View the PSF used. */
  int vnoconv;			/* View the not convolved image. */
  int vconv;			/* View the convolved image. */
  char *initcomments;		/* Comments of the input table. */
  double *profileparams;	/* Table of profile parameters. */
  size_t numppcols;		/* Number of columns in the above. */
  size_t nummock;		/* Number of mock profiles. */
};





/*Inputs of oneprofile(): 
   x_c and y_c: Only the fractional part matter, the 
                central pixel will be used. 
   p1, p2:    See above.
   pa_d:      Position angle in degrees.
   q:         Axis ratio.
   trunc:     See above.
   integaccu: Accuracy to do integration. This is the 
              point where integration will stop and the
              pixel centers will be used.
   av0_tot1:  Does "flux" (below) correspond to average or sum?
   flux:      A measure of the flux of the result.
   s0_m1_g2:  Sersic (0), Moffat (1) or Gaussian (2).
   mock:      Address of the pointer keeping the array to be made.
   x_w, y_w:  Number of rows and columns respectively of the output.
   numpixs:   Number of pixels in the object.*/
void
oneprofile(float x_c, float y_c, float p1, float p2, float pa_d, 
        float q, float trunc, float integaccu, int av0_tot1, 
        float flux, float s0_m1_g2, float **mock, size_t *x_w, 
        size_t *y_w, size_t *numpixs);




/* Make random profile parameters. */
void
setprflprms(double **prflprms, size_t numprflprms, 
	    size_t nummock, int size1, int size2);




/* Make any number of profiles inside an image. The inputs are
   discribed above struct mockparams. */
void
mockimg(struct mockparams *p);

#endif
