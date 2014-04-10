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
#ifndef MOCK_H
#define MOCK_H

#define EMPTYIMAGE 0
#define ONLYONEPROFILE 1
#define ONGRIDNTIMESTEN 40


/* 
   The parameters:

   Sersic: p1: n. p2: re. trunc*re;
   Moffat: p1: beta. p2: (input)fwhm, (processing)later alpha.
           (FWHM/2)*trunc.
   Gaussian: p1: sigma. p2: nothing (you can set it to zero).
           sigma*trunc.    
*/

struct mockparams
{
  double   xl;    /* lower  x boundary */
  double   xh;    /* higher x boundary */
  double    y;    /* y value when integrating over x.*/
  double   yl;    /* lower  y boundary */
  double   yh;    /* higher y boundary */
  double    c;    /* Cosine of the position angle. */
  double    s;    /* Sine of the position angle. */
  double    q;    /* axis ratio of the position angle.*/
  double pa_r;    /* Profile position angle in radians.*/
  double   p1;    /* Parameter 1, see above. */
  double   p2;    /* Parameter 2, see above. */
  double   co;    /* The constant in any profile. */
  double (*profile)(double, double, double);
};

/* 
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
   numpixs:   Number of pixels in the object.
*/
void
oneprofile(float x_c, float y_c, float p1, float p2, float pa_d, 
        float q, float trunc, float integaccu, int av0_tot1, 
        float flux, float s0_m1_g2, float **mock, size_t *x_w, 
        size_t *y_w, size_t *numpixs);

void
mockimg(size_t s0, size_t s1, float sky, size_t nummock, 
        double *prflprms, float psf_fwhm, float psf_beta, 
	int vpsf, int vnoconv, int vconv, int vhist, 
	float histmin, float histmax, float **mock, 
	char *outname, char *infoname);

#endif
