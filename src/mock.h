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
#ifndef MOCK_H
#define MOCK_H





struct uiparams
{
  char        *infoname;      /* Name of file with galaxy info.     */
  char         *outname;      /* Name of output FITS file.          */
  char         *psfname;      /* Name of PSF FITS name.             */
  int        dontdelete;      /* ==1, don't delete an existing file.*/
  int     removenamedir;      /* Remove output name dir information.*/
};



/* Inputs into mockimg(). */
struct mockparams
{
  /* ui.c parameters: */
  struct uiparams    up;      /* Parameters only for ui.c.          */

  /* Operating modes: */
  int              verb;      /* Verbatim mode (if ==1).            */
  size_t     numthreads;      /* Number of CPU threads.             */
  
  /* Profiles: */
  float      truncation;      /* Truncation radius of the profiles. */
  float       tolerance;      /* Accuracy to stop integration.      */

  /* PSF: */
  float            *psf;      /* Point Spread Function.             */
  char      *outpsfname;      /* Output PSF name.                   */
  char       *sepsfname;      /* SExtractor's .conv PSF.            */
  int       psffunction;      /* PSF Moffat or Gaussian.            */
  float          psf_p1;      /* First paramr of PSF (FWHM).        */
  float          psf_p2;      /* Second param of PSF (Moffat beta). */
  float           psf_t;      /* PSF truncation radius.             */
  size_t         psf_s0;      /* Side length of psf along axis 0.   */
  size_t         psf_s1;      /* Side length of psf along axis 0.   */
  int           onlypsf;      /* Only view the PSF.                 */

  /* Output */
  size_t             s0;      /* C standard axis 0 size.            */
  size_t             s1;      /* C standard axis 1 size.            */
  float      background;      /* Sky value in the image.            */
  float       zeropoint;      /* Magnitude zero point.              */
  int        viewnoconv;      /* View the not convolved image.      */
  int          viewconv;      /* View the convolved image.          */
  char        *fitsname;      /* Output FITS image name.            */
  char         *catname;      /* Output catalog name.               */

  /* Internal parameters: */
  double *profileparams;      /* Table of profile parameters.       */
  size_t      numppcols;      /* Number of columns in the above.    */
  size_t        nummock;      /* Number of mock profiles.           */
};





/*Inputs of oneprofile(): 
   x_c and y_c: Only the fractional part matter, the 
                central pixel will be used. 
   p1, p2:     See above.
   pa_d:       Position angle in degrees.
   q:          Axis ratio.
   truncation: See above.
   tolerance:  Accuracy to do integration. This is the 
               point where integration will stop and the
               pixel centers will be used.
   av0_tot1:   Does "flux" (below) correspond to average or sum?
   flux:       A measure of the flux of the result.
   s0_m1_g2:   Sersic (0), Moffat (1) or Gaussian (2).
   mock:       Address of the pointer keeping the array to be made.
   x_w, y_w:   Number of rows and columns respectively of the output.
   numpixs:    Number of pixels in the object.*/
void
oneprofile(float x_c, float y_c, float p1, float p2, float pa_d, 
        float q, float truncation, float tolerance, int av0_tot1, 
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
