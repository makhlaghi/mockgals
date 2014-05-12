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

#include <gsl/gsl_sf_gamma.h>	 /* For total Sersic flux. */

#include "profiles.h"

/****************************************************************
 *****************         Profiles:         ********************
 ****************************************************************/
/* The Gaussian function at a point: junk is only used here to make it
   conform to the general format of functions in this file: like
   Sersic(). */
double
Gaussian(double r, double junk, double a)
{
  junk=1;
  return exp( a*r*r );
}





/* This function will find the moffat function alpha value based on
   the explantions here:
 
   http://labs.adsabs.harvard.edu/adsabs/abs/2001MNRAS.328..977T/

   alpha=(FWHM/2)/(2^(1.0/beta)-1)^(0.5). Then the moffat
   function at r is: (1.0 + (r/alpha)^2.0)^(-1.0*beta)*/
double
moffat_alpha(double fwhm, double beta)
{
    return (fwhm/2)/pow((pow(2, 1/beta)-1), 0.5f);
}





/* Find the Moffat profile for a certain radius. 
 
   rda=r/alpha     and nb=-1*b.

   This is done before hand to speed up the process. Junk is only used
   here to make it conform to the general format of functions in this
   file: like Sersic(). */
double
Moffat(double rda, double nb, double junk)
{
  junk=1; /* So it can keep the same functional form as the others. */
  return pow(junk+rda*rda, nb);
}





/* This approximation of b(n) for n>0.35 is taken from McArthur,
   Courteau and Holtzman 2003:
   http://adsabs.harvard.edu/abs/2003ApJ...582..689 */
double
sersic_b(double n)
{
  if(n<=0.35f) 
    {
      printf("\n\n   ##ERROR in mock.c's sersic_b()\n");
      printf("\t Sersic n (=%f) must be smaller than 0.35\n\n", n);
      exit(EXIT_FAILURE);
    }
  return 2*n-(1/3)+(4/(405*n))+(46/(25515*n*n))+
    (131/(1148175*n*n*n)-(2194697/(30690717750*n*n*n*n)));
}





/* Find the Sersic profile for a certain radius. rdre=r/re, inv_n=1/n,
   nb= -1*b.  */
double
Sersic(double rdre, double inv_n, double nb)
{
  return exp( nb*( pow(rdre,inv_n)-1 ) );
}




/* Find the total flux in a Sersic profile. From equation 4 in Peng
   2010. This assumes the surface brightness at the effective radius
   is 1.*/
double
totsersic(double n, double re, double b, double q)
{
  return (2*M_PI*re*re*exp(b)*n*pow(b, -2*n)*q*
	  gsl_sf_gamma(2*n));
}





/* For a point source. */
double
Point(double j1, double j2, double j3)
{
  j1=j2=j3;
  return 1;
}
