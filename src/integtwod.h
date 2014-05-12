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

#ifndef INTEGTWOD_H
#define INTEGTWOD_H

/* Integration relative error, look in GSL manual: */
#define EPSREL_FOR_INTEG 2 

/* Inputs into the 2D integration routine: integ2d();
   The parameters for the various profiles:
   - Sersic:   p1: re.   p2: n.       trunc*re;
   - Moffat:   p1: FWHM. p2: beta.    (FWHM/2)*trunc.
   - Gaussian: p1: FWHM. p2: Nothing. sigma*trunc.    
   - Point: none matter! */
struct integparams
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


double
integ2d(struct integparams *params);

void
setintegparams(int s0_m1_g2_p3, float p1, float p2, float pa_d, 
	       float q, float trunc, float *trunc_r, 
	       char *profletter, struct integparams *p);

#endif
