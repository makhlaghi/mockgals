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

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_integration.h> /* gsl_integration_qng*/

#include "profiles.h"
#include "integtwod.h"

/****************************************************************
 *****************      2D integration:      ********************
 ****************************************************************/
/* Calculate a profile on position (i,j) in a cartesian
   grid. While the ellipse of the profile is rotated and has axis
   ratio q. i,j: positions on the grid. */
double
rot_ell(double x, struct integparams *p)
{
  double c=p->c, s=p->s;
  double r;
  r=sqrt((x*c+p->y*s)*(x*c+p->y*s)+
	 ((p->y*c-x*s)*(p->y*c-x*s)/p->q/p->q));
  return p->profile(r/p->p1, p->p2, p->co);
}





/* This function is used in the integration of a profile. It
   assumes a fixed y and integrates over a range of x values.  */
double
twod_over_x(double x, void *params)
{
  struct integparams *p;
  p=(struct integparams *)params;

  return rot_ell(x, params);
}





/* Find the 2d integration over the region. */
double
twod_over_xy(double y, void *params)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  struct integparams *p;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result;

  F.function = &twod_over_x;
  F.params = params;

  p=(struct integparams *)params;
  p->y=y;
  gsl_integration_qng(&F, p->xl, p->xh, epsabs, epsrel, 
		      &result, &abserr, &neval);
  return result;
}




/* 2D integration of a profile.*/
double
integ2d(struct integparams *params)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result; 

  F.function = &twod_over_xy;
  F.params = params;   
  gsl_integration_qng(&F, params->yl, params->yh, epsabs, 
		      epsrel, &result, &abserr, &neval);
  return result;
}





/* Set the general parameters of struct integparams. In all profiles,
   p1 is the parameter that the radius is divided by. */
void
setintegparams(int s0_m1_g2_p3, float p1, float p2, float pa_d, 
	       float q, float trunc, float *trunc_r, char *profletter, 
	       struct integparams *p)
{
  switch(s0_m1_g2_p3)
    {
    case 0: /* Sersic: p1: re, p2: n */
      *trunc_r=p1*trunc;	/* Trunctation in units of re. */
      p->p1=p1; 		/* re. */
      p->p2=1/p2;		/* n is used as 1/n in Sersic. */
      p->co=-1*sersic_b(p2);	/* Constant in the power. */
      p->profile=&Sersic;
      *profletter='s';
      break;
    case 1: /* Moffat: p1: fwhm(input) -> alpha. p2: beta*/
      *trunc_r=(p1/2)*trunc;	/* Truncation in units of FWHM/2. */
      p->p1=moffat_alpha(p1, p2);
      p->p2=-1*p2;              /* Beta is always negative!  */
      p->co=0;	                /* No constant terms needed */
      p->profile=&Moffat;
      *profletter='m';
      break;
    case 2: /* Gaussian: p1: FWHM(input) -> sigma */
      p1/=2.35482;              /* Convert FWHM to sigma. */
      *trunc_r=p1*trunc;	/* Truncation in units of sigma */
      p->co=-1.0f/(2.0f*p1*p1);	/* Constant to multiply */
      p->p1=1;			/* r gets divided by p->p1! */
      p->p2=0;			/* Not needed! */
      p->profile=&Gaussian;   
      *profletter='g';
      break;
    case 3: /* Point:*/
      *trunc_r=2;
      p->p1=p->p2=0;
      p->co=1;
      p->profile=&Point;
      *profletter='p';
      break;
    default:
      printf("\n\n\tError: s0_m1_g2 must be 0, 1 or 2\n");
      printf("\t\tIt is: %d\n\n", s0_m1_g2_p3);
      exit(EXIT_FAILURE);
      break;
    }
  p->pa_r=pa_d*M_PI/180;
  p->c=cos(p->pa_r);
  p->s=sin(p->pa_r);
  p->q=q;
}
