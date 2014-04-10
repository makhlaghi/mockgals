/*********************************************************************
raddist - radial distance of all pixels in an array.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

raddist is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

raddist is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with raddist. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "raddist.h"

/* Any ellipse can be enclosed into a rectangular box. The purpose of
   this function is to give the height and width of that box. The
   logic behind it is this: All the points on the circumference of an
   ellipse that is aligned on the x axis can be written as:

   (acos(t),bsin(t)) where 0<t<2\pi.  (1)
  
   But when we rotate the ellipse by \theta, the points can be
   characterized by:

   ( acos(t)cos(\theta)+bsin(t)sin(\theta),        (2)
    -acos(t)sin(\theta)+bsin(t)cos(\theta) )

   To find the maximum and minimum points of this function you just
   have to take the derivative of each with respect to "t" and set it
   to zero. This will give you the "t" that maximizes both x and the
   "t" that maximizes y.  Once you do that, you will get:

   For x: tan(t)=(b/a)tan(\theta)                  (3)
   For y: tan(t)=(-b/a)cot(\theta)

   Once you find the "t", put it in (2) for the respective coordinate
   and you will find the distance (about the center of the ellipse
   that encloses the whole ellipse.  */
void
encloseellipse(double a, double b, double theta_rad, 
        size_t *x_w, size_t *y_w)
{
  double t_x, t_y, max_x, max_y;
  t_x=atan(b/a*tan(theta_rad));
  t_y=atan(-1*b/a/tan(theta_rad));

  max_x=a*cos(t_x)*cos(theta_rad)+b*sin(t_x)*sin(theta_rad);
  max_y=-1*a*cos(t_y)*sin(theta_rad)+b*sin(t_y)*cos(theta_rad);
 
  /* max_x and max_y are calculated from the 
     center of the ellipse. We want the final height
     and width of the box enclosing the ellipse. So 
     we have to multiply them by two, then take one 
     from them (for the center).*/
  *x_w=2*( (size_t)fabs(max_x)+1 ) - 1;
  *y_w=2*( (size_t)fabs(max_y)+1 ) - 1;
}





float 
elraddist(struct elraddistp *e, double x, double y)
{
  double x_n, y_n;		/* For readability */
  x_n = (e->xc-x)*e->cos + (e->yc-y)*e->sin;
  y_n = (e->yc-y)*e->cos - (e->xc-x)*e->sin;
  return sqrt(x_n*x_n+y_n*y_n/e->q/e->q);
}





/* Find the correct size of the output image of an elliptical profile
   and malloc the array to keep the elliptical profile. Don't forget
   to free it later.  It is made for radial profiles, so it will fill
   the array with the distance to the center (x_c,y_c).  this will
   make applying a radial profile really easy.  This function is
   defined separately so that it can be used for a wide variety of
   profiles.
 
   Input parameters:
   x_c, y_c: The center of the profile in a final array 
             The final array is not necessarily the one 
             we are making here. The final array might
             contain several mock profiles. Here we are
             just making one of them. In short only the 
             pixel fractions are important here.
   a,b:      The semi-major and semi-minor axis of the 
             ellipse
   theta_rad:The position angle of the ellipse in 
             units of radians.
  
   Output parameters (pointers):
   x_w, y_w:        width of the final image in the x and 
                    y directions.
   m_i_c, m_j_c:    The profile's central point in this grid.
   mock:            Pointer to the pointer keeping the 
                    allocated space for the mock array.  */
void
makecanvas(double a, double b, double theta_rad, float x_c, 
        float y_c, float *m_i_c, float *m_j_c, size_t *x_w, 
        size_t *y_w, float **mock)
{
  float *pmock; /* to speed up the job. */
  int i,j,i_w,j_w;
  double x_n, y_n, diff_x, diff_y;
  register double co, si, i_c, j_c, q;

  q=b/a;

  encloseellipse(a, b, theta_rad, x_w, y_w);

  diff_x=x_c-(int)x_c;
  if(diff_x>0.5) {x_c++; diff_x=1-diff_x;}
  diff_y=y_c-(int)y_c;
  if(diff_y>0.5) {y_c++; diff_y=1-diff_y;}

  i_c=*x_w/2+diff_x;
  j_c=*y_w/2+diff_y;

  *mock=malloc( (*x_w) * (*y_w) * sizeof(float));
  assert(*mock!=NULL);

  pmock=*mock;
  co=cos(theta_rad);
  si=sin(theta_rad);
  i_w=*x_w; j_w=*y_w;
  for(i=0;i<i_w;i++)
    for(j=0;j<j_w;j++)
      {
	x_n=(i_c-i)*co+(j_c-j)*si;
	y_n=(j_c-j)*co-(i_c-i)*si;
	pmock[i*j_w+j]=sqrt(x_n*x_n+y_n*y_n/q/q);
      }

  *m_i_c=i_c;
  *m_j_c=j_c;
}





/* Find the position of a smaller image in a larger one when we know
   that the central pixel of the smaller is in position (xc, yc) of
   the larger one.  In order to undertand the outputs, consider the
   two points, we want their coordinates in both images: (1) The
   bottom left corner of the smaller array.  (2) The top right corner
   of the smaller array.  Don't forget that we are assuming the array
   element [0] is positioned on the bottom left.

   l1,l2:     The number of rows and columns in the larger.
   s1,l2:     The number of rows and columns in the smaller.
   xc,yc:     The position of the central pixel of the 
                 smaller array in the larger one. 
   x1l, y1l:  Coordinates of point (1) in the larger array.
   x2l, y2l:  Coordinates of point (2) in the larger array.
   x1s, y1s:  Coordinates of point (1) in the smaller array.
   x2s, y2s:  Coordinates of point (2) in the smaller array.*/
#define CHECKSMALLINLARGE 0
void
smallinlarge(int l1, int l2, int s1, int s2, int xc, 
        int yc, int *x1l, int *y1l, int *x2l, int *y2l, 
        int *x1s, int *y1s, int *x2s, int *y2s)
{
  assert(s1%2==1 && s2%2==1);

  *x1l=xc-s1/2;    *y1l=yc-s2/2;
  *x2l=xc+s1/2;    *y2l=yc+s2/2;

  *x1s=0;     *y1s=0;
  *x2s=s1;    *y2s=s2;

  if(CHECKSMALLINLARGE)
    {
      printf("l1: %d, l2: %d\n", l1, l2);
      printf("s1: %d, s2: %d\n", s1, s2);
      printf("xc: %d, yc: %d\n", xc, yc);
      printf("\nInitial positions:\n");
      printf("(x1l,y1l)=(%d,%d)\n", *x1l, *y1l);
      printf("(x1s,y1s)=(%d,%d)\n", *x1s, *y1s);
      printf("(x2l,y2l)=(%d,%d)\n", *x2l, *y2l);
      printf("(x2s,y2s)=(%d,%d)\n", *x2s, *y2s);
    }


  if(*x1l<0)  {*x1s=-1*(*x1l); *x1l=0;}
  if(*y1l<0)  {*y1s=-1*(*y1l); *y1l=0;}
  if(*x2l>l1) {*x2s=s1-(*x2l-l1)-1; *x2l=l1;}
  if(*y2l>l2) {*y2s=s2-(*y2l-l1)-1; *y2l=l2;}

  if(CHECKSMALLINLARGE)
    {
      printf("\nFinal positions:\n");
      printf("(x1l,y1l)=(%d,%d)\n", *x1l, *y1l);
      printf("(x1s,y1s)=(%d,%d)\n", *x1s, *y1s);
      printf("(x2l,y2l)=(%d,%d)\n", *x2l, *y2l);
      printf("(x2s,y2s)=(%d,%d)\n", *x2s, *y2s);
    }
}
