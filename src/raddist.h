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
#ifndef RADDIST_H
#define RADDIST_H

/* This is used to find the radial distance of a point in position
   (x,y) from a point in position (m,n) in the elliptical space
   defined by cos(theta_rad), sin(theta_rad) and q. theta_rad is the
   position angle in radians.*/
struct elraddistp
{
  double xc;			/* Center x position */
  double yc;			/* Center y position */
  double t;			/* Theta(in radians) */
  double cos;			/* cos(theta) */
  double sin;			/* sin(theta) */
  double q;			/* axis ratio */
};

void
encloseellipse(double a, double b, double theta_rad, 
        size_t *x_w, size_t *y_w);

float 
elraddist(struct elraddistp *e, double x, double y);

void
makecanvas(double a, double b, double theta_rad, float x_c, 
        float y_c, float *m_i_c, float *m_j_c, size_t *x_w, 
        size_t *y_w, float **mock);

void
smallinlarge(int l1, int l2, int s1, int s2, int xc, 
        int yc, int *x1l, int *y1l, int *x2l, int *y2l, 
        int *x1s, int *y1s, int *x2s, int *y2s);

#endif
