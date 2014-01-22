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
#ifndef RADDIST_H
#define RADDIST_H

void
encloseellipse(double a, double b, double theta_rad, 
        size_t *x_w, size_t *y_w);

void
makecanvas(double a, double b, double theta_rad, float x_c, 
        float y_c, float *m_i_c, float *m_j_c, size_t *x_w, 
        size_t *y_w, float **mock);

void
smallinlarge(int l1, int l2, int s1, int s2, int xc, 
        int yc, int *x1l, int *y1l, int *x2l, int *y2l, 
        int *x1s, int *y1s, int *x2s, int *y2s);

#endif
