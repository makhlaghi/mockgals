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

#ifndef PROFILES_H
#define PROFILES_H

double
Gaussian(double r, double junk, double a);

double
moffat_alpha(double fwhm, double beta);

double
Moffat(double rda, double nb, double junk);

double
sersic_b(double n);

double
Sersic(double rdre, double inv_n, double nb);


double
totsersic(double n, double re, double b, double q);

double
Point(double j1, double j2, double j3);

#endif