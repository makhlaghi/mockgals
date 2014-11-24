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
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "mock.h"
#include "ui.h" 		/* Needs mock.h, includes time.h and
				   sys/time.h */

int
main(int argc, char *argv[])
{
  time_t rawtime;
  struct timeval t0;
  struct mockparams p;

  time(&rawtime);	
  gettimeofday(&t0, NULL);

  setparams(&p, argc, argv, &rawtime); 

  mockimg(&p);

  freeandreporttime(&p, &t0);

  return 0;
}
