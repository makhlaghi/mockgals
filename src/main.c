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
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mock.h"
#include "attaavv.h"
#include "ui.h"





int
main(int argc, char *argv[])
{
  time_t rawtime;
  struct mockparams p;
  struct timeval t0, t1;

  time(&rawtime);
  gettimeofday(&t0, NULL);

  setdefaultoptions(&p);

  getsaveoptions(&p, argc, argv);

  if(p.verb)
    {
      printf("\n\n--------------------------\n");
      printf("mockgals started on %s\n", ctime(&rawtime));
    }

  readinputinfo(&p);

  checkremoveoutimage(p.outname);

  mockimg(&p);
  
  free(p.profileparams);

  if(p.verb)
    {
      gettimeofday(&t1, NULL);
      printf("mockgals finished in %.4f (seconds)\n",
	     ((double)t1.tv_sec+(double)t1.tv_usec/1e6) - 
	     ((double)t0.tv_sec+(double)t0.tv_usec/1e6));
      printf("--------------------------\n\n");
    }

  return 0;
}
