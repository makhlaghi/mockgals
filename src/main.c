/*********************************************************************
mockgals - Detect and deblend objects.

Make any number of mock profiles in an array.

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
#include <ctype.h>
#include "attaavv.h"
#include "mock.h"

/*
 * Make 100 random profiles. 
 */
int
main(int argc, char *argv[])
{
    float *img;
    double *params=NULL;
    struct ArrayInfo intable;
    float psf_p1=3, psf_p2=3, sky=10000;
    size_t size1=200, size2=200, nummock=20;

    if(argc>1){
        if(!isdigit(argv[1][0])){
            readasciitable(argv[1], &intable);
            nummock=intable.s0;
            params=intable.d;
            if(intable.s1!=10){
                printf("\n\n\tERROR: %s doesn't have 9 columns.\n", argv[1]);
                printf("\t\tAborted\n\n");
                exit(EXIT_FAILURE);
            }
        }
        else nummock=atoi(argv[1]);
    }

    mockimg(size1, size2, sky, nummock, params, psf_p1, psf_p2, &img);

    freeasciitable(&intable);
    free(img);
    return 0;
}
