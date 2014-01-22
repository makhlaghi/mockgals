/*********************************************************************
arraymanip - Simple functions on arrays.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

arraymanip is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

arraymanip is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with arraymanip. If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*********************************************************************
 **********************       Make copy         **********************
 *********************************************************************/
void
floatarrcpy(float *in, size_t size, float **out)
{
    float *ipt, *opt, *fpt;
    ipt=in;
    *out=malloc(sizeof(float)*size);
    assert(*out!=NULL);
    fpt=*out+size;
    opt=*out;
    while(opt<fpt) 
        *opt++=*ipt++;
}

void
uchararrcpy(unsigned char *in, size_t size, unsigned char **out)
{
    unsigned char *ipt, *opt, *fpt;
    ipt=in;
    *out=malloc(size*sizeof(unsigned char));
    assert(*out!=NULL);
    fpt=*out+size;
    opt=*out;
    while(opt<fpt) 
        *opt++=*ipt++;
}


/*********************************************************************
 **********************    Total Sum of array   **********************
 *********************************************************************/
float 
floatarrsum(float *in, size_t size)
{
    float *pt, *fpt;
    double sum=0;
    fpt=in+size;
    pt=in;
    while(pt<fpt) 
        sum+=*pt++;
    return sum;
}

float 
floatarrsumsquared(float *in, size_t size)
{
    float *pt, *fpt;
    double sum=0;
    fpt=in+size;
    for(pt=in;pt<fpt;pt++) 
        sum+=*pt * *pt;
    return sum;
}

/* Sum over all elements of the array that are not 
covered by a mask. Any non-zero masked pixel is 
considered to be a masked pixel. */
float 
floatarrsummask(float *in, long *mask, size_t size)
{
    float *pt, *fpt;
    double sum=0;
    fpt=in+size;
    pt=in;
    while(pt<fpt) 
        if(mask==0) 
            sum+=*pt++;
    return sum;
}

/*********************************************************************
 **********************   Multiply or Sum with  **********************
 *********************************************************************/
void
floatarrmwith(float *in, size_t size, float a)
{
    float *pt, *fpt;
    fpt=in+size;
    pt=in;
    while(pt<fpt) 
        *pt++ *= a;
}

void
floatarrswith(float *in, size_t size, float a)
{
    float *pt, *fpt;
    fpt=in+size;
    pt=in;
    while(pt<fpt) 
        *pt++ += a;
}


/*********************************************************************
 **********************      Change values:     **********************
 *********************************************************************/
void
floatsetbelowtozero(float *in, size_t size, float min)
{
    float *fpt, *lfpt;
    fpt=in; 
    lfpt=in+size;
    while(fpt<lfpt) {
       if(*fpt<min) *(fpt++)=0;
       else fpt++;
    }
}

void
floatsetabovetomax(float *in, size_t size, float max)
{
    float *fpt, *lfpt;
    fpt=in; 
    lfpt=in+size;
    while(fpt<lfpt) {
       if(*fpt>max) *(fpt++)=max;
       else fpt++;
    }
}

/*********************************************************************
 **********************      shrink_array       **********************
 *********************************************************************/
/* We have a large array of size (size1*size2). We want to 
    shrink this array, such that (x1,y1) comes down to point
    (0,0) and the new array now only extends to the old 
    (x2,y2). So the size of the new array is: w1*w2 where
    w1=x2-x1 and w2=y2-y1.  */
void
floatshrinkarray(float **in, int size1, int size2,
        int x1, int y1, int x2, int y2)
{
    int temp;
    float *fpt, *start;
    size_t i, ux1, uy1, us2, numahead, w1, w2;

    assert(x1!=x2);                     /* Conditions to make sure   */
    assert(y1!=y2);                     /* the function will operate */
    if(x2<x1){temp=x1;x1=x2;x2=temp;}   /* even if bad input is fed. */
    if(y2<y1){temp=y1;y1=y2;y2=temp;}
    if(x1<0) x1=0;    if(y1<0) y1=0;
    if(x2>size1) x2=size1;
    if(y2>size2) y2=size2;
    if(x1==0 && y1==0 && x2==size1 && y2==size2) return;

    w1=x2-x1;  w2=y2-y1;
    ux1=x1; uy1=y1; us2=size2;   /* The inputs are int (can be negative,
                                     which is allowed: will become zero).
                                     but pointers are unsigned, so to 
                                     faciliate the process in the loop, 
                                     they are converted to size_t. */
    numahead=ux1*us2+uy1;         /* For the first row. */
    for(i=0;i<w1;i++){
        fpt=start=*in+i*w2;
        for(;fpt<start+w2;fpt++)
            *fpt=*(fpt+numahead);
        numahead+=us2-w2;         /* For the next rows numback increases. */
    }   
    *in=(float *)realloc(*in, w1*w2*sizeof(float));
    assert(*in!=NULL);    
}

/*********************************************************************
 **********************      Merge arrays       **********************
 *********************************************************************/
/*
 * a and b are two column arrays with "numrows" rows.
 * It is assumed that their first columns are equal.
 * This function will take the second column of b and
 * add it to a.
 */
#define NUMOUTCOLS 3
void
floatvmerge(float *a, float *b, size_t numrows, float **out)
{
   size_t i;
   float *tmp;
   
   tmp=malloc(NUMOUTCOLS*numrows*sizeof(float));
   assert(tmp!=NULL);
  
   for(i=0;i<numrows;i++){
       tmp[i*NUMOUTCOLS]=a[i*2];
       tmp[i*NUMOUTCOLS+1]=a[i*2+1];
       if(b[i*2]!=a[i*2]){
           printf("\n\n\tError in floatvmerge() (arraymanip.c)\n");
           printf("\t\tThe first columns are not equal on row %lu\n\n", i);
           exit(1);
       }
       else tmp[i*NUMOUTCOLS+2]=b[i*2+1];
   }
   *out=tmp;
}
