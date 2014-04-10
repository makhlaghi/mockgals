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

#include "stats.h"

/****************************************************************
 *****************   Initialize an array     ********************
 ****************************************************************/
void
longinit(long *in, size_t size, const long v)
{
  long *pt=in, *fpt;
  fpt=in+size;
  while(pt<fpt) *pt++=v;
}





void
floatinit(float *in, size_t size, const float v)
{
  float *pt=in, *fpt;
  fpt=in+size;
  while(pt<fpt) *pt++=v;
}





void
ucharinit(unsigned char *in, size_t size, 
          const unsigned char v)
{
  unsigned char *pt=in, *fpt;
  fpt=in+size;
  while(pt<fpt) *pt++=v;
}




















/****************************************************************
 *****************       Change values       ********************
 ****************************************************************/
void 
changefvalueinarray(float *in, size_t size, 
        float from, float to)
{
  float *pt=in, *fpt;
  fpt=in+size;    
  for(;pt<fpt;pt++)
    if(*pt==from) *pt=to;
}





void 
truncfarray(float *in, size_t size, float thresh)
{
  float *pt=in, *fpt;
  fpt=in+size;    
  for(;pt<fpt;pt++)
    if(*pt>thresh) *pt=thresh;
}





void 
setabovetozerofarray(float *in, size_t size, float thresh)
{
  float *pt=in, *fpt;
  fpt=in+size;    
  for(;pt<fpt;pt++)
    if(*pt>thresh) 
      *pt=0;
}





void 
changelvalueinarray(long *in, size_t size, 
                    long from, long to)
{
  long *pt=in, *fpt;
  fpt=in+size;    
  for(;pt<fpt;pt++)
    if(*pt==from) *pt=to;
}





void
floatinvert(float *in, size_t size, float base)
{
    float *pt=in;
    for(;pt<in+size;pt++)
        *pt=base-*pt;
}




















/*********************************************************************
 **********************      Change type:       **********************
 *********************************************************************/
void
convertftd(float *f, size_t size, double **d)
{
  double *td, *dp;
  float *fp=f, *lfp;

  td=malloc(size*sizeof*td);
  assert(td!=NULL);
 
  lfp=fp+size;  dp=td;
  while(fp<lfp) *dp++=*fp++;

  *d=td;
}




















/*********************************************************************
 **********************      Remove masked      **********************
 *********************************************************************/
void
removemasked(float **in, size_t size, unsigned char *mask, 
             size_t *nsize)
{
  float *tmp;
  size_t i, counter=0;

  tmp=*in;        
  for(i=0;i<size;i++)
    if(mask[i]==0)
      tmp[counter++]=tmp[i];
  *nsize=counter;

  *in=(float *)realloc(*in, counter*sizeof(float));
  assert(*in!=NULL);
}




















/*********************************************************************
 **********************       Make copy         **********************
 *********************************************************************/
void
floatarrcpynomalloc(float *in, size_t size, float **out)
{
  float *ipt, *opt, *fpt;
  ipt=in;
  fpt=*out+size;
  opt=*out;
  while(opt<fpt) 
    *opt++=*ipt++;
}





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





/* Copy the array without copying the masked elements. */
void
floatarrcpymask(float *in, size_t size, 
        unsigned char *mask, size_t *nsize, 
        float **out)
{
  floatarrcpy(in, size, out);

  if(mask==NULL)
    {
      *nsize=size; 
      return;
    }
  else removemasked(out, size, mask, nsize);
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





void
longarrcpy(long *in, size_t size, long **out)
{
  long *ipt, *opt, *fpt;
  ipt=in;
  *out=malloc(size*sizeof(long));
  assert(*out!=NULL);
  fpt=*out+size;
  opt=*out;
  while(opt<fpt) 
    *opt++=*ipt++;
}




















/*********************************************************************
 **********************      Mask an array      **********************
 *********************************************************************/
void
maskfarray(float *in, unsigned char *mask, size_t size, 
        unsigned char f1_b0)
{
  size_t index;
  float min, *fpt=in, *ffpt;
  unsigned char *mpt=mask;

  fmin_withindex(in, size, &min, &index);

  ffpt=fpt+size;
  for(;fpt<ffpt;fpt++,mpt++)
    if(*mpt==f1_b0) *fpt=min;
}





void
masklfarray(float *in, long *mask, size_t size, long f1_b0)
{
  long *mpt=mask;
  float min, *fpt=in, *ffpt;

  floatmin(in, size, &min);

  ffpt=fpt+size;
  if(f1_b0==0)
    {
    for(;fpt<ffpt;fpt++,mpt++)
      if(*mpt==0) *fpt=min;
    }
  else if(f1_b0==1)
    for(;fpt<ffpt;fpt++,mpt++)
      if(*mpt>0) *fpt=min;
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
  while(fpt<lfpt) 
    {
      if(*fpt<min) *(fpt++)=0;
      else fpt++;
    }
}





void
floatsetbelowtomin(float *in, size_t size, float min)
{
  float *fpt, *lfpt;
  fpt=in; 
  lfpt=in+size;
  while(fpt<lfpt) 
    {
      if(*fpt<min) *(fpt++)=min;
      else fpt++;
    }
}





void
floatsetabovetomax(float *in, size_t size, float max)
{
  float *fpt, *lfpt;
  fpt=in; 
  lfpt=in+size;
  while(fpt<lfpt) 
    {
      if(*fpt>max) *(fpt++)=max;
      else fpt++;
    }
}




















/*********************************************************************
 **********************      Shrink Array       **********************
 *********************************************************************/
/* Check to see if a box defined by the two points (x1,y1) and (x2,y2)
   is inside an array of size size1 and size2. If it doesn't overlap,
   then x1=x2 and y1=y2.*/
void
checkifinarray(int *x1, int *y1, int *x2, int *y2, int s0, int s1)
{
  int temp;
  if(*x1==*x2 && *y1==*y2) return;        

  if(*x2<*x1){temp=*x1;*x1=*x2;*x2=temp;} 
  if(*y2<*y1){temp=*y1;*y1=*y2;*y2=temp;}

  if(*x1<0) *x1=0;    if(*x1>s0) *x1=s0;
  if(*y1<0) *y1=0;    if(*y1>s1) *y1=s1;
  if(*x2<0) *x2=0;    if(*x2>s0) *x2=s0;
  if(*y2<0) *y2=0;    if(*y2>s1) *y2=s1;
}





/* We have a large array of size (size1*size2). We want to shrink this
    array, such that (x1,y1) comes down to point (0,0) and the new
    array now only extends to the old (x2,y2). So the size of the new
    array is: w1*w2 where w1=x2-x1 and w2=y2-y1. 

    If the desired region is totally out of the array, a NULL pointer
    is returned.*/
void
floatshrinkarray(float **in, int size1, int size2,
        int x1, int y1, int x2, int y2)
{
  float *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2) 		/* The required region does not */
    {				/* overlap with the array. */
      free(*in);
      *in=NULL;
      return;
    }
  /* The region covers the whole image, no need for the next step. */
  if(x1==0 && y1==0 && x2==size1 && y2==size2) return;

  w1=x2-x1;  w2=y2-y1;
  ux1=x1; uy1=y1; us2=size2;  /* The inputs are int (can be negative,
				 which is allowed: will become zero).
				 but pointers are unsigned, so to 
				 faciliate the process in the loop, 
				 they are converted to size_t. */
  for(i=0;i<w1;i++)
    {
      ofpt=rowstart=*in+i*w2;
      ifpt=*in+(ux1+i)*us2+uy1;
      while(ofpt<rowstart+w2)
	*ofpt++=*ifpt++;
 
    }      
  *in=(float *)realloc(*in, w1*w2*sizeof(float));
  assert(*in!=NULL);    
}





/* similar to floatshrinkarray(), but the old array is kept untouched
   and a new one is created to keep the cropped region. */
void
floatshrinkarraytonew(float *in, int size1, int size2,
		      int x1, int y1, int x2, int y2, 
		      float **out)
{
  float *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2) 	
    {			      
      *out=NULL;
      return;
    }
  if(x1==0 && y1==0 && x2==size1 && y2==size2)
    {
      floatarrcpy(in, size1*size2, out);
      return;
    }

  w1=x2-x1;  w2=y2-y1;
  *out=malloc(w1*w2*sizeof **out);
  assert(*out!=NULL);

  ux1=x1; uy1=y1; us2=size2; 
  for(i=0;i<w1;i++)
    {
      ofpt=rowstart=*out+i*w2;
      ifpt=in+(ux1+i)*us2+uy1;
      while(ofpt<rowstart+w2)
	*ofpt++=*ifpt++; 
    }      
}





void
longshrinkarraytonew(long *in, int size1, int size2,
		      int x1, int y1, int x2, int y2, 
		      long **out)
{
  long *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2) 	
    {			 
      *out=NULL;
      return;
    }
  if(x1==0 && y1==0 && x2==size1 && y2==size2)
    {
      longarrcpy(in, size1*size2, out);
      return;
    }

  w1=x2-x1;  w2=y2-y1;
  *out=malloc(w1*w2*sizeof **out);
  assert(*out!=NULL);

  ux1=x1; uy1=y1; us2=size2;

  for(i=0;i<w1;i++)
    {
      ofpt=rowstart=*out+i*w2;
      ifpt=in+(ux1+i)*us2+uy1;
      while(ofpt<rowstart+w2)
	*ofpt++=*ifpt++;
 
    }      
}




















/*********************************************************************
 **********************      Merge arrays       **********************
 *********************************************************************/
/* a and b are two column arrays with "numrows" rows. It is assumed
   that their first columns are equal. This function will take the
   second column of b and add it to a.  */
#define NUMOUTCOLS 3
void
floatvmerge(float *a, float *b, size_t numrows, float **out)
{
  size_t i;
  float *tmp;
   
  tmp=malloc(NUMOUTCOLS*numrows*sizeof(float));
  assert(tmp!=NULL);
  
  for(i=0;i<numrows;i++)
    {
      tmp[i*NUMOUTCOLS]=a[i*2];
      tmp[i*NUMOUTCOLS+1]=a[i*2+1];
      if(b[i*2]!=a[i*2])
	{
	  printf("\n\n\tError in floatvmerge() (arraymanip.c)\n");
	  printf("\t\tThe first columns are not equal on row %lu\n\n", 
		 i);
	  exit(1);
	}
      else tmp[i*NUMOUTCOLS+2]=b[i*2+1];
    }
  *out=tmp;
}




















/*********************************************************************
 **********************       Print array       **********************
 *********************************************************************/
/* Print a 2D array into an ascii file. Note that s0 and s1 are one larger than the actual number of columns and number of rows.  */
void
printfarray(float *array, size_t s0, size_t s1, 
	    char *comment, char *filename, int decimals)
{
  FILE *fp;
  double *d; 
  size_t i,j;
  char fmt[20];

  sprintf(fmt, "%%10.%df", decimals);
   
  convertftd(array, s0*s1, &d);	/* For printing larger  */
  fp=fopen(filename, "w"); 	/* than 6 decimals. */

  fprintf(fp, "%s", comment);
  for(i=0;i<s0;i++)
    {
      for(j=0;j<s1;j++) 
	{
	  fprintf(fp, fmt, d[i*s1+j]);
	}
      fprintf(fp, "\n");
    }
  free(d);
  fclose(fp);
}
