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
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_sort.h>	 /* gsl_sort2_float */

#include "pix.h"
#include "sll.h"
#include "mock.h"
#include "stats.h"
#include "attaavv.h"
#include "raddist.h"
#include "convolve.h"
#include "profiles.h"
#include "integtwod.h"
#include "arraymanip.h"
#include "fitsarrayvv.h"
#include "macrofunctions.h"

#include "rand.h"		/* Needs integtwod! */
#include "ui.h"			/* Needs mock.h */





/****************************************************************
 *****************  Fill in the mock values  ********************
 ****************************************************************/
/* Until calling this function, findstartingpixel(), looks on a
   cartesian space. Therefore it will be give a wrong result (by one
   pixel) when the point is near the boundaries and the profile is
   extremely elliptical. So this function will search the */
void
fsp_checkell(size_t *ngbs, size_t s1, struct elraddistp *e, 
	     size_t *p)
{
  float r, min_r;
  size_t n, nrow, numngb, min_p;

  min_r=elraddist(e, *p/s1, *p%s1);
  min_p=*p;

  nrow=*p*NGBSCOLS;
  numngb=1;
  do
    {
      n=ngbs[nrow+numngb];
      if( (r=elraddist(e, n/s1, n%s1))<min_r )
	{
	  min_r=r;
	  min_p=n;
	}
    }
  while(ngbs[nrow+ ++numngb]!=NONINDEX);

  *p=min_p;
}


/* Find the first pixel in the image to begin building the profile.

   The input sizes and positions are based on the FITS standard,
   But in main(), we reversed the sizes to fits the C standard and
   when calling this function, we reversed the positions to fit
   the C standard. So by the time they get here, the inputs are
   all in the C standard.*/
void
findstartingpixel(size_t *ngbs, size_t s0, size_t s1, float truncr, 
		  struct elraddistp *e, size_t *p)
{
  float rmin, x_c, y_c, fs0, fs1;
  int is0, is1, i, j, x, y, x1, y1, x2, y2;
  size_t x_w, y_w, xmin=NONINDEX, ymin=NONINDEX;

  is0=s0; is1=s1;
  x_c=e->xc; y_c=e->yc;

  /* Find the central pixel, this will be needed if it is inside the
     image or outside it. */
  if(x_c-(int)x_c>0.5) x=x_c+1;
  else x=x_c;
  if(y_c-(int)y_c>0.5) y=y_c+1;
  else y=y_c;

  /* The central pixel is in the image, set the pixel and return. */
  fs0=s0;			/* Just to make things easier! */
  fs1=s1;
  if(x_c>=0 && x_c<fs0 && y_c>=0 && y_c<fs1)
    {
      *p=x*s1+y;
      fsp_checkell(ngbs, s1, e, p);
      return;
    }

  /* The center is out of the image. Use encloseellipse() from
   raddist.c to see if any of the pixels within the truncation
   radius fit into the image.*/
  encloseellipse(truncr, e->q*truncr, e->t, &x_w, &y_w);
  x1=x-x_w/2;
  y1=y-y_w/2;
  x2=x+x_w/2;
  y2=y+y_w/2;
  
  /* Check if any of the four corners of the box inclosing the
     profile are in the mock image. If they are not, x1=x2 or
     y1=y2*/
  checkifinarray(&x1, &y1, &x2, &y2, s0, s1);
  if(x1==x2 || y1==y2)	    
    {			     /* The profile's region is */
      *p=NONINDEX;	     /* Completely out of the image. */
      return;		     /* Return NULL. */
    }
  else			     /* The profile and the image overlap */
    {			     /* Find the point on the side of the */
      rmin=1e10;	     /* image with the smallest radius. */
      if(x1==0)		     /* This is important, because the  */
	for(j=y1;j<y2;j++)   /* first check later will be  */
	  if(elraddist(e, 0, j)<rmin)
	    {		     /* integration and we want to be sure  */
              xmin=0; ymin=j; /* we start with the smallest radius. */
	      rmin=elraddist(e, 0, j); /* in the image. */
	    }
      if(x2==is0)		       
	for(j=y1;j<y2;j++)
	  if(elraddist(e, s0-1, j)<rmin)
	    {
              xmin=s0-1; ymin=j;
	      rmin=elraddist(e, s0-1, j);
	    }
      if(y1==0)		       
	for(i=x1;i<x2;i++)
	  if(elraddist(e, i, 0)<rmin)
	    {
              xmin=i; ymin=0;
	      rmin=elraddist(e, i, 0);
	    }
      if(y2==is1)		       
	for(i=x1;i<x2;i++)
	  if(elraddist(e, i, s1-1)<rmin)
	    {
              xmin=i; ymin=s1-1;
	      rmin=elraddist(e, i, s1-1);
	    }
      if(rmin<truncr && xmin!=NONINDEX && ymin!=NONINDEX)
	*p=xmin*s1+ymin;
      else *p=NONINDEX;
    }

  fsp_checkell(ngbs, s1, e, p);
}





/* Make a profile in an array.

   The logic: 'byt' is an array the same size as the image, that 
   should be completely zero upon making each profile. If it is
   zero at first, it will be zero once this function is finished
   with it. It is used to mark which pixels have been checked in 
   the image. 'bytind[]' will keep the index (1D) of those pixels that
   have been marked so after the job is finished, it can set them
   all back to zero and not bother with resetting the whole array!

   We begin with the pixel in the image that is closest to the desired
   center of the profile. A pixel queue (a FIFO or simple linked list)
   will keep all the neighbors of all the pixels in order to check
   them all. 

   Since some profiles will fall on the sides of the image, there is
   no way we can calculate the whole flux by summing of the pixels and
   then setting them to the desired value. So we have to use
   integration, the profiles were all made with their constant set to
   1. We integrate over the profile to infinity over a surface and
   consider that as the total flux. Note that when the truncation
   radius is small, this theoretical total flux will be much larger
   (>10%) than the actual total flux of what is actually put in the
   image if it all fits in.  Since that total flux will be set to the
   desired value, then if the truncation radius is too small, the real
   total flux in the image is slightly lower than the desired
   value. This problem did not exist in the makeprofile_old() function
   which actually allocated an array for each profile, summed over it
   to find the total flux and then only used the intersection with the
   main image to put the profile's image into the main image. But that
   was too slow, especially if a large number of profiles were
   needed. It was removed on April 12th, 2014. So if you are
   interested, you can see it in the commits before this date.

   For all pixels, first it is checked if the pixel is within the
   truncation radius. If it isn't, its neighbors will not be added to
   the queue, but it's byt value will be set to one and a value of
   zero will be put into its D[i]->v value.
 */
int
makeprofile(float *img, unsigned char *byt, size_t *bytind, 
	    size_t *ngbs, size_t s0, size_t s1, float trunc, 
	    float integaccu, float s0_m1_g2_p3, float x_c, float y_c, 
	    double p1, double p2, float pa_d, float q, float avflux, 
	    double *totflux)
{
  float t_i, t_j;		/* 2D position from 1D. */
  char profletter;
  struct elraddistp e;
  struct ssll *Q=NULL;
  struct integparams ip;
  float accurate, approx;
  double sum=0, area=0, co;
  float r, truncr, multiple=0;
  struct tossll *lQ=NULL, *sQ;	/* lQ: Largest. sQ: Smallest in queue */
  double (*func)(double, double, double);
  size_t i, numngb, nrow, counter=0, p, tp;
  int userandpoints=1, outofaccurateloop=0;
  
  e.q=q;
  e.t=M_PI*pa_d/180;
  e.xc=x_c;          e.yc=y_c;
  e.cos=cos(e.t);    e.sin=sin(e.t);

  setintegparams(s0_m1_g2_p3, p1, p2, pa_d, q, trunc, 
		 &truncr, &profletter, &ip);

  findstartingpixel(ngbs, s0, s1, truncr, &e, &p);
  if(p==NONINDEX)
    return 0;	      /* Profile is completely out of image. */

  area=M_PI*truncr*truncr*q;
  if(s0_m1_g2_p3==0)
    sum=totsersic(p2, p1, sersic_b(p2), q);
  else if(s0_m1_g2_p3==1)
    sum=totmoffat(ip.p1, p2, q); /* Note that p->p2=-1/p2!. In Moffat*/
  else if(s0_m1_g2_p3==2)	 /* we want p2 (beta), not p->p2! */
    sum=totgaussian(q);
  else if(s0_m1_g2_p3==3)
    {
      img[p]=avflux;
      *totflux=avflux;
      return 1;			/* Successful. */
    }
  else
    {
      printf("\n\ns0_m1_g2_p3=%.0f is not recognized.\n\n", 
	     s0_m1_g2_p3);
      exit(EXIT_FAILURE);
    }
  *totflux=avflux*area;
  multiple=*totflux/sum;

  add_to_tossll_end( &lQ, &sQ, p, elraddist(&e, p/s1, p%s1) );
  byt[p]=1;
  bytind[counter++]=p;

  co=ip.co;
  func=ip.profile;
  p1=ip.p1;
  p2=ip.p2;

  while(sQ)
    {
      /* In case you want to see the status of the twosided ordered
	 queue, increasing and decreasing side by side, uncomment this
	 line. Note that there will be a lot of lines printed! */
      /*print_tossll(lQ, sQ);*/

      pop_from_tossll_start(&lQ, &sQ, &p, &r);

      if(r>truncr) continue;
      
      /* Find the value for this pixel: */
      t_i=p/s1-x_c;        t_j=p%s1-y_c;
      ip.xl=t_i-0.5;       ip.xh=t_i+0.5;
      ip.yl=t_j-0.5;       ip.yh=t_j+0.5;

      if(userandpoints)	       /* See when 10e5 randomly chosen */
	{		       /* points are no longer needed. */
	  accurate=randompoints(&ip);	
	  approx=integ2d(&ip);		
	  img[p]+=accurate*multiple;
	  /*printf("%-5d", userandpoints);*/
	  if (fabs(accurate-approx)/accurate<integaccu) 
	    userandpoints=0;
	}
      else			 /* See when integration is not */
	{			 /* needed and you can switch */
	  accurate=integ2d(&ip); /* to the pixel center. */
	  approx=func(r/p1, p2, co); 
	  img[p]+=accurate*multiple;
	  /*printf("%-5d", userandpoints);*/
	
	  if (fabs(accurate-approx)/accurate<integaccu) 
	    outofaccurateloop=1;
	 
	}

      /*array_to_fits("tmp.fits", NULL, "", FLOAT_IMG, img, s0, s1);*/
      /*printf(" (%lu, %lu): %-10.4f %-12.8f %-10.8f \n", 
	     p%s1+1, p/s1+1, r, fabs(accurate-approx)/accurate, 
	     log10(accurate/sum));*/
      

      if(outofaccurateloop) break;
      
      /* It is very important to go over the 8 connected neighbors,
	 because we are dealing with elliptical radii, not a cartesian
	 space. So the diagonal neighbor might be closer (in
	 elliptical radii) then the 4 connected pixels. */
      nrow=p*NGBSCOLS;
      numngb=1;	   /* all 8-connected neighbors will be checked.*/
      do
	if(byt[ tp=ngbs[nrow+numngb] ]==0)
	  {
	    byt[tp]=1;
	    bytind[counter++]=tp;
	    add_to_tossll_end( &lQ, &sQ, tp, 
			       elraddist(&e, tp/s1, tp%s1) );
	  }
      while(ngbs[nrow+ ++numngb]!=NONINDEX);
    }

  /* All the pixels that required integration are now done, so we
     don't need an ordered array any more! */
  tossll_into_ssll(lQ, &Q);

  while(Q)
    {
      pop_from_ssll(&Q, &p);

      r=elraddist(&e, p/s1, p%s1);
      if(r>truncr) continue;
      
      /* Find the value for this pixel: */
      img[p]+=func(r/p1, p2, co)*multiple; 

      /*array_to_fits("tmp2.fits", NULL, "", FLOAT_IMG, img, s0, s1);*/

      /* Here it is best to use 4 connectivity, because the 8
      connected ones have most probably been done before. */
      nrow=p*NGBSCOLS;
      numngb=ngbs[nrow];	/* Check only 4 connected neighbors */
      do
	if(byt[ tp=ngbs[nrow+numngb] ]==0)
	  {
	    byt[tp]=1;
	    bytind[counter++]=tp;
	    add_to_ssll(&Q, tp);
	  }
      while(ngbs[nrow+ ++numngb]!=NONINDEX);
    }

  /* Clean the byt image for the next profile. */
  for(i=0;i<counter;i++)
    byt[bytind[i]]=0;
  return 1;
}




















/****************************************************************
 *****************    Read or make the PSF   ********************
 ****************************************************************/
void
makepsf(float **psf, size_t *psf_s0, size_t *psf_s1, 
	      struct mockparams *p)
{
  int suc;
  double t;
  float pa_d=0, q=1;
  unsigned char *byt;
  float *tpsf, trunc_r;
  size_t w, *bytind, *ngbs;

  trunc_r=p->psf_t * p->psf_p1/2;
  w=2*trunc_r+1;
  if(w%2==0) w--;		/* To make sure the width is odd. */

  assert( (tpsf=calloc(w*w, sizeof *tpsf))!=NULL );
  assert( (byt=calloc(w*w, sizeof *byt))!=NULL );
  assert( (bytind=calloc(w*w, sizeof *bytind))!=NULL );

  imgngbs(w, w, &ngbs);

  /* We are going to fix the total flux here, so the average flux is
     junk! So is t! The total flux comes from integration, which is
     not what we want. The truncation radius might be too small. */
  suc=makeprofile(tpsf, byt, bytind, ngbs, w, w, p->psf_t, 
		  p->integaccu, p->psf_mg, w/2, w/2, p->psf_p1,
		  p->psf_p2, pa_d, q, 1, &t);
  
  floatarrmwith(tpsf, w*w, 1/floatsum(tpsf, w*w));

  free(byt);
  free(bytind);

  *psf=tpsf;
  *psf_s0=*psf_s1=w;
}





void
readormakepsf(float **psf, size_t *psf_s0, size_t *psf_s1, 
	      struct mockparams *p)
{
  void *tmp;
  int  bitpix;
  float psfsum;
  FILE *tmpfile;
  char fitspsfname[]="PSF.fits";
  char asciipsfname[]="PSF.conv";

  if(strlen(p->psfname))	/* A PSF file name was given. */
    {
      if ((tmpfile = fopen(p->psfname, "r")) != NULL) 
	{			/* The file exists! */
	  fclose(tmpfile);
	  fits_to_array(p->psfname, 0, &bitpix, &tmp, psf_s0, psf_s1);

          /* Check if the sides are odd */
	  if(*psf_s0%2==0)
	    {
	      printf("\n\n\tError: NAXIS2 of %s (PSF) must be odd!\n\n",
		     p->psfname);
	      exit(EXIT_FAILURE);
	    }
	  if(*psf_s1%2==0)
	    {
	      printf("\n\n\tError: NAXIS1 of %s (PSF) must be odd!\n\n",
		     p->psfname);
	      exit(EXIT_FAILURE);
	    }

	  /* Everything is fine, make sure the output is in float. */
	  convertanytofloat(tmp, *psf_s0 * *psf_s1, bitpix, psf,
			    BYTE_IMG, SHORT_IMG, LONG_IMG, FLOAT_IMG, 
			    DOUBLE_IMG);

	  /* Check if the sum of the PSF is unity. */
	  psfsum=floatsum(*psf, *psf_s0 * *psf_s1);
          if(psfsum!=1)
	    floatarrmwith(*psf, *psf_s0 * *psf_s1, 1/psfsum);
	}
      else
	{
	  printf("\n\n\tError: %s could not be opened\n\n",p->psfname);
	  exit(EXIT_FAILURE);
	}
    }
  else
    {
      makepsf(psf, psf_s0, psf_s1, p);
      if(p->vpsf || p->ovpsf)
	{
	  /* Save the PSF in a FITS image: */
	  if ((tmpfile = fopen(fitspsfname, "r")) != NULL) 
	    {			/* The file exists! */
	      fclose(tmpfile);
	      assert(remove(fitspsfname)==0);
	      printf("\nWARNING:  %s already existed, was deleted\n",
		     fitspsfname);
	    }
	  array_to_fits(fitspsfname, NULL, "PSF", FLOAT_IMG, 
			*psf, *psf_s0, *psf_s1);

	  /* Save the PSF as an ASCII .conv file */
	  if ((tmpfile = fopen(asciipsfname, "r")) != NULL) 
	    {		       
	      fclose(tmpfile);
	      assert(remove(asciipsfname)==0);
	      printf("\nWARNING:  %s already existed, was deleted\n\n",
		     asciipsfname);
	    }
	  printfarray(*psf, *psf_s0, *psf_s1, "CONV NONORM\n", 
		      asciipsfname, 4, 'e');
	}
    }
}

















/****************************************************************
 *****************     Main output program   ********************
 ****************************************************************/
/* Put one or more mock profiles into and image, convolve it and add
   noise to the final result.  The convolution is going to make the
   sides darker.  So the actual image where the galaxies will be
   placed is going to be larger than the desired image.  After
   convolution those sides will be trimed and the pixels on the sides
   of the result will be equally convolved as the central pixels.

   Note on MINFLOAT: We are using float values here, due to roundoff
   errors, after convolution we might have pixel values less than
   NUMMOCK, which are pure error, so they will all be set to zero. */
#define MINFLOAT 1e-8f
void
mockimg(struct mockparams *p)
{
  double *pp, ss;
  unsigned char *byt;
  int extcounter=0, suc;
  float *conv, *psf, *img;
  float *nonoisehist, *preconv, trunc=10;
  size_t i, psf_s0, psf_s1, ns0, ns1, *ngbs;
  size_t nc, nsize, size, hs0, hs1, *bytind;

  readormakepsf(&psf, &psf_s0, &psf_s1, p);
  if(p->ovpsf) {free(psf); return;}

  hs0=psf_s0/2;          hs1=psf_s1/2;
  ns0=p->s0+2*hs0;       ns1=p->s1+2*hs1;
  size=p->s0*p->s1;	 /* Shorter name ;-). */
  nsize=ns0*ns1;	 /* Shorter name ;-). */
  nc=p->numppcols;	 /* Shorter name ;-). */
  pp=p->profileparams;	 /* Shorter name ;-). */
  ss=sqrt(p->sky);	 /* Shorter name ;-). */

  assert( (img=calloc(nsize,sizeof *img))!=NULL );
  assert( (byt=calloc(nsize,sizeof *byt))!=NULL );
  assert( (bytind=malloc(nsize*sizeof *bytind))!=NULL );

  imgngbs(ns0, ns1, &ngbs);

  if(p->verb)
    {
      printf("(x,y): profile position.\n");
      printf("\t Y : At least part of it was in the image.\n");
      printf("\t-*-: Profile was not in the image.\n");
    }

  for(i=0;i<p->nummock;i++)
    {
      suc=makeprofile(img, byt, bytind, ngbs, ns0, 
		      ns1, trunc, p->integaccu, 
		      pp[i*nc+1],	  /* Profile function. */
		      pp[i*nc+3]+hs0-1,   /* x_c (C format) */
		      pp[i*nc+2]+hs1-1,   /* y_c (C format) */
		      pp[i*nc+4],	  /* p1: sersic re. */
		      pp[i*nc+5],	  /* p2: sersic n. */
		      90-pp[i*nc+6],	  /* position angle. */
		      pp[i*nc+7],	  /* axis ratio. */
		      ss*pp[i*nc+8],	  /* average flux.*/
		      &pp[i*nc+9]);	  /* Total flux of profile*/ 
      if(p->verb)
	printf(" -(%-.2f,%-.2f),n=%.2f, re=%.2f, pa=%.2f, q=%.2f\t%s\n",
	       pp[i*nc+2], pp[i*nc+3], pp[i*nc+5], pp[i*nc+4], 
	       pp[i*nc+6], pp[i*nc+7], suc ? " Y" : "-*-");
    }
  if(p->verb)
    printf("\n\n");

  /* The user wants to see the unconvoved image. Since the image to
     convolve is larger than the final image, first crop out the
     sides, then save that central part. */
  if(p->vnoconv)
    {
      floatshrinkarraytonew(img, ns0, ns1, hs0, hs1, 
			    p->s0+hs0, p->s1+hs1, &preconv);
      array_to_fits(p->outname, NULL, "NOCONV", FLOAT_IMG, 
		    preconv, p->s0, p->s1);
      free(preconv);
      if(p->verb)
	printf("- Pre-convolved profiles saved in '%s' (ext %d)\n",
	       p->outname, extcounter++);
    }

  convolve(img, ns0, ns1, psf, psf_s0, psf_s1, &conv);

  floatshrinkarray(&conv, ns0, ns1, hs0, hs1, p->s0+hs0, p->s1+hs1);

  floatsetbelowtozero(conv, size, MINFLOAT);

  if(p->vconv)
    {
      array_to_fits(p->outname, NULL, "NONOISE", FLOAT_IMG, conv, 
		    p->s0, p->s1);
      if(p->verb)
	printf("- Convolved profiles saved in '%s' (ext %d)\n",
	       p->outname, extcounter++);
    }
   
  if(p->vhist)
    histogram(conv, size, p->vhist, &p->histmin, 
	      &p->histmax, &nonoisehist, 1, 0, 0);

  addnoise(conv, size, p->sky);

  array_to_fits(p->outname, NULL, "WITHNOISE", FLOAT_IMG, conv, 
		p->s0, p->s1);

  if(p->verb)
    printf("- Noised image saved in '%s' (ext %d)\n",
	   p->outname, extcounter++);
   
  if(p->vhist)
    printmockhist(conv, size, p->vhist, p->histmin, 
		  p->histmax, nonoisehist);

  savemockinfo(p);
  if(p->verb)
    printf("- Profile info saved in '%s'\n\n", p->infoname);

  free(psf);
  free(img);
  free(byt);
  free(conv);
  free(ngbs);
  free(bytind);
}
