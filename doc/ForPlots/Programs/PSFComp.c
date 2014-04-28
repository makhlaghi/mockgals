#include <math.h>
#include <stdio.h>

/* Functions taken from src/mock.c  */
double
Gaussian(double r, double junk, double a)
{
  junk=1;
  return exp( a*r*r );
}

double
moffat_alpha(double fwhm, double beta)
{
    return (fwhm/2)/pow((pow(2, 1/beta)-1), 0.5f);
}

double
Moffat(double rda, double nb, double junk)
{
  junk=1; 
  return pow(junk+rda*rda, nb);
}

int
main(void)
{
  FILE *out;
  double co_g, co_m;
  float r, fwhm=5, sigma;
  double p1_g, p1_m2, p1_m3, p1_m4, p1_m5;
  double p2_g, p2_m2, p2_m3, p2_m4, p2_m5;

  out=fopen("MoffatGaussianComp.txt", "w");

  p1_g=1; p2_g=0;
  sigma=fwhm/2.35482;
  co_g=-1.0f/(2.0f*sigma*sigma);

  p1_m2=moffat_alpha(fwhm, 2);    p2_m2=-2.0f;
  p1_m3=moffat_alpha(fwhm, 3);    p2_m3=-3.0f;
  p1_m4=moffat_alpha(fwhm, 4);    p2_m4=-4.0f;
  p1_m5=moffat_alpha(fwhm, 5);    p2_m5=-5.0f;
  co_m=0;

  for(r=0.0f; r<25; r+=0.1)
    fprintf(out, "%-10.1f %-10f %-10f %-10f %-10f %-10f\n", r,
	    Gaussian(r/p1_g, p2_g, co_g),
	    Moffat(r/p1_m2, p2_m2, co_m),
	    Moffat(r/p1_m3, p2_m3, co_m),
	    Moffat(r/p1_m4, p2_m4, co_m),
	    Moffat(r/p1_m5, p2_m5, co_m));

  fclose(out);
}

