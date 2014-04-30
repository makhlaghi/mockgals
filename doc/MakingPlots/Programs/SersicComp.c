#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Functions taken from src/mock.c  */
double
sersic_b(double n)
{
  if(n<=0.35f) 
    {
      printf("\n\n   ##ERROR in mock.c's sersic_b()\n");
      printf("\t Sersic n (=%f) must be smaller than 0.35\n\n", n);
      exit(EXIT_FAILURE);
    }
  return 2*n-(1/3)+(4/(405*n))+(46/(25515*n*n))+
    (131/(1146175*n*n*n)-(2194697/(30690717750*n*n*n*n)));
}
double
Sersic(double rdre, double inv_n, double nb)
{
  return exp( nb*( pow(rdre,inv_n)-1 ) );
}


int
main(void)
{
  FILE *out;
  double r, p1=5;
  double p2_05, p2_1, p2_25, p2_4, p2_6;
  double co_05, co_1, co_25, co_4, co_6;
  double s0_05, s0_1, s0_25, s0_4, s0_6;
	
  out=fopen("SersicComp.txt", "w");

  p2_05=1.0f/0.5f;   co_05=-1*sersic_b(0.5f);
  p2_1=1;            co_1=-1*sersic_b(1);
  p2_25=1/2.5f;      co_25=-1*sersic_b(2.5);
  p2_4=1.0f/4.0f;    co_4=-1*sersic_b(4);
  p2_6=1.0f/6.0f;    co_6=-1*sersic_b(6);

  s0_05=Sersic(0, 0.5, co_05);
  s0_1=Sersic(0, 1, co_1);
  s0_25=Sersic(0, 2.5, co_25);
  s0_4=Sersic(0, 4, co_4);
  s0_6=Sersic(0, 6, co_6);

  for(r=0.0f; r<25.2; r+=0.1)
    fprintf(out, "%-15.1f %-15.8f %-15.8f %-15.8f %-15.8f %-15.8f\n", r,
            Sersic(r/p1, p2_05, co_05)/s0_05,
            Sersic(r/p1, p2_1, co_1)/s0_1,
            Sersic(r/p1, p2_25, co_25)/s0_25,
            Sersic(r/p1, p2_4, co_4)/s0_4,
            Sersic(r/p1, p2_6, co_6)/s0_6);

  fclose(out);
}

