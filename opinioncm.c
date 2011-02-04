#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

double tau1h(double J, double q, double epsilon, double h) {

  // taxa de transição dependente do campo h, com os paramestros J,
  // q e episilon

  if(h<q) return epsilon;
  else {
    if((q<=h)&&(h<=(1-q))) return (1.0/(1.0+exp(-2*J*(2*h-1))));
    else {
      if(h>(1-q)) return (1.0-epsilon);
      else return 0;
    }
  }

}

double fc(double k, double c, double J, double q, double epsilon) {

  // Map de recorrência para uma rede regular linear de coordenação
  // k, para p modelo descrito por tau1h

  int i;
  double x, sum_x, Cnm;

  sum_x=0.0;
  for(i=0;i<=k;i++) {
    Cnm=gsl_sf_choose(k,i);
    x=Cnm*pow(c,i)*pow(1-c,k-i)*tau1h(J,q,epsilon,i/(double)k);
    sum_x+=x;
  }

  return (double)sum_x;

}

main (int argc, char *argv[]) {

  if(argc<7) {
    printf("Use %s k q J epsilon c0 N\n", argv[0]);
    exit(-1);
  }

  // parâmetros do modelo
  int i, k=atoi(argv[1]), T=atoi(argv[6]);
  double c=atof(argv[5]), J=atof(argv[3]), q=atof(argv[2]);
  double epsilon=atof(argv[4]);

  printf("#k=%d; q=%1.2lf; J=%1.2lf; epsilon=%1.2lf; c=%1.2lf; N=%d\n",k,q,J,epsilon,c,T);

  // Função em solução de campo médio
  printf("%lf\n",c);
  for(i=0;i<T;i++) {
    c=fc(k,c,J,q,epsilon);
    printf("%lf\n",c);
  }

}
