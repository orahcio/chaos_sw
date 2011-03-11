#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

// Nesta versão podemos escolher a densdade de sítios ativos

int main(int argc, char *argv[]) {

  if(argc<12) {
    //                                9                      10                                 11
    printf("Use %s <dados de entrada> <saida estado inicial> <saida pra semente do gerador GSL> <densidade de sitios ativos>\n",argv[0]);
    printf("Entrada: <n sítios>"); // 1
    printf(" <k vizinhos pra frente>"); // 2
    printf(" <p WS-model>"); // 3
    printf(" <semente>"); // 4
    printf(" <semente da rede>"); // 5
    printf(" <epsilon>"); // 6
    printf(" <J>"); // 7
    printf(" <q>"); // 8
    printf("Veja o cabeçalho do fonte para as definições dos parametros\n");
    exit(-1);
  }

  int i, n=atoi(argv[1]), k=atoi(argv[2]), T=atoi(argv[4]), sem=atoi(argv[4]), sem_rede=atoi(argv[5]), N=atoi(argv[10]);
  double p=atof(argv[3]),eps=atof(argv[6]), J=atof(argv[7]), q=atof(argv[8]), dens=atof(argv[11]);
  FILE *out=fopen(argv[9],"w"), *g=fopen(argv[10],"w");
  gsl_rng *r;

  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r,sem_rede);

  // Cabeçalho da saída
  fprintf(out, "# N=%d\tk=%d\tp=%lf\tsem=%u\tJ=%lf\teps=%lf\tq=%lf\tSR=%u\n",
	  n,k,p,sem,J,eps,q,sem_rede);

  fprintf(out, "# Semente: %u\n", sem);

  for(i=0;i<n;i++) {
    if(gsl_rng_uniform(r)<dens) fprintf(out,"1\n");
    else fprintf(out,"0\n");
  }


  gsl_rng_set(r,sem);
  gsl_rng_fwrite(g,r);

  gsl_rng_free(r);
  fclose(out);

  return 0;

}
