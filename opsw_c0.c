#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define I (sizeof(int))
#define Ip (sizeof(int *))
#define D (sizeof(double))

typedef struct parametros {
  
  // Modelo de opinião
  double eps, J, q;
  int T_mcs, D_mcs;

}par;

// Variaveis globais
int n, K;
double p;

void makesw(int **rede, unsigned long int sem) {

  // Constrói a rede small-world com base numa lista de vizinhança, em que o primeiro elemento é a conectividade do sítio

  int i, j, nex, nrand;
  double ran;
  gsl_rng *r;

  // iniciar semente
  r=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,sem);
  for(i=0;i<n;i++) {
    for(j=1;j<=K;j++) {
      ran=gsl_rng_uniform(r);
      if(ran<p) {
	nrand=(i+K+1+gsl_rng_get(r)%(n-K-1))%n;
	// Aloca para poder ter espaço, i-->nrand
	rede[i][0]++; // Incrementa a conectividade
	rede[i]=(int *)realloc(rede[i],(rede[i][0]+2)*I);
	rede[i][rede[i][0]]=nrand;
	// Aloca para poder ter espaço, nrand-->i
	rede[nrand][0]++; // Incrementa a conectividade
	rede[nrand]=(int *)realloc(rede[nrand],(rede[nrand][0]+2)*I);
	rede[nrand][rede[nrand][0]]=i;
	
      }
      else {
	nex=(i+j)%n;
	// Aloca para poder ter espaço, i-->nex
	rede[i][0]++;
	rede[i]=(int *)realloc(rede[i],(rede[i][0]+2)*I);
	rede[i][rede[i][0]]=nex;
	// Aloca para poder ter espaço, nex-->i
	rede[nex][0]++;
	rede[nex]=(int *)realloc(rede[nex],(rede[nex][0]+2)*I);
	rede[nex][rede[nex][0]]=i;
      }
    }
  }
 
  gsl_rng_free(r);

}

double *make_tau(int n, double J, double q, double epsilon) {
  // --------------------------------------
  // Define as taxas de transica do modelo.
  // --------------------------------------

  int s, k;
  double *tau, h;

  tau=(double *)calloc((n+1)*(n+2)/2,D);
  if(tau==NULL) {
    puts("Nao foi possivel alocar as taxas");
    exit(-1);
  }

  for(k=1;k<=n;k++)
    for(s=0;s<=k;s++) {
      h=(double)s/(double)k;
      if(h<q)                                                                     // 1a condição
	tau[s-1+(k+1)*k/2]=epsilon;
      else {
	if(h>(1-q)) tau[s-1+(k+1)*k/2]=1.0-epsilon;                               // 2a condição
	else {
	  if((q<=h)&&(h<=(1-q))) tau[s-1+(k+1)*k/2]=1.0/(1.0+exp(-2*J*(2*h-1)));  // 3a condição
	  else {
	    printf("Erro, condição inválida encontrada.\n");
	    exit(-1);
	  }
	}
      }
    }

  return(tau);

}

int zero(int sem, double rho_0) {

  // Constrói um estado incial com uma concentração rho_0 de sítios ativos

  int i, *s_0;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_mt19937);;

  // Iniciando a semente de números aleatórios
  gsl_rng_set(r,sem_rede);

  // Alocando o vetor de estado
  s_0=(int *)calloc(n,I);

  for(i=0;i<n;i++) {
    if(gsl_rng_uniform(r)<rho_0) s_0[i]=1;
  }

  return s_0;

}

void opiniao(int **rede, par *P, int s[], unsigned long int sem, unsigned long int sem_rede, char out[], char fran[]) {
  // ---------------------------------------------------
  // Aqui é a parte da simulação da dinâmica de opinião.
  // ---------------------------------------------------

  // Obtendo os parâmetros
  int *h, k_max=0; // vetor de conectividade e com o estado dos sitios;
  double epsilon=P->eps, J=P->J, q=P->q;
  int i, w, t, tt, T=P->T_mcs, nt=P->D_mcs; // passo monte-carlo e iterador de sítios da rede
  double *tau, rho=0; // campo local de cada sítio e taca de transição
  FILE *f=fopen(out,"w"), *g, *ff=fopen(fran,"r");
  char ou[255];
  // Variáveis do gerador
  gsl_rng *r;

  // Arquivo de saída de estado
  strcpy(ou,out);
  strcat(ou,"~");

  // Iniciando a galera
  h=(int *)calloc(n,I);
  // Pegar o grau máximo
  for(i=0;i<n;i++)
    if(rede[i][0]>k_max) k_max=rede[i][0];
  
  // Constrói as taxas de k=0 até o grau máximo
  tau=make_tau(k_max,J,q,epsilon);
  // Iniciando a semente de números aleatórios
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_fread(ff,r);
  fclose(ff);

  for(i=0;i<n;i++)
    rho+=s[i];
  rho/=((double)n);
  // Primeiro ponto na saida
  fprintf(f, "%lf\n", rho);

  // O monte-carlo começa aqui
  for(t=0;t<T;t+=nt) {
    g=fopen(ou,"w");
    ff=fopen(fran,"w");
    fprintf(g,"# N=%d\tk=%d\tp=%lf\tsem=-1\tJ=%lf\teps=%lf\tq=%lf\tSR=%lu\n",n,K,p,J,epsilon,q,sem_rede);
    fprintf(g, "# Semente: -1\n");
    for(tt=0;tt<nt;tt++) {
            
      // Calcular o vetor h dos sítios
      for(i=0;i<n;i++) {
	h[i]=0;
	for(w=1;w<=rede[i][0];w++) {
	  h[i]+=s[rede[i][w]];
	}
      }

      // Atualiza os sítios
      for(i=0;i<n;i++) {
	if(gsl_rng_uniform(r)<tau[h[i]-1+(rede[i][0]+1)*rede[i][0]/2]) {
	  if(s[i]==0) s[i]=1;
	}
	else
	  if(s[i]==1) s[i]=0;
      }

      // Calcula a densidade de estados
      rho=0;
      for(i=0;i<n;i++)
	rho+=s[i];
      rho/=((double)n);
      
      // Imprime saída
      fprintf(f, "%lf\n", rho);
      
    }
    // Salvando o estado
    for(i=0;i<n;i++)
      fprintf(g, "%d\n", s[i]);
    gsl_rng_fwrite(ff,r);
    fclose(g);
    fclose(ff);
  }
  
  free(tau);
  free(h);
  // Libera memória do gerador
  gsl_rng_free (r);
  
}

int main(int argc, char *argv[]) {

  if(argc<6) {
    printf("Use %s <1. estado inicial> <2. T mcs> <3. T mcs para salvar> ",argv[0]);
    printf("<4. TS para densidade> <5. bin do gerador>\n");
    exit(-1);
  }

  unsigned long int sem, sem_rede;
  int i, **rede, *s;
  par P;
  char out[255], fran[255];
  FILE *in=fopen(argv[1],"r");

  // int *hist;

  if(in==NULL) {
    printf("Arquivo de estado inicial inválido.\n");
    exit(-1);
  }

  // Coletando informações a partir do arquivo de estado inicial
  fscanf(in,"# N=%d\tk=%d\tp=%lf\tsem=%lu\tJ=%lf\teps=%lf\tq=%lf\tSR=%lu\n",&n,&K,&p,&sem,&P.J,&P.eps,&P.q,&sem_rede);
  fscanf(in, "# Semente: %lu\n", &sem);
  P.T_mcs=atoi(argv[2]); P.D_mcs=atoi(argv[3]);

  // Obtendo estado inicial
  s=(int *)calloc(n,I);
  for(i=0;i<n;i++) {
    fscanf(in,"%d",&s[i]);
  }
  fclose(in);

  // Nomes para saída do modelo e entrada do gerador
  strcpy(out,argv[4]);
  strcpy(fran,argv[5]);

  // Mostra na tela os parâmetros do problema
  printf("\n# N=%d\tk=%d\tp=%lf\tsem=%lu\tJ=%lf\teps=%lf\tq=%lf\tSR=%lu\n",
  	  n,K,p,sem,P.J,P.eps,P.q,sem_rede);

  // Alocando a estrutura de armazenamento da rede
  rede=(int **)malloc(n*Ip);
  for(i=0;i<n;i++) {
    rede[i]=(int *)calloc(1,I);
  }
  
  makesw(rede,sem_rede);
  opiniao(rede,&P,s,sem,sem_rede,out,fran);


  /* hist=(int *)calloc(n,I); */
  /* for(i=0;i<n;i++) { */
  /*   hist[rede[i][0]]++; */
  /* } */

  /* for(i=0;i<n;i++) */
  /*   if(hist[i]!=0) printf("%d %d\n",i,hist[i]); */

  /* for(i=0;i<n;i++) { */
  /*   free(rede[i]); */
  /* } */
  /* free(hist); */

  free(rede);

  //printf("A semente foi %lu.\n", sem);

  return 0;

}
