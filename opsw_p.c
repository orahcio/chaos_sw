#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define I (sizeof(int))
#define Ip (sizeof(int *))
#define D (sizeof(double))
#define IDX(h,k) ((h)-1+((k)+1)*(k)/2)

typedef struct parametros {
  
  // Modelo de opinião
  double eps, J, q;
  int T_mcs, D_mcs;

}par;

// Variaveis globais
int n, K;
double p;

void makesw(int **rede, gsl_rng *r) {

  // Constrói a rede small-world com base numa lista de vizinhança, em que o primeiro elemento é a conectividade do sítio

  int i, j, nex, nrand;
  double ran;
  
  // Conectividades nulas inicialmente
  for(i=0;i<n;i++) rede[i][0]=0;

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
 
}

double *make_tau(int kmax, double J, double q, double epsilon) {
  // --------------------------------------
  // Define as taxas de transica do modelo.
  // --------------------------------------

  int s, k;
  double *tau, h;

  tau=(double *)calloc((kmax+1)*(kmax+2)/2,D);
  if(tau==NULL) {
    puts("Nao foi possivel alocar as taxas");
    exit(-1);
  }

  for(k=1;k<=kmax;k++)
    for(s=0;s<=k;s++) {
      h=(double)s/(double)k;
      if(h<q)                                                                     // 1a condição
	tau[IDX(s,k)]=epsilon;
      else {
	if(h>(1-q)) tau[IDX(s,k)]=1.0-epsilon;                               // 2a condição
	else {
	  if((q<=h)&&(h<=(1-q))) tau[IDX(s,k)]=1.0/(1.0+exp(-2*J*(2*h-1)));  // 3a condição
	  else {
	    printf("Erro, condição inválida encontrada.\n");
	    exit(-1);
	  }
	}
      }
    }

  return(tau);

}

void zero(int s_0[], gsl_rng *r, double rho_0) {

  // Constrói um estado incial com uma concentração rho_0 de sítios ativos

  int i;

  for(i=0;i<n;i++) {
    if(gsl_rng_uniform(r)<rho_0) s_0[i]=1;
    else s_0[i]=0;
  }

}

void opiniao(int **rede, par *P, int s[], gsl_rng *r, char out[]) {
  // ---------------------------------------------------
  // Aqui é a parte da simulação da dinâmica de opinião.
  // ---------------------------------------------------

  // Obtendo os parâmetros
  int *h, k_max=0; // vetor de conectividade e com o estado dos sitios;
  double epsilon=P->eps, J=P->J, q=P->q;
  int i, w, t, tt, T=P->T_mcs, nt=P->D_mcs; // passo monte-carlo e iterador de sítios da rede
  double *tau, rho=0; // campo local de cada sítio e taca de transição
  FILE *f=fopen(out,"w"), *g;
  char ou[255];

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

  for(i=0;i<n;i++)
    rho+=s[i];
  rho/=((double)n);
  // Primeiro ponto na saida
  fprintf(f, "%lf\n", rho);

  // O monte-carlo começa aqui
  for(t=0;t<T;t+=nt) {
    g=fopen(ou,"w");
    fprintf(g,"# N=%d\tk=%d\tp=%lf\tsem=-1\tJ=%lf\teps=%lf\tq=%lf\tSR=%lu\n",n,K,p,J,epsilon,q,gsl_rng_get(r));
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
	if(gsl_rng_uniform(r)<tau[IDX(h[i],rede[i][0])]) {
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
    fclose(g);
  }
  
  fclose(f);

  free(tau);
  free(h);
  
}

int main(int argc, char *argv[]) {

  // ----------------------------------------------------------------
  // Nesta versão os estados iniciais são gerados para diversas redes
  // com o meso valor de p.
  // ----------------------------------------------------------------

  if(argc<16) {
    printf("Use %s <Dados de entrada:> <13. T mcs> <14. T mcs para salvar>  <15. saida s/ extensao>\n",argv[0]);
    printf("1. Dados de entrada: <n sítios>\n");
    printf("2. <k vizinhos pra frente>\n");
    printf("3. <semente>\n");
    printf("4. <epsilon>\n");
    printf("5. <J>\n");
    printf("6. <q>\n");
    printf("7. <densidade inicial>\n");
    printf("8. <densidade final>\n");
    printf("9. <# de condicoes inicias>");
    printf("10. <p inicial>\n");
    printf("11. <p final>\n");
    printf("12. <# de p's>");
    
    printf("Veja o cabeçalho do fonte para as definições dos parametros\n");
    exit(-1);
  }

  unsigned long int sem, sem_rede;
  int i, **rede, *s, *s0, num_c=atoi(argv[9]), num_p=atoi(argv[12]);
  double c_0i=strtod(argv[7],NULL), c_0f=strtod(argv[8],NULL), c_0, delta_c=(c_0f-c_0i)/(double)num_c, rho=0;
  double p_i=strtod(argv[10],NULL), p_f=strtod(argv[11],NULL), p, delta_p=(p_f-p_i)/(double)num_p;
  par P;
  char out[255], ai[255];
  gsl_rng *r;

  // Nomes para saída do modelo e entrada do gerador
  strcpy(out,argv[15]);

  // Copia os parâmetros da lista de argumetos
  P.T_mcs=atoi(argv[13]); P.D_mcs=atoi(argv[14]);
  n=atoi(argv[1]);
  K=atoi(argv[2]);
  sem=strtoul(argv[3],NULL,10);
  P.J=strtod(argv[5],NULL); P.eps=strtod(argv[4],NULL); P.q=strtod(argv[6],NULL);

  // Iniciando a semente de números aleatórios
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r,sem);

  // Mostra na tela os parâmetros do problema
  printf("\n# N=%d\tk=%d\tsem=%lu\tJ=%lf\teps=%lf\tq=%lf\n",
  	  n,K,sem,P.J,P.eps,P.q);

  // Alocando a estrutura de armazenamento da rede
  rede=(int **)malloc(n*Ip);
  for(i=0;i<n;i++) {
    rede[i]=(int *)calloc(1,I);
  }
  s=(int *)calloc(n,I);
  
  i=0;
  for(c_0=c_0i;c_0<c_0f;c_0+=delta_c) {
    for(p=p_i;p<p_f;p+=delta_p) {
      strcpy(ai,out);
      sprintf(ai,"%s_p=%g_%d.dat",ai,p,i); // para diferenciar os arquivos de saída
      makesw(rede,r);
      zero(s,r,c_0);
      opiniao(rede,&P,s,r,ai);
    }
    i++; // i indica a condição inicial
  }

  // Liberando memória
  for(i=0;i<n;i++) free(rede[i]);
  free(rede);
  free(s);
  gsl_rng_free(r);

  //printf("A semente foi %lu.\n", sem);

  return 0;

}
