#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

typedef struct rede {

  int *st;

  struct rede *prox;
  struct rede *ante;

} vertx;

// Estrutura de parâmetros do modelo
typedef struct parametros {
  // Rede
  int n, k, *s;
  vertx **rede;
  // Modelo de opinião
  double p, eps, J, q;
  int T_mcs, D_mcs;
  // Sementes de números aleatórios
  unsigned sem, sem_rede;
  // Entrada e saída de dados
  char out[255], fran[255];
}par;

int *alocavi(int n) {
  
  int *v;
  
  v=calloc(n,sizeof(int));
  if(v==NULL) {
    printf("Memoria insuficiente");
    exit(-1);
  }
  
  return (v);
  
}

vertx **alocar_rede(int n) {
  // ---------------------------------------
  // Alocar a estrutura que abrigará a rede.
  // ---------------------------------------

  int i;
  vertx **v;

  // Alocando:  
  if(n<1) {
    printf("\n**Parametros invalidos**\n");
    exit(-1);
  }
  v=calloc(n,sizeof(vertx));
  if(v==NULL) {
    printf("Memoria insuficiente");
    exit(-1);
  }
  for(i=0;i<n;i++) {
    v[i]=NULL;
    if(v==NULL) {
      printf("Memoria insuficiente");
      exit(-1);
    }
  }
  // Fim da alocação de memória

  return(v);

}

double *alocavd(int n) {
  
  double *v;
  
  v=calloc(n,sizeof(double));
  if(v==NULL) {
    printf("Memoria insuficiente");
    exit(-1);
  }
  
  return (v);
  
}

void regular1D(vertx **rede, int *x, int n, int k) {
  // -----------------------------------------------------------------
  // Essa função retorna um vetor unidimensional de n nós para um anel 
  // unidimensional com conectividade 2k.
  // -----------------------------------------------------------------

  int i, j;
  vertx *novo, *p; // novo vertice e atual vertice

  for(i=0;i<n;i++) {
    rede[i]=malloc(sizeof(vertx));
    rede[i]->st=&x[i]; // Linka o endereço de estado
    rede[i]->prox=NULL;
    rede[i]->ante=NULL;
    p=rede[i];
    for(j=0;j<k;j++) { // ligações pra frente
      //while(p->prox!=NULL) p=p->prox;
      novo=malloc(sizeof(vertx));
      novo->st=&x[(i+j+1)%n]; // a operação % coloca a condicao de contorno periodica
      novo->prox=NULL;
      novo->ante=p;
      p->prox=novo;
      p=p->prox; // salta pra proximo vertice
    }
    for(j=0;j<k;j++) { // ligações para trás
      //while(p->prox!=NULL) p=p->prox;
      novo=malloc(sizeof(vertx));
      novo->st=&x[(n-k+i+j)%n]; // a operação % coloca a condicao de contorno periodica
      novo->prox=NULL;
      novo->ante=p;
      p->prox=novo;
      p=p->prox; // salta pra proximo vertice
    }
  }

}

void free_rede(vertx **rede, int n) {
  // ---------------------------------
  // Libera o espaço do vetor de rede.
  // ---------------------------------

  int i;

  for(i=0;i<n;i++)
    free(rede[i]);

  free(rede);

}

void reagrupar(vertx **rede, int x[], int i, int j, int k) {
  // ------------------------------------------------
  // Essa função quebra a ligação i-j e constrói i-k.
  // ------------------------------------------------

  vertx *p, *novo;
  int a=0;
   
  // se i-k ja existe não ocorre nada. Versão Nova!
  for(p=rede[i];p!=NULL;p=p->prox)
    if(p->st==&x[k]) return;

  // no sitio i apenas o rotulo de j muda para apontar pra k
  for(p=rede[i];p!=NULL;p=p->prox)
    if(p->st==&x[j]) {
      p->st=&x[k];
      a=1;
      break;
    }
  if(a==0) {
    printf("%d-%d nao existe\n",i,j);
    exit(-1);
  }
  a=0;

  // no sitio j a ligacao com i é perdida
  if(rede[j]->prox!=NULL) {
    for(p=rede[j]->prox;p!=NULL;p=p->prox) {
      if(p->st==&x[i]) {
	if(p->ante!=NULL) (p->ante)->prox=p->prox;
	if(p->prox!=NULL) (p->prox)->ante=p->ante;
	free(p);
	a=1;
	break;
      }
    }
  }
  if(a==0) {
    printf("%d-%d nao existe\n",j,i);
    exit(-1);
  }
  // k ganha a ligação de i
  if(rede[k]->prox!=NULL) { // k ganha mais um vizinho
    p=rede[k]->prox;
    while(p->prox!=NULL) p=p->prox;
    novo=malloc(sizeof(vertx));
    novo->st=&x[i];
    novo->prox=NULL;
    novo->ante=p;
    p->prox=novo;
    p=p->prox;
    }

}

void doSW(vertx **rede, int x[], int n, int k, double p, unsigned sem) {

  int i, j, l;

  //srand(234556431);
  srand(sem);
  for(i=0;i<n;i++)
    for(j=0;j<k;j++) {
      if(rand()/(double)RAND_MAX<p) {
	//print_rede(rede,n); getchar();
	l=(i+k+1+rand()%(n-k-1))%n;
	reagrupar(rede,x,i,(i+j+1)%n,l);
	//print_rede(rede,10); getchar();
      }
    }
    
}

int *grau_vertices(vertx **rede, int n) {
  // -------------------------------------
  // Conta o grau de cada vértice da rede.
  // -------------------------------------

  int i, *k;
  vertx *p;

  k=alocavi(n);

  for(i=0;i<n;i++)
    for(p=rede[i]->prox;p!=NULL;p=p->prox)
      k[i]++;

  return(k);

}

double *make_tau(int n, double J, double q, double epsilon) {
  // --------------------------------------
  // Define as taxas de transica do modelo.
  // --------------------------------------

  int s, k;
  double *tau, h;

  tau=alocavd((n+1)*(n+2)/2);

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

void salvar_estado(int n, int *state, FILE *out) {

  int i;

  for(i=0;i<n;i++)
    fprintf(out, "%d\n", state[i]);

}

void opiniao(par *P) {
  // ---------------------------------------------------
  // Aqui é a parte da simulação da dinâmica de opinião.
  // ---------------------------------------------------

  // Obtendo os parâmetros
  int n=P->n, *k, k_max, *s=P->s, *h; // vetor de conectividade e com o estado dos sitios;
  vertx **rede=P->rede;
  double epsilon=P->eps, J=P->J, q=P->q;
  int i, t, tt, T=P->T_mcs, nt=P->D_mcs; // passo monte-carlo e iterador de sítios da rede
  unsigned long int nr; // Variável para receber o aleatório inteiro.
  double *tau, prob, rho=0; // campo local de cada sítio e taca de transição
  vertx *p; // iterador de sítios da rede
  FILE *f=fopen(P->out,"w"), *g, *ff=fopen(P->fran,"r");
  char out[255], fran[255];
  // Variáveis do gerador
  gsl_rng *r;

  strcpy(out,P->out);
  strcpy(fran,P->fran);
  strcat(out,"~");

  // Iniciando a galera
  h=alocavi(n);
  k=grau_vertices(rede,n);
  // Pegar o grau máximo
  k_max=0;
  for(i=0;i<n;i++)
    if(k[i]>k_max) k_max=k[i];
  
  // Constrói as taxas com apenas o grau máximos
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
    g=fopen(out,"w");
    ff=fopen(fran,"w");
    fprintf(g,"# N=%d\tk=%d\tp=%lf\tsem=%lu\tJ=%lf\teps=%lf\tq=%lf\tSR=%u\n",n,P->k,P->p,nr,J,epsilon,q,P->sem_rede);
    fprintf(g, "# Semente: %lu\n", nr);
    for(tt=0;tt<nt;tt++) {
            
      // Calcular o vetor h dos sítios
      for(i=0;i<n;i++) {
	h[i]=0;
	for(p=rede[i]->prox;p!=NULL;p=p->prox) {
	  h[i]+=*p->st;
	}
      }

      // Atualiza os sítios
      for(i=0;i<n;i++) {
	nr=rand();
	prob=gsl_rng_uniform(r);
	if(prob<tau[h[i]-1+(k[i]+1)*k[i]/2]) {
	  if(s[i]==0) s[i]=1;
	}
	else
	  if(s[i]==1) s[i]=0;
      }

      // Calcula a densidade
      rho=0;
      for(i=0;i<n;i++)
	rho+=s[i];
      rho/=((double)n);
      
      // Imprime saída
      fprintf(f, "%lf\n", rho);
      
    }
    salvar_estado(n,s,g);
    gsl_rng_fwrite(ff,r);
    fclose(g);
    fclose(ff);
  }
  
  free(tau);
  free(h);
  free(k);
  // Libera memória do gerador
  gsl_rng_free (r);
  
}

int main(int argc, char *argv[]) {

  if(argc<6) {
    //             1                2       3               4                 5
    // 6 é um arquivo oculto .<saida densidade>
    printf("Use %s <estado inicial> <T mcs> <N para salvar> <saida densidade> <arquivo com estado do gerador gsl>\n", argv[0]);
    printf("Veja o cabeçalho do fonte para as definições dos parametros\n");
    exit(-1);
  }

  int i, *s, si;
  vertx **sw;
  par P;
  FILE *in=fopen(argv[1],"r");

  if(in==NULL) {
    printf("Arquivo de estado inicial inválido.\n");
    exit(-1);
  }

  // Coletando informações a partir do arquivo de dados
  fscanf(in,"# N=%d\tk=%d\tp=%lf\tsem=%u\tJ=%lf\teps=%lf\tq=%lf\tSR=%u\n",&P.n,&P.k,&P.p,&P.sem,&P.J,&P.eps,&P.q,&P.sem_rede);
  fscanf(in, "# Semente: %u\n", &P.sem);

  // Coletando estado inicial
  s=alocavi(P.n);
  for(i=0;i<P.n;i++) {
    fscanf(in,"%d",&s[i]);
  }

  // Entrada a partir da linha de comando
  P.D_mcs=atoi(argv[3]);
  P.T_mcs=atoi(argv[2]);
  strcpy(P.out,argv[4]);
  strcpy(P.fran,argv[5]);
  
  fclose(in);

  // Mostra na tela os parâmetros do problema
  //printf("\n# N=%d\tk=%d\tp=%lf\tsem=%d\tJ=%lf\teps=%lf\tq=%lf\tSR=%u\n",
  //	  P.n,P.k,P.p,P.sem,P.J,P.eps,P.q,P.sem_rede);

  // Montar substrato
  sw=alocar_rede(P.n);
  regular1D(sw,s,P.n,P.k);
  doSW(sw,s,P.n,P.k,P.p,P.sem_rede);

  P.rede=sw; // Apontar para a rede
  P.s=s; // Estado da rede

  //for(i=0;i<10;i++)
  //  printf("%d ",P.s[i]);
  //printf("\n");
  
  // Experimento
  opiniao(&P);

  free_rede(sw,P.n);
  free(s);
  
  return 0;


}
