#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ---------------------------------------------------------------
// Modelo de opinião em redes small-world, beseado no preprint:
// arXiv:0909.0117v3 [physics.soc-ph] 23 Dec 2009, Bagnoli et al.
//
// q) é o parâmetro de pressão social para os estados absorventes
//    0 para 1-q e 1 para q.
// J) é a interação ferro antiferro que estabiliza ou não os esta-
//    dos absorventes.
// epsilon) é o termo de discretização do espaço, dirá o quão dis-
// tantes são os sítios um dos outros.
// ---------------------------------------------------------------


// Estrutura que definirá o substrato
struct vertice {
  int rotulo; // esse é o rotulo do dito cujo
  struct vertice *prox; // a lista de adjacencia começa aqui
  struct vertice *ante; // para que a lista possa ter elementos do meio excuidos
};

typedef struct vertice no;

void print_rede(no **rede, int n) {

  int i;
  no *p;

  for(i=0;i<n;i++)
    for(p=rede[i]->prox;p!=NULL;p=p->prox)
      printf("%d %d\n",i,p->rotulo);

}

no **alocar_rede(int n) {
  // ---------------------------------------
  // Alocar a estrutura que abrigará a rede.
  // ---------------------------------------

  int i;
  no **v;

  // Alocando:  
  if(n<1) {
    printf("\n**Parametros invalidos**\n");
    exit(-1);
  }
  v=calloc(n,sizeof(no));
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

void free_rede(no **rede, int n) {
  // ---------------------------------
  // Libera o espaço do vetor de rede.
  // ---------------------------------

  int i;

  for(i=0;i<n;i++)
    free(rede[i]);

  free(rede);

}

int *alocavi(int n) {
  
  int *v;
  
  v=calloc(n,sizeof(int));
  if(v==NULL) {
    printf("Memoria insuficiente");
    exit(-1);
  }
  
  return (v);
  
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

void regular1D(no **rede, int n, int k) {
  // -----------------------------------------------------------------
  // Essa função retorna um vetor unidimensional de n nós para um anel 
  // unidimensional com conectividade 2k.
  // -----------------------------------------------------------------

  int i, j;
  no *novo, *p; // novo vertice e atual vertice

  for(i=0;i<n;i++) {
    rede[i]=malloc(sizeof(no));
    rede[i]->rotulo=i;
    rede[i]->prox=NULL;
    rede[i]->ante=NULL;
    p=rede[i];
    for(j=0;j<k;j++) { // ligações pra frente
      //while(p->prox!=NULL) p=p->prox;
      novo=malloc(sizeof(no));
      novo->rotulo=(i+j+1)%n; // a operação % coloca a condicao de contorno periodica
      novo->prox=NULL;
      novo->ante=p;
      p->prox=novo;
      p=p->prox; // salta pra proximo vertice
    }
    for(j=0;j<k;j++) { // ligações para trás
      //while(p->prox!=NULL) p=p->prox;
      novo=malloc(sizeof(no));
      novo->rotulo=(n-k+i+j)%n; // a operação % coloca a condicao de contorno periodica
      novo->prox=NULL;
      novo->ante=p;
      p->prox=novo;
      p=p->prox; // salta pra proximo vertice
    }
  }

}

void reagrupar(no **rede, int i, int j, int k) {
  // ------------------------------------------------
  // Essa função quebra a ligação i-j e constrói i-k.
  // ------------------------------------------------

  no *p, *novo;
  int a=0;
  
  // se i-k ja existe não ocorre nada.
  for(p=rede[i];p!=NULL;p=p->prox)
    if(p->rotulo==k) return;

  // no sitio i apenas o rotulo de j muda para apontar pra k
  for(p=rede[i];p!=NULL;p=p->prox)
    if(p->rotulo==j) {
      p->rotulo=k;
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
      if(p->rotulo==i) {
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
    novo=malloc(sizeof(no));
    novo->rotulo=i;
    novo->prox=NULL;
    novo->ante=p;
    p->prox=novo;
    p=p->prox;
    }

}

void doSW(no **rede, int n, int k, double p) {

  int i, j, l;

  srand(234556431);
  for(i=0;i<n;i++)
    for(j=0;j<k;j++) {
      if(rand()/(double)RAND_MAX<p) {
	//print_rede(rede,n); getchar();
	l=(i+k+1+rand()%(n-k-1))%n;
	reagrupar(rede,i,(i+j+1)%n,l);
	//print_rede(rede,10); getchar();
      }
    }
    
}

int *grau_vertices(no **rede, int n) {
  // -------------------------------------
  // Conta o grau de cada vértice da rede.
  // -------------------------------------

  int i, *k;
  no *p;

  k=alocavi(n);

  for(i=0;i<n;i++)
    for(p=rede[i]->prox;p!=NULL;p=p->prox)
      k[i]++;

  return(k);

}

double tau(double J, double epsilon, double q, double h) {
  // --------------------------------------
  // Define as taxas de transica do modelo.
  // --------------------------------------

  if(h<q) return (epsilon);
  if(h>(1-q)) return (1-epsilon);
  if((q<=h)&&(h<=1-q)) return (1.0/(1.0+exp(-2*J*(2*h-1))));

  return epsilon;

}

void opiniao(int n, no **rede, int T, double semente, double J, double epsilon, double q, FILE *out) {
  // ---------------------------------------------------
  // Aqui é a parte da simulação da dinâmica de opinião.
  // ---------------------------------------------------

  int t=0, i; // passo monte-carlo e iterador de sítios da rede
  int *k, *s; // vetor de conectividade e com o estado dos sitios;
  double *h, prob, rho=0; // campo local de cada sítio e taca de transição
  no *p; // iterador de sítios da re
  
  // Iniciando a galera
  s=alocavi(n);
  h=alocavd(n);
  k=grau_vertices(rede,n);
  srand(semente);

  // Estado inicial aleatório
  for(i=0;i<n;i++)
    s[i]=rand()%2;
  
  for(i=0;i<n;i++)
    rho+=s[i];
  rho/=((double)n);
  // Primeiro ponto na saida
  fprintf(out, "%d %lf\n", t, rho);

  // O monte-carlo começa aqui
  for(t=0;t<T;t++) {

    rho=0;

    // Calcular o vetor h dos sítios
    for(i=0;i<n;i++) {
      for(p=rede[i]->prox;p!=NULL;p=p->prox)
	h[i]+=s[p->rotulo];
      h[i]=h[i]/(double)k[i];
    }

      // Atualiza os sítios
    for(i=0;i<n;i++) {
      prob=rand()/(double)RAND_MAX;
      if(prob<tau(J,epsilon,q,h[i])) {
	if(s[i]==0) s[i]=1;
      }
      else
	if(s[i]==1) s[i]=0;
    }

    // Calcula a densidade
    for(i=0;i<n;i++)
      rho+=s[i];
    rho/=((double)n);

    // Imprime saída
    fprintf(out, "%d %lf\n", t, rho);

  }

  free(s);
  free(h);
  free(k);
  
}

int main(int argc, char *argv[]) {

  if(argc<10) {
    printf("Use %s <n sítios> <k vizinhos pra frente> <p WS-model> <T passos MC> <semeste> <epsilon> <J> <q> <saida>\n",argv[0]);
    printf("Veja o cabeçalho do fonte para as definições dos parametros\n");
    exit(-1);
  }

  no **sw;
  int n=atoi(argv[1]), k=atoi(argv[2]), T=atoi(argv[4]), sem=atoi(argv[5]);
  double p=atof(argv[3]),eps=atof(argv[6]), J=atof(argv[7]), q=atof(argv[8]);
  FILE *f=fopen(argv[9],"w");

  // Cabeçalho da saída
  printf("# N=%d\tk=%d\tp=%lf\tsem=%d\t=J%lf\teps=%lf\tq=%lf\t\n",
	  n,k,p,sem,J,eps,q);
  fprintf(f, "# N=%d\tk=%d\tp=%lf\tsem=%d\t=J%lf\teps=%lf\tq=%lf\t\n",
	  n,k,p,sem,J,eps,q);

  // Montar substrato
  sw=alocar_rede(n);
  regular1D(sw,n,k);
  doSW(sw,n,k,p);

  // Experimento
  opiniao(n,sw,T,sem,J,eps,q,f);

  free(sw);
  fclose(f);  
  

  return 0;

}
