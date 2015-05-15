/*
 * pmedian
 * autor: Fabricio Olivetti de França
 * ra: 027783
 *
 * para compilar em ambiente com gcc/make
 * instalados, basta digitar 'make'
 *
 * modo de uso:
 * pmedian [nome_arquivo_problemas] [tipo de algo] [tipo relatorio] [iteracoes]
 * 
 * [tipo de algo]: 1 - heur. construtiva, 2 - busca local, 3 - busca tabu, 4 - alg. genetico
 * [tipo de relatorio]: 1 - relatorio simples, 2 - relatorio completo
 * [iteracoes]: quantas iteracoes utilizar na heur. construtiva
 */

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "pmedian.h"

#define CLOCKS_PER_SECOND 1000

typedef unsigned long long int64;
int64 clockCycles;
inline void CountCycles()
{
asm volatile ("rdtsc":"=A"(clockCycles));
}

int main(int argc, char **argv){

	int i,j,k,n_probs,algo,dummy,best,vertex,medians,cap,T,osman;
	int rel=1,max_it=1;
	FILE *arq_dados;
	costumer *costum;
	double **dist, answer, *rpc,soma;
	int **x, *y;
	char tipo;

	unsigned long long cyc_start,cyc_stop;

	//clock_t time_start,time_stop;

	/* inicializa parametros defaults caso necessario */
	srand((int)time(NULL));
	if(argc < 2){
		printf("Uso: pmedian arq [algo] [rel] [max_it]\n\nOnde:\n\tarq e´ o nome do arquivo de dados\n\t \
				talgo:\n\t\t1 - heuristica construtiva\n\t\t2 - heuristica busca local\n\t\t \
				t3 - algoritmo genetico\n\t\t4 - busca tabu\n\t	\
				rel:\n\t\t1 - relatorio simples\n\t\t2 - relatorio completo\n\tmax_it - qtde. de \
				iteracoes para heuristica construtiva.");
		exit(0);
	}else if(argc == 2){
		algo = 1;
	}
	if(argc > 2){
		sscanf(argv[2],"%d",&algo);
	}
	if(argc > 3){
		sscanf(argv[3],"%d",&rel);
	}
	if(argc > 4){
		sscanf(argv[4],"%d",&max_it);
	}
	if(argc > 5){
		sscanf(argv[5],"%d",&T);
	}

	if((arq_dados = fopen(argv[1],"r")) == NULL){
		perror("fopen - Arquivo inexistente\n");
		exit(-1);
	}

	if(strcmp(strchr(argv[1],'.'),".pmd")==0){
		osman=0;
	}else{
		osman=1;
	}
	
	if(!osman)
		n_probs=1;
	else
		fscanf(arq_dados,"%d",&n_probs);

	rpc = (double *)malloc(sizeof(double)*n_probs);
	/* para cada problema ... */
	for(i=0;i<n_probs;i++){

		if(!osman){
			/* rodrigo */

			fscanf(arq_dados,"%*s\n%*s %d",&best);
			fscanf(arq_dados,"%*s %c",&tipo);
			fscanf(arq_dados,"%*s %d",&vertex);
			fscanf(arq_dados,"%*s %d",&medians);
			fscanf(arq_dados,"%*s %d",&cap);
			fscanf(arq_dados,"%*s");
		
			costum = (costumer *) malloc(sizeof(costumer)*vertex);
	
			for(j=0;j<vertex;j++)
				fscanf(arq_dados,"%d",&costum[j].demand);
			fscanf(arq_dados,"%*s");
	
			dist = (double **) malloc(sizeof(double *)*vertex);
			for(j=0;j<vertex;j++){
				dist[j] = (double *) malloc(sizeof(double)*vertex);
				for(k=0;k<vertex;k++){
					fscanf(arq_dados,"%lf",&dist[j][k]);
				}
			}

	
		}else{
			/* osman e lorena */
		
			/* le as variaveis necessarias */
			fscanf(arq_dados,"%d %d",&dummy,&best);
			fscanf(arq_dados,"%d %d %d",&vertex,&medians,&cap);
			/*aloca espaco necessario */
			costum = (costumer *) malloc(sizeof(costumer)*vertex);
		
			for(j=0;j<vertex;j++){
				fscanf(arq_dados,"%d %d %d %d",&dummy,&costum[j].x,&costum[j].y,&costum[j].demand);
			}
	
			/* calcula as distancias e coloca em uma tabela arredondando para baixo como em Osman */
			dist = (double **) malloc(sizeof(double *)*vertex);
			for(j=0;j<vertex;j++){
				dist[j] = (double *) malloc(sizeof(double)*vertex);
				for(k=0;k<vertex;k++){
					dist[j][k] = floor(sqrt(pow((costum[j].x-costum[k].x),2)+pow((costum[j].y-costum[k].y),2)));
				}
			}
		}
	

		x = (int **)malloc(sizeof(int *)*vertex);
		y = (int *)malloc(sizeof(int)*vertex);
		for(j=0;j<vertex;j++)
			x[j] = (int *)malloc(sizeof(int)*vertex);

		for(k=0;k<1;++k){
		/* manda resolver o problema */
		CountCycles();
		cyc_start = clockCycles; 
		if(algo==1){
			cheur(costum,vertex,dist,medians,cap,x,y,max_it);
		}else if(algo==2){
                        cheur(costum,vertex,dist,medians,cap,x,y,1);
			lsheur2(costum, vertex, dist, medians, cap, x, y, max_it);
		}else if(algo==3){
                        cheur(costum,vertex,dist,medians,cap,x,y,1);
			ts(costum, vertex, dist, medians, cap, x, y, T,max_it);
			//ts2(costum, vertex, dist, medians, cap, x, y, T,max_it);
			//ts4(costum, vertex, dist, medians, cap, x, y, T,max_it);
			
		}else if(algo==4){
			ag(costum, vertex, dist, medians, cap, x, y, max_it);
		}else if(algo==5){
				aco(costum, vertex, dist, medians, cap, x, y, max_it);
		}else if(algo==6){
                        cheur(costum,vertex,dist,medians,cap,x,y,1);
			geo(costum, vertex, dist, medians, cap, x, y, max_it);
		}
		
		CountCycles();
		cyc_stop = clockCycles; 

		printf("tempo em seg.: %f\n",(float)(cyc_stop-cyc_start)/2000000000.0);
		printf("factivel: %d\n",factivel(costum,vertex,medians,cap,x,y));
		/* calcula e imprime o resultado */
		answer = calculate(costum,vertex,dist,medians,cap,x,y);
		
		/* imprime um relatorio de acordo com o parametro (1=simples, 2=completo) */
		switch(rel){
			case 1:
				rpc[i] = print_short(i,best,answer);
				break;
			case 2:
				rpc[i] = print_report(i,best,answer,costum,vertex,dist,medians,cap,x,y);
				break;
		}
		}

		/* limpa a memoria utilizada */
		free(costum);
		free(y);
		for(j=0;j<vertex;j++){
			free(dist[j]);
			free(x[j]);
		}
		free(dist);
		free(x);
	}
	fclose(arq_dados);

	for(i=0,soma=0;i<n_probs;++i) soma += rpc[i];
	printf("\n\nRPC medio: %f\n",soma/n_probs);
	free(rpc);

	return 0;
}

double calculate(costumer *costum, int vertex, double **dist, int medians, int cap, int **x, int *y){

	double answer=0.0;
		int i,j; //,soma,med;

	/* para cada vertice...*/
	for(i=0;i<vertex;++i){
		/* ...para os outros vertices.. */
		for(j=0;j<vertex;++j){
			/* ... se for cliente->mediana, soma a distancia*/
			if(x[i][j]){
				answer += (dist[i][j]);
			}
		}
	}

	
	return answer;

}

double print_report(int it,int best,double answer,costumer *costum, int n, double **dist, int medians, int cap, int **x, int *y){

	/*respostas do Osman*/
	static int hoc[] = {786,816,972,891,804,882,968,945,752,1017,1761,1567,1847,1635,1517,1780,1665,1345,1634,1872};
	int i,j,primeiro=0,soma,med=-1;
	double rpc;

	rpc = 100.0*(answer - (double)best)/(double)best;
	
	printf("\nExperimento %d:\n ",it+1);
	for(i=0;i<60;++i) printf("-");
	printf("\n| medianas: {");

	for(i=0;i<n;++i){
		if(y[i]){
			if(primeiro) printf(", %d",i);
			else{ printf("%d",i);primeiro=1;}
		}
	}
	printf("}\n ");
	for(i=0;i<60;++i) printf("-");
	printf("\n|alocacoes:\n|\n|");
	for(i=0;i<n;++i){
		if(y[i]){
			printf("mediana %d: {",i);
			primeiro=0;
			for(j=0;j<n;++j){
				if(x[j][i]){
					if(primeiro){ printf(", %d",j);
					}else{ printf("%d",j);primeiro=1;}
				}
			}
			printf("}\n|");
		}
	}
	printf("\n ");
	for(i=0;i<60;++i) printf("-");
	printf("\n|Uso das medianas: \n|\t\tcapac.\tusado\n|");
	for(i=0;i<medians;++i){
		soma=0;
		for(j=med+1;j<n;++j){
			if(y[j]){
				med=j;
				break;
			}
		}
		for(j=0;j<n;++j){
			if(x[j][med]){
				soma+=costum[j].demand;
			}
		}
		printf("mediana %d:\t%d\t%d\n|",med,cap,soma);
	}
	printf("\n ");
	for(i=0;i<60;++i) printf("-");
	printf("\n|dist. total\tH.OC (Osman)\tMelhor\tRPC\n|%.0f\t\t%d\t\t%d\t\t%f",answer,hoc[it],best,rpc);
	printf("\n ");
	for(i=0;i<60;++i) printf("-");
	printf("\n");

	return rpc;
}

double print_short(int it,int best,double answer){
	
	static int hoc[] = {786,816,972,891,804,882,968,945,752,1017,1761,1567,1847,1635,1517,1780,1665,1345,1634,1872};
	int i;
	double rpc;

	rpc = 100.0*(answer - (double)best)/(double)best;

	printf("\nExperimento %d:\n ",it+1);
	for(i=0;i<60;++i) printf("-");
	printf("\n|dist. total\tH.OC (Osman)\tMelhor\tRPC\n|%.0f\t\t%d\t\t%d\t\t%f",answer,hoc[it],best,rpc);
	printf("\n ");
	for(i=0;i<60;++i) printf("-");
	printf("\n");

	return rpc;
	
}
