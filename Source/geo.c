#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pmedian.h"

typedef struct s_geo{
	int entra;
	int sai;
	double fitness;
	int rank;
}s_geo;

void geo(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it){

	s_geo *sols;
	int it, i,j,k, *meds,l,m, velho;
	int *ordem;
	double tal=1.5,answer;

	sols = (s_geo *)malloc(sizeof(s_geo)*(n-medians)*medians);
	ordem = (int *)malloc(sizeof(int)*(n-medians)*medians);
	meds = (int *)malloc(sizeof(int)*medians);

	/* info do best */
	int **best_x, *best_y;
	double best_answer;

	best_x = (int **)malloc(sizeof(int *)*n);
	best_y = (int *)malloc(sizeof(int)*n);
	for(i=0;i<n;++i)
		best_x[i] = (int *)malloc(sizeof(int)*n);

	for(i=0,j=0;i<n;++i)
		if(y[i])
			meds[j++]=i;

	for(l=0;l<n;++l){
		best_y[l]=y[l];
		for(m=0;m<n;++m)
			best_x[l][m]=x[l][m];
	}
	best_answer = calculate(costu, n, dist, medians, cap, x, y);

	/* a partir de uma solução inicial */
	
	for(it=0;it<max_it;++it){
		
		/* 1- avalia o fitness que se obtem para cada caso de troca */
		for(i=0,k=0;i<medians;++i){
			for(j=0;j<n;++j){
				if(!y[j]){
					sols[k].entra = j;
					sols[k].sai = i;
					
					y[meds[i]]=0;
					velho = meds[i];
					y[j]=1;
					meds[i]=j;
					for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
					designar(costu, medians,n,dist,cap,x,meds);
					
					sols[k++].fitness = calculate(costu, n, dist, medians, cap, x, y);
					
					meds[i] = velho;
					y[velho] = 1;
					y[j] = 0;
				}
			}
		}
		
		/* 2 - rankeia as solucoes */
		for(i=0;i<(n-medians)*medians;++i){
			ordem[i]=i;
		}
		for(i=0;i<(n-medians)*medians;++i){
			for(j=i+1;j<(n-medians)*medians;++j){
				if(sols[ordem[i]].fitness > sols[ordem[j]].fitness){
					ordem[i]^=ordem[j]^=ordem[i]^=ordem[j];
				}
			}
		}
		for(i=0;i<(n-medians)*medians;++i){
			sols[ordem[i]].rank=i;
		}

		/* 3 - para cada solucao sorteia um numero aleatorio que, 
		 * se for <= ao rank efetua-se a troca e volta para 1 
		 * caso contrario, repete 3 
		 * */
//		for(i=0;i<(n-medians)*medians;++i){
		while(1){

			i = (int) floor((double)rand()*(((n-medians)*medians)-1)/RAND_MAX);
//			printf("%f\n",(double)rand()/RAND_MAX);
			if(pow(sols[i].rank,-tal) >= (double)rand()/RAND_MAX){
				y[meds[sols[i].sai]]=0;
				y[sols[i].entra]=1;
				meds[sols[i].sai]=sols[i].entra;
				for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
				designar(costu, medians,n,dist,cap,x,meds);
					
				answer = calculate(costu, n, dist, medians, cap, x, y);
				if(answer < best_answer && factivel(costu,n,medians,cap,x,y)){
					for(l=0;l<n;++l){
						best_y[l]=y[l];
						for(m=0;m<n;++m)
							best_x[l][m]=x[l][m];
					}
					best_answer = answer;
					printf("%d - %f\n",it,answer);
				}
//				printf("ans: %f - %d/%d\n",answer,sols[i].sai,sols[i].entra);

				break;
			}
		}
	}
	for(l=0;l<n;++l){
		y[l]=best_y[l];
		for(m=0;m<n;++m)
			x[l][m] = best_x[l][m];
	}
	lsheur2(costu,n,dist,medians,cap,x,y,5);

	free(sols);
	free(meds);
	free(ordem);
	free(best_y);
	for(i=0;i<n;++i)
		free(best_x[i]);
	free(best_x);
}
