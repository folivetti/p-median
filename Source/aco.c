#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pmedian.h"

typedef struct ants{
	int *sol;
	double fitness;
	int **x;
	int *y;
}ants;

void density(costumer *costu, int n, double **dist, int cap, double *dens, int *cids);
void ant_copy(ants *ant,int **x,int *y,int n);

/*
 *
 * alpha
 * beta
 * phero
 * ni = density
 * 
 * */
void aco(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it){

	double *dens,soma;
	double *phero, *dtal, ro=0.1;
	double p,rand_p;
	double tal_max, tal_min;
	int i,j,k,it,n_ants=10;
	int *cids, *sol;
	ants *ant, best;
	int melhor_i;
	double melhor;
	double alpha=1, beta=0, soma_pher, *pa,soma_a,m_dist,soma_fit;
	int *meds, its_n=0;
	double soma_p,soma_p_novo;

	dens = (double *)malloc(sizeof(double)*n);
	cids = (int *)malloc(sizeof(int)*n);
	phero = (double *)malloc(sizeof(double)*n);
	meds =	(int *)malloc(sizeof(int)*medians);
	sol =	(int *)malloc(sizeof(int)*n);
	dtal =	(double *)malloc(sizeof(double)*n);
	pa =	(double *)malloc(sizeof(double)*n);
	ant = (ants *)malloc(sizeof(ants)*n_ants);

	m_dist=0;
	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			if(m_dist<dist[i][j])
				m_dist=dist[i][j];

	//tal_max = (1/(1-ro))*(1/(n*m_dist));
	//tal_min = tal_max/(2*n);
	tal_max = 0.9999;
	tal_min = 0.0001;

	for(i=0;i<n_ants;++i){
		ant[i].sol = (int *)malloc(sizeof(int)*medians);
		ant[i].y = (int *) malloc(sizeof(int)*n);
		ant[i].x = (int **) malloc(sizeof(int *)*n);
		for(j=0;j<n;++j){
			ant[i].x[j] = (int *) malloc(sizeof(int)*n);
		}
	}
	best.sol = (int *)malloc(sizeof(int)*medians);
	best.y = (int *) malloc(sizeof(int)*n);
	best.x = (int **) malloc(sizeof(int *)*n);
	for(j=0;j<n;++j){
		best.x[j] = (int *) malloc(sizeof(int)*n);
	}

	best.fitness = 9e+10;
	
	density(costu, n, dist, cap, dens,cids);


	for(i=0;i<n;++i)
		//phero[i]=tal_max;
		phero[i]=0.5;
soma_p=0.0;
soma_p_novo=(double)(n-medians)*tal_min*10.0;
	for(it=0;it<max_it;++it){

		//if(its_n>50*max_it/100){
		//printf("%f\n",(double)(n-medians)*tal_min*10);
		if(soma_p<((double)medians+soma_p_novo) && soma_p>((double)medians-soma_p_novo)){
		//if( ((soma_p < ((double)medians+0.3)) && (soma_p > ((double)medians-0.3))) && abs(soma_p_novo - soma_p) < 0.01){
		//	printf("%f %f\n",soma_p,soma_p_novo);
			for(i=0;i<n;++i)
				//phero[i]=tal_max;
				phero[i]=0.5;	
			its_n=0;
			printf("%d change\n",it);
		}
		//soma_p_novo=soma_p;
		soma_pher=0.0;
		for(i=0;i<n;++i){
			soma_pher += pow(dens[i],beta)*pow(phero[i],alpha);
		}

		/* para cada formiga */
		for(i=0;i<n_ants;++i){

			for(j=0;j<n;++j){
				y[j]=0;
				for(k=0;k<n;++k){
					x[j][k]=0;
				}
			}
			for(j=0;j<n;++j)
				sol[j]=0;

			/* constroi uma solucao nova usando a densidade e feromonio */
			for(j=0;j<medians;++j){

				soma_a = 0.0;
				for(k=0;k<n;++k){
					pa[k] = pow(dens[k],beta)*pow(phero[k],alpha)/soma_pher;
					if(!sol[k])
						soma_a += pa[k];
				}
				rand_p = (double)rand()/RAND_MAX;
				for(k=0,p=0.0;k<n && p<rand_p;++k){
					if(sol[k]) continue;
					//printf("%f %f %f\n",pa[k],soma_a,pa[k]/soma_a);
					p += pa[k]/soma_a;
				}
				--k;
				k=k%n;
				while(sol[k]) --k;
				if(k<0) k=0;
				else if(k>n) k=n-1;
				ant[i].sol[j] = k;
				y[k]=1;
				meds[j]=k;
				sol[k]=1;
				soma_pher -=pow(dens[k],beta)*pow(phero[k],alpha);
			}
			
			/* avalia cada solucao */
			designar(costu, medians,n,dist,cap,x,meds);
			
			lsheur2(costu,n,dist,medians,cap,x,y,1);
			ant_copy(&ant[i],x,y,n);
			
			ant[i].fitness = calculate(costu,n,dist,medians,cap,x,y);
			for(j=0,k=0;j<n && k<medians;++j)
				if(y[j])
					ant[i].sol[k++]=j;
			
			soma = 0;
			for(j=0;j<n;++j)
				for(k=0;k<n;++k)
					soma+=x[j][k];
			soma = n-soma;
			ant[i].fitness += soma*m_dist;		
		}
		
		/* verifica a melhor formiga da solucao */
		melhor = ant[0].fitness;
		//soma_fit = 1.0/ant[0].fitness + 1.0/best.fitness;
		melhor_i = 0;
		for(i=1;i<n_ants;++i){
			if(ant[i].fitness < ant[melhor_i].fitness){
				melhor = ant[i].fitness;
				melhor_i = i;
			}
			//soma_fit += 1.0/ant[i].fitness;
		}
		
		/* zera delta_tal */
		for(i=0;i<n;++i)
			dtal[i] = 0.0;
		
		/* atualiza melhor solucao global */
		if(ant[melhor_i].fitness < best.fitness){
			for(i=0;i<medians;++i){
				best.sol[i] = ant[melhor_i].sol[i];
			}
			for(i=0;i<n;++i){
				best.y[i] = ant[melhor_i].y[i];
				for(j=0;j<n;++j){
					best.x[i][j] = ant[melhor_i].x[i][j];
				}
			}
			best.fitness = ant[melhor_i].fitness;
			printf("%d melhorou: %f %f\n",it,best.fitness,soma_p);
			its_n=0;
		}else{
			++its_n;
		}

		/* calula delta_tal */
		for(i=0;i<medians;++i){
			//dtal[best.sol[i]] += 1/best.fitness;
			//dtal[ant[melhor_i].sol[i]] += 1/ant[melhor_i].fitness;
			//dtal[best.sol[i]] += (1.0/best.fitness)/soma_fit;
			//printf("%d %f\n",i,(1.0/best.fitness)/soma_fit);
			//dtal[ant[melhor_i].sol[i]] += (1.0/ant[melhor_i].fitness)/soma_fit;
			dtal[best.sol[i]] += 0.7*(best.fitness - ant[melhor_i].fitness)/best.fitness;
			dtal[ant[melhor_i].sol[i]] += 1 - (best.fitness - ant[melhor_i].fitness)/best.fitness;
		}

		/* atualiza trilha de feromonio */
		for(i=0;i<n;++i){
			phero[i] = phero[i] + ro*(dtal[i]-phero[i]);
			//printf("%d %f\n",i,phero[i]);
			//phero[i] = ro*phero[i] + dtal[i];
			if(phero[i] < tal_min)
				phero[i] = tal_min;
			else if(phero[i] > tal_max)
				phero[i] = tal_max;
		}
		
		for(i=0,soma_p=0.0;i<n;++i){
			soma_p += phero[i];
		}
		//printf("it %d - %f\n",it,soma_p);
	}


	for(j=0;j<n;++j){
		y[j]=best.y[j];
		for(k=0;k<n;++k){
			x[j][k]=best.x[j][k];
		}
	}
	
	for(j=0;j<medians;++j){
		meds[j] = best.sol[j];
	}
	//designar(costu, medians,n,dist,cap,x,meds);
	lsheur2(costu,n,dist,medians,cap,x,y,10);

	free(dens);
	free(phero);
	free(cids);
	free(meds);
	free(sol);
	free(dtal);
	for(i=0;i<n_ants;++i){
		free(ant[i].sol);
		free(ant[i].y);
		for(j=0;j<n;++j){
			free(ant[i].x[j]);
		}
		free(ant[i].x);
	}
	free(ant);
	free(pa);
	free(best.sol);
	free(best.y);
	for(i=0;i<n;++i)
		free(best.x[i]);
	free(best.x);

}

void ant_copy(ants *ant,int **x,int *y,int n){

	int i,j;

	for(i=0;i<n;++i){
		ant->y[i]=y[i];
		for(j=0;j<n;++j)
			ant->x[i][j]=x[i][j];
	}
}

void density(costumer *costu, int n, double **dist, int cap, double *dens, int *cids){

	int i,j,k,soma_cids;
	int *ord, caps,tot_cids;
	double soma_dens;

	ord = (int *)malloc(sizeof(int)*n);
	for(i=0;i<n;++i){

		
		/* ordena por ordem de distancia do ponto "i" */
		for(j=0;j<n;++j) ord[j]=j;
		for(j=0;j<n-1;++j){
			for(k=j+1;k<n;++k){
				if(dist[ord[k]][i] < dist[ord[j]][i]){
					ord[k]^=ord[j]^=ord[k]^=ord[j];
				}
			}
		}

		/* soma a distancia até "encher" */
		caps = cap;
		dens[i] = 0;
		tot_cids=0;
		for(j=0;j<n;++j){
			if(costu[ord[j]].demand <= caps){
					dens[i] += dist[ord[j]][i];
					caps -= costu[ord[j]].demand;
					++tot_cids;
			}
		}
		cids[i] = tot_cids;
		
	}
	
	soma_cids=0;
	soma_dens=0.0;
	for(j=0;j<n;++j){
		soma_cids+=cids[j];
		soma_dens+=dens[j];
	}

	for(j=0;j<n;++j){
		dens[j]=((double)cids[j]/soma_cids)/(dens[j]/soma_dens);
	}


	free(ord);
}
