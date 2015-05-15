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
	double tal_max, tal_min,cf=0.0;
	int i,j,k,it,n_ants=10,bs_update=0;
	int *cids, *sol;
	ants *ant, best, best_r;
	int melhor_i;
	double melhor,soma_fit;
	double alpha=3, beta=1, soma_pher, *pa,soma_a,m_dist;
	int *meds;

	dens = (double *)malloc(sizeof(double)*n);
	cids = (int *)malloc(sizeof(int)*n);
	phero = (double *)malloc(sizeof(double)*n);
	meds =	(int *)malloc(sizeof(int)*medians);
	sol =	(int *)malloc(sizeof(int)*n);
	dtal =	(double *)malloc(sizeof(double)*n);
	pa =	(double *)malloc(sizeof(double)*n);
	ant = (ants *)malloc(sizeof(ants)*n_ants);


	tal_max = 0.999; 
	tal_min = 0.001;

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
	
	best_r.sol = (int *)malloc(sizeof(int)*medians);
	best_r.y = (int *) malloc(sizeof(int)*n);
	best_r.x = (int **) malloc(sizeof(int *)*n);
	for(j=0;j<n;++j){
		best_r.x[j] = (int *) malloc(sizeof(int)*n);
	}

	best_r.fitness = 9e+10;

	density(costu, n, dist, cap, dens,cids);


	for(i=0;i<n;++i)
		phero[i]=0.5;

	for(it=0;it<max_it;++it){

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
				soma_pher -= pow(dens[k],beta)*pow(phero[k],alpha);
			}
			
			/* avalia cada solucao */
			designar(costu, medians,n,dist,cap,x,meds);
			
			lsheur2(costu,n,dist,medians,cap,x,y,1);
			ant_copy(&ant[i],x,y,n);
			
			ant[i].fitness = calculate(costu,n,dist,medians,cap,x,y);
			if(!factivel(costu,n,medians,cap,x,y)){
				ant[i].fitness *= 1.5;
			}
			
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
		melhor_i = 0;
		for(i=1;i<n_ants;++i){
			if(ant[i].fitness < ant[melhor_i].fitness){
				melhor = ant[i].fitness;
				melhor_i = i;
			}			
		}

		
		
		/* atualiza melhor solucao global e local */
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
			printf("%d melhorou: %f\n",it,best.fitness);
		}

		if(ant[melhor_i].fitness < best_r.fitness){
			for(i=0;i<medians;++i){
				best_r.sol[i] = ant[melhor_i].sol[i];
			}
			for(i=0;i<n;++i){
				best_r.y[i] = ant[melhor_i].y[i];
				for(j=0;j<n;++j){
					best_r.x[i][j] = ant[melhor_i].x[i][j];
				}
			}
			best_r.fitness = ant[melhor_i].fitness;
		}

		for(i=0,soma_fit=0.0;i<n_ants;++i){
			soma_fit+=ant[i].fitness;
		}
		soma_fit += (best.fitness + best_r.fitness);
		
		/* zera delta_tal */
		for(i=0;i<n;++i)
			dtal[i] = 0.0;
		
		/* calula delta_tal */
		if(bs_update==0){
			if(cf < 0.4){
				for(i=0;i<medians;++i){
					dtal[ant[melhor_i].sol[i]] += 1.0*ant[melhor_i].fitness/soma_fit;
				}
			}else if(cf >= 0.4 && cf < 0.6){
				for(i=0;i<medians;++i){
					dtal[ant[melhor_i].sol[i]] += (2.0/3.0)*ant[melhor_i].fitness/soma_fit;
					dtal[best_r.sol[i]] += (1.0/3.0)*best_r.fitness/soma_fit;
				}
			}else if(cf >=0.6 && cf < 0.8){
				for(i=0;i<medians;++i){
					dtal[ant[melhor_i].sol[i]] += (1.0/3.0)*ant[melhor_i].fitness/soma_fit;
					dtal[best_r.sol[i]] += (2.0/3.0)*best_r.fitness/soma_fit;
				}
			}else if(cf >= 0.8){
				for(i=0;i<medians;++i){
					dtal[best_r.sol[i]] += 0.7*best_r.fitness/soma_fit;
					dtal[best.sol[i]] += 0.3*best.fitness/soma_fit;
				}
			}
		}else{
			for(i=0;i<medians;++i){
				dtal[best.sol[i]] += 1.0*best.fitness/soma_fit;
			}
		}

		/* atualiza trilha de feromonio */
		for(i=0;i<n;++i){
			phero[i] = phero[i] + ro*(dtal[i]-phero[i]);
			if(phero[i] < tal_min)
				phero[i] = tal_min;
			else if(phero[i] > tal_max)
				phero[i] = tal_max;
		}
		
		

		for(i=0,cf=0.0;i<n;++i){
			if((tal_max-phero[i])>(phero[i]-tal_min))
				cf += tal_max-phero[i];
			else
				cf += phero[i]-tal_min;
		}
		cf = cf/(n*(tal_max-tal_min));
		cf = 2*(cf-0.5);
		//printf("cf = %f\n",cf);

		if(cf>0.90){
			if(bs_update == 1){
				for(i=0;i<n;++i)
					phero[i]=0.5;
				best_r.fitness = 9e+10;
				bs_update = 0;
			}else{
				bs_update = 1;
				printf("update now: %d\n",it);
			}
		}
		
		
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
