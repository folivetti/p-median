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

void density(costumer *costu, int n, double **dist, int cap, double *dens, int *cids,int medians, int **x, int *y);
void densp(costumer *costu, int n, double **dist, int cap, double *dens, int *cids, int medians, int **x, int *y);
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

	double *dens,soma_dens,soma;
	double *phero, *dtal, ro=0.8;
	double p,rand_p;
	double tal_max, tal_min;
	int i,j,k,it,n_ants=10;
	int *cids, *sol;
	ants *ant, best;
	int melhor_i;
	double melhor;
	double alpha=0, beta=1, soma_pher, *a,soma_a,m_dist;
	int *meds, its_n=0;

	dens = (double *)malloc(sizeof(double)*n);
	cids = (int *)malloc(sizeof(int)*n);
	phero = (double *)malloc(sizeof(double)*n);
	meds =	(int *)malloc(sizeof(int)*medians);
	sol =	(int *)malloc(sizeof(int)*n);
	dtal =	(double *)malloc(sizeof(double)*n);
	a =	(double *)malloc(sizeof(double)*n);
	ant = (ants *)malloc(sizeof(ants)*n_ants);

	m_dist=0;
	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			if(m_dist<dist[i][j])
				m_dist=dist[i][j];

	tal_max = (1/(1-ro))*(1/(n*m_dist));
	tal_min = tal_max/(2*n);

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
	
	//densp(costu, n, dist, cap, dens,cids,medians,x,y);
	density(costu, n, dist, cap, dens,cids,medians,x,y);

	soma_dens=0;
	for(i=0;i<n;++i)
		soma_dens+=dens[i];

	for(i=0;i<n;++i)
		phero[i]=tal_max;

	for(it=0;it<max_it;++it){

		if(its_n==50){
			for(i=0;i<n;++i)
				phero[i]=tal_max;
			its_n=0;
		}
			
		soma_pher=0.0;
		for(i=0;i<n;++i)
			soma_pher += pow(dens[i],beta)*pow(phero[i],alpha);

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
					a[k] = pow(dens[k],beta)*pow(phero[k],alpha)/soma_pher;
					if(!sol[k])
						soma_a += a[k];
				}
				rand_p = (double)rand()/RAND_MAX;
				for(k=0,p=0.0;k<n && p<rand_p;++k){
					if(sol[k]) continue;
					//printf("%f\n",pow(dens[k],alpha)*pow(phero[k],beta)/soma_pher);
					p += a[k]/soma_a;
				}
				--k;
				k=k%n;
				while(sol[k] && k>0) --k;
				//printf("ant[%d].sol[%d]=%d  %f %f\n",i,j,k,rand_p,p);
				ant[i].sol[j] = k;
				y[k]=1;
				meds[j]=k;
				sol[k]=1;
				soma_pher -= pow(dens[k],beta)*pow(phero[k],alpha);
			}
			/* avalia cada solucao */
			designar(costu, medians,n,dist,cap,x,meds);
			//lsheur2(costu,n,dist,medians,cap,x,y,1);
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
		melhor_i = 0;
		for(i=1;i<n_ants;++i){
			if(ant[i].fitness < ant[melhor_i].fitness){
				melhor = ant[i].fitness;
				melhor_i = i;
			}			
		}
		
		/* zera delta_tal */
		for(i=0;i<n;++i)
			dtal[i] = 0;
		
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
			printf("%d melhorou: %f\n",it,best.fitness);
			its_n=0;
		}else{
			++its_n;
		}

		/* calula delta_tal */
		for(i=0;i<medians;++i){
			dtal[best.sol[i]] = 1/best.fitness;
		}

		/* atualiza trilha de feromonio */
		for(i=0;i<n;++i){
			phero[i] = ro*phero[i] + dtal[i];
			if(phero[i] < tal_min)
				phero[i] = tal_min;
			else if(phero[i] > tal_max)
				phero[i] = tal_max;
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
	//lsheur2(costu,n,dist,medians,cap,x,y,10);

	free(dens);
	free(phero);
	free(cids);
	free(meds);
	free(sol);
	free(dtal);
	for(i=0;i<n_ants;++i){
		free(ant[i].sol);
		free(ant[i].y);
		for(j=0;j<n;++j)
			free(ant[i].x[j]);
		free(ant[i].x);
	}
	free(ant);
	free(a);
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

void density(costumer *costu, int n, double **dist, int cap, double *dens, int *cids, int medians, int **x, int *y){

	int i,j,k,soma_cids;
	int *ord, *meds, caps,soma_dens;

	ord = (int *)malloc(sizeof(int)*n);
	meds = (int *)malloc(sizeof(int)*medians);

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
		for(j=0;j<n;++j){
			if(costu[ord[j]].demand <= caps){
				dens[i] += dist[ord[j]][i];
				caps -= costu[ord[j]].demand;
			}else{
				break;
			}
		}
		//printf("dens[%d] = %d/%f = %f\n",i,j,dens[i],j/dens[i]);
		cids[i] = j;
		dens[i] = dens[i];	
		
	}
	
	soma_cids=0;
	soma_dens=0.0;
	for(j=0;j<n;++j){
		soma_cids+=cids[j];
		soma_dens+=dens[j];
	}

	for(j=0;j<n;++j){
//		printf("%d: %f %f %f\n",j,(double)cids[ord[j]]/soma_cids,dens[ord[j]]/soma_dens,((double)cids[ord[j]]/soma_cids)/(dens[ord[j]]/soma_dens));
		dens[j]=((double)cids[ord[j]]/soma_cids)/(dens[ord[j]]/soma_dens);
	}

/*
	for(j=0;j<n;++j) ord[j]=j;
	for(j=0;j<n-1;++j){
		for(k=j+1;k<n;++k){
			if(dens[ord[k]] < dens[ord[j]]){
				ord[k]^=ord[j]^=ord[k]^=ord[j];
			}
		}
	}
	
for(i=0;i<n;++i) y[i]=0;
	for(i=0;i<medians;++i){
		y[ord[i]]=1;
		meds[i]=ord[i];
	}
	
	for(i=0;i<n;++i) for(j=0;j<n;++j) x[i][j]=0;
	designar(costu, medians,n,dist,cap,x,meds);
*/
	free(ord);
	free(meds);		
}


/* algoritmo construtivo utilizando densidade */

void densp(costumer *costu, int n, double **dist, int cap, double *dens, int *cids, int medians, int **x, int *y){

	int i,j,k,soma_cids,it,tot_cids;
	int *ord, *meds,caps;
	double soma_dens;
	int **t_x, *s_x;

	ord = (int *)malloc(sizeof(int)*n);
	s_x = (int *)malloc(sizeof(int)*n);
	t_x = (int **)malloc(sizeof(int *)*n);
	for(i=0;i<n;++i)
		t_x[i] = (int *)malloc(sizeof(int)*n);
	meds = (int *)malloc(sizeof(int)*medians);
	
	for(i=0;i<n;++i) y[i]=0;
	for(i=0;i<n;++i) s_x[i]=0;
	for(i=0;i<n;++i) for(j=0;j<n;++j) x[i][j]=0;
	for(i=0;i<n;++i) for(j=0;j<n;++j) t_x[i][j]=0;
	
for(it=0;it<medians;++it){
	
	for(i=0;i<n;++i){

		if(y[i]){
			dens[i]=1;
			cids[i]=0;
			continue;
		}
		
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
				if(!s_x[ord[j]]){
					dens[i] += dist[ord[j]][i];
					caps -= costu[ord[j]].demand;
					t_x[ord[j]][i] = 1;
					++tot_cids;
				}
			}else{
				break;
			}
		}
		//printf("dens[%d] = %d/%f = %f\n",i,j,dens[i],j/dens[i]);		
		cids[i] = tot_cids;
		dens[i] = dens[i];	
		
	}
	//printf("\n");
	
	soma_cids=0;
	soma_dens=0.0;
	for(j=0;j<n;++j){
		if(!y[j]){
		soma_cids+=cids[j];
		soma_dens+=dens[j];
		}
	//	printf("soma_cids: %d\tsoma_dens: %f\n",soma_cids,soma_dens);
	}

	for(j=0;j<n;++j){
		dens[j]=((double)cids[j]/soma_cids)/(dens[j]/soma_dens);
	//	printf("%d: %f\n",j,dens[j]);
	}


	for(j=0;j<n;++j) ord[j]=j;
	for(j=0;j<n-1;++j){
		for(k=j+1;k<n;++k){
			if(dens[ord[k]] > dens[ord[j]]){
				ord[k]^=ord[j]^=ord[k]^=ord[j];
			}
		}
	}

	y[ord[0]]=1;
	for(k=0;k<n;++k){
		x[k][ord[0]]=t_x[k][ord[0]];
		if(x[k][ord[0]])
			s_x[k] = 1;
	}

	for(j=0;j<n;++j)
		for(k=0;k<n;++k)
			t_x[j][k]=0;

	}
	
	for(i=0,j=0;i<n;++i){
		if(y[i])
			meds[j++]=i;
	}
	
	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			x[i][j]=0;
	designar(costu, medians,n,dist,cap,x,meds);
	
	free(ord);
	free(s_x);
	for(j=0;j<n;++j)
		free(t_x[j]);
	free(t_x);
	free(meds);		
}
