#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "pmedian.h"
#include <time.h>
#include <math.h>

typedef struct ag_pmc{

	int *cromo;
	double fitness;
	int **x;
	int *y;
	
}ag_pmc;

void gen_pop(ag_pmc *pop, int pops, int n, int medians, double **dist);
void eval_fit(ag_pmc *pop,int pops,int n,int medians,int **x,int *y,double **dist,int cap,costumer *costu);
void sort_pop(ag_pmc *pop,int pops,int medians,int n);
void hyper(ag_pmc *pop, int pops, double ph, int n, int medians, int **x, int *y, double **dist, int cap, costumer *costu,int max_it);
void hypertroca(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds,double answer);
void ag_select(int *sels,int pais,ag_pmc *pop,int pops);
void tournament(int *sels,int pais,ag_pmc *pop,int pops);
int bern(double per);
void mutation(int *sels,int pais,ag_pmc *pop,int pops,int n,int medians, double **dist, double pm);
void crossover(int *sels,int pais,ag_pmc *filho,ag_pmc *pop,int pops,int n,int medians,double pc);
void reproduce(ag_pmc *filho,int pais,ag_pmc *pop,int pops,int medians,int n);
int ag_fact(ag_pmc ind,int **x, int *y, int n, int cap, int medians, costumer *costu, double **dist);

/*
 * Variaveis:
 *
 * pc 	- probabilidade de crossover
 * pm 	- probabilidade de mutacao
 * ph 	- probabilidade de hypermutacao
 * pops	- no. de individuos
 * 
 */
void ag(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it){

	int i,j,pops=50,it;
	ag_pmc *pop, *filho,best;
	double ph=0.2;
	double kp = 5;
	int pais=20;
	int *sels;
	double pm=0.1;
	double pc=1;
	double passo_m = (0.2-0.01)/max_it;
	double passo_c = (1-0.5)/max_it;

	pops = floor(kp*(double)n/medians);
	pais = 2;//10*pops/100;
	if(pais%2)
		++pais;
	
	pop = (ag_pmc *)malloc(sizeof(ag_pmc)*pops);
	for(i=0;i<pops;++i){
		pop[i].cromo = (int *)malloc(sizeof(int)*medians);
		pop[i].y = (int *)malloc(sizeof(int)*n);
		pop[i].x = (int **)malloc(sizeof(int *)*n);
		for(j=0;j<n;++j)
			pop[i].x[j] = (int *)malloc(sizeof(int)*n);
	}
	filho = (ag_pmc *)malloc(sizeof(ag_pmc)*pais);
	for(i=0;i<pais;++i){
		filho[i].cromo = (int *)malloc(sizeof(int)*medians);
		filho[i].y = (int *)malloc(sizeof(int)*n);
		filho[i].x = (int **)malloc(sizeof(int *)*n);
		for(j=0;j<n;++j)
			filho[i].x[j] = (int *)malloc(sizeof(int)*n);
	}
	best.cromo = (int *)malloc(sizeof(int)*medians);
	best.y = (int *)malloc(sizeof(int)*n);
	best.x = (int **)malloc(sizeof(int *)*n);
	for(j=0;j<n;++j)
		best.x[j] = (int *)malloc(sizeof(int)*n);
	sels = (int *)malloc(sizeof(int)*pais);

	/* Gera solucoes aleatorias iniciais */
	gen_pop(pop,pops,n,medians,dist);

	/* Avalia Fitness */
	eval_fit(pop,pops,n,medians,x,y,dist,cap,costu);
	
	/* ordena a populacao */
	sort_pop(pop,pops,medians,n);

	/* Realiza Busca local em uma porcentagem da populacao */
//	hyper(pop,pops,ph,n,medians,x,y,dist,cap,costu);

	/* ordena a populacao */ 
//	sort_pop(pop,pops,medians);

	best.fitness = pop[0].fitness;
	for(i=0;i<medians;++i)
		best.cromo[i]=pop[0].cromo[i];
	for(i=0;i<n;++i){
		best.y[i] = pop[0].y[i];
		for(j=0;j<n;++j){
			best.x[i][j] = pop[0].x[i][j];
		}
	}

	for(it=0;it<max_it;++it){

		/* crossover & mutacao */
		ag_select(sels,pais,pop,pops);

		mutation(sels,pais,pop,pops,n,medians,dist,pm);
		eval_fit(pop,pops,n,medians,x,y,dist,cap,costu);
		sort_pop(pop,pops,medians,n);

		if(pop[0].fitness < best.fitness && ag_fact(pop[0],x, y, n, cap, medians, costu,dist)){
			best.fitness = pop[0].fitness;
			for(i=0;i<medians;++i)
				best.cromo[i]=pop[0].cromo[i];
			printf("melhorou: %d %f\n",it,best.fitness);
			for(i=0;i<n;++i){
				best.y[i] = pop[0].y[i];
				for(j=0;j<n;++j){
					best.x[i][j] = pop[0].x[i][j];
				}
			}
		}

		crossover(sels,pais,filho,pop,pops,n,medians,pc);
		eval_fit(filho,pais,n,medians,x,y,dist,cap,costu);

		reproduce(filho,pais,pop,pops,medians,n);
		sort_pop(pop,pops,medians,n);

		//hyper(pop,pops,ph,n,medians,x,y,dist,cap,costu);
		
		if((it%20)==0){
			hyper(pop,pops,ph,n,medians,x,y,dist,cap,costu,1);
			sort_pop(pop,pops,medians,n);
		}
		

		if(pop[0].fitness < best.fitness && ag_fact(pop[0],x, y, n, cap, medians, costu,dist)){
			best.fitness = pop[0].fitness;
			for(i=0;i<medians;++i)
				best.cromo[i]=pop[0].cromo[i];
			printf("melhorou: %d %f\n",it,best.fitness);
			for(i=0;i<n;++i){
				best.y[i] = pop[0].y[i];
				for(j=0;j<n;++j){
					best.x[i][j] = pop[0].x[i][j];
				}
			}
		}
		pm+=passo_m;
		pc-=passo_c;
	}

	hyper(pop,pops,1,n,medians,x,y,dist,cap,costu,1);
	sort_pop(pop,pops,medians,n);

		if(pop[0].fitness < best.fitness && ag_fact(pop[0],x, y, n, cap, medians, costu,dist)){
			best.fitness = pop[0].fitness;
			for(i=0;i<medians;++i)
				best.cromo[i]=pop[0].cromo[i];
			printf("melhorou: %d %f\n",it,best.fitness);
			for(i=0;i<n;++i){
				best.y[i] = pop[0].y[i];
				for(j=0;j<n;++j){
					best.x[i][j] = pop[0].x[i][j];
				}
			}
		}
	/* avalia primeiro individuo e coloca em x e y */
	/*
	printf("%f\n",best.fitness);
	eval_fit(&best,-1,n,medians,x,y,dist,cap,costu);
	printf("%f\n",best.fitness);
	*/
	hyper(&best,-1,1,n,medians,x,y,dist,cap,costu,15);
	printf("%f\n",best.fitness);
	//eval_fit(&best,-1,n,medians,x,y,dist,cap,costu);
	
	for(i=0;i<n;++i){
		y[i] = best.y[i];
		for(j=0;j<n;++j){
			x[i][j] = best.x[i][j];
		}
	}

	for(i=0;i<pops;++i){
		free(pop[i].cromo);
		free(pop[i].y);
		for(j=0;j<n;++j)
			free(pop[i].x[j]);
		free(pop[i].x);
	}
	free(pop);
	for(i=0;i<pais;++i){
		free(filho[i].cromo);
		free(filho[i].y);
		for(j=0;j<n;++j)
			free(filho[i].x[j]);
		free(filho[i].x);
	}
	free(best.cromo);
	free(filho);
	free(sels);
}

int ag_fact(ag_pmc ind,int **x, int *y, int n, int cap, int medians, costumer *costu, double **dist){

	int j,k,fact,*meds=NULL;

	meds = (int *)malloc(sizeof(int)*medians);
	
	for(j=0;j<n;++j){
		y[j]=0;
	for(k=0;k<n;++k)
		x[j][k]=0;
	}
	
	for(j=0;j<medians;++j){
		meds[j]=ind.cromo[j];
		y[meds[j]]=1;
	}
	designar(costu, medians, n, dist, cap,x, meds);
	fact = factivel(costu,n,medians,cap,x,y);
		
	free(meds);

	return fact;

}

/* gera individuos aleatorios */
void gen_pop(ag_pmc *pop, int pops, int n, int medians, double **dist){

	int *sorteio, rnd, i, j,k,l,temp;

	sorteio = (int *)malloc(sizeof(int)*n);

	for(i=0;i<n;++i){
		sorteio[i]=i;
	}

	for(j=0;j<pops;j+=n/medians){
		
		for(i=0;i<n;++i){
			rnd = floor((double)rand()*n/RAND_MAX);
			rnd = rnd%n;
			temp = sorteio[i];
			sorteio[i]=sorteio[rnd];
			sorteio[rnd]=temp;
		}

		for(i=j,l=0;i<j+(n/medians) && i<pops;++i,l+=medians){
			for(k=0;k<medians;++k)
				pop[i].cromo[k]=sorteio[k+l];
		}
	}

	free(sorteio);
}
	
void eval_fit(ag_pmc *pop,int pops,int n,int medians,int **x,int *y,double **dist,int cap,costumer *costu){

	int i,j,k;
	static double m_dist;
	int soma;
	int *meds;

	/* retirar isso daqui pra aumentar desempenho */
		for(i=0;i<n;++i)
			for(j=0;j<n;++j)
				if(m_dist<dist[i][j])
					m_dist=dist[i][j];
	
	meds = (int *)malloc(sizeof(int)*medians);
if(pops>-1){	
	for(i=0;i<pops;++i){
		for(j=0;j<n;++j){
			y[j]=0;
			for(k=0;k<n;++k)
				x[j][k]=0;
		}
		
		for(j=0;j<medians;++j){
			meds[j]=pop[i].cromo[j];
			y[meds[j]]=1;
		}
		
		designar(costu, medians, n, dist, cap,x, meds);
		pop[i].fitness = calculate(costu,n,dist,medians,cap,x,y);
		soma = 0;
		for(j=0;j<n;++j)
			for(k=0;k<n;++k)
				soma+=x[j][k];
		soma = n-soma;
		pop[i].fitness += soma*m_dist;

		for(j=0;j<n;++j){
			pop[i].y[j] = y[j];
			for(k=0;k<n;++k){
				pop[i].x[j][k]=x[j][k];
			}
		}
		//printf("%d: %f\n",i,pop[i].fitness);
	}
}else{
		for(j=0;j<n;++j){
			y[j]=0;
			for(k=0;k<n;++k)
				x[j][k]=0;
		}
		
		for(j=0;j<medians;++j){
			meds[j]=pop->cromo[j];
			y[meds[j]]=1;
		}
		
		designar(costu, medians, n, dist, cap,x, meds);
		pop->fitness = calculate(costu,n,dist,medians,cap,x,y);
		for(j=0;j<n;++j){
			pop->y[j] = y[j];
			for(k=0;k<n;++k){
				pop->x[j][k]=x[j][k];
			}
		}
}

	free(meds);
}
	
void sort_pop(ag_pmc *pop,int pops,int medians,int n){

	int c_temp, i,j,k,l;
	double f_temp;

	for(i=0;i<pops-1;++i){

		for(j=i+1;j<pops;++j){
			
			if(pop[i].fitness > pop[j].fitness){

				f_temp = pop[i].fitness;
				pop[i].fitness = pop[j].fitness;
				pop[j].fitness = f_temp;

				for(k=0;k<medians;++k){
					c_temp = pop[i].cromo[k];
					pop[i].cromo[k] = pop[j].cromo[k];
					pop[j].cromo[k] = c_temp;
				}
				for(l=0;l<n;++l){
					c_temp = pop[i].y[l];
					pop[i].y[l] = pop[j].y[l];
					pop[j].y[l] = c_temp;
					for(k=0;k<n;++k){
						c_temp = pop[i].x[l][k];
						pop[i].x[l][k] = pop[j].x[l][k];
						pop[j].x[l][k] = c_temp;
					}
				}
			}
		}
	}

	/*
	for(i=0;i<pops;++i)
		printf("%f\t",pop[i].fitness);
	printf("\n");
	*/
}


void hyper(ag_pmc *pop, int pops, double ph, int n, int medians, int **x, int *y, double **dist, int cap, costumer *costu,int max_it){

	int i,j,k,soma;
	int *meds;
	static double m_dist;
	
	/* retirar isso daqui pra aumentar desempenho */
		for(i=0;i<n;++i)
			for(j=0;j<n;++j)
				if(m_dist<dist[i][j])
					m_dist=dist[i][j];

	meds = (int *)malloc(sizeof(int)*medians);
if(pops>-1){	
	for(i=0;i<pops;++i){

		if((double)rand()/RAND_MAX <= ph){
			for(j=0;j<n;++j){
				y[j]=0;
				for(k=0;k<n;++k)
					x[j][k]=0;
			}

			for(j=0;j<medians;++j){
				meds[j]=pop[i].cromo[j];
				y[meds[j]]=1;
			}
		
			designar(costu, medians, n, dist, cap,x, meds);
			//hypertroca(costu, n, dist, medians, cap, x, y,meds,pop[i].fitness);
			lsheur2(costu, n, dist, medians, cap, x, y, max_it);
			pop[i].fitness = calculate(costu,n,dist,medians,cap,x,y);

			soma = 0;
			for(j=0;j<n;++j)
				for(k=0;k<n;++k)
					soma+=x[j][k];
			soma = n-soma;
			pop[i].fitness += soma*m_dist;
	
			for(j=0,k=0;j<n && k<medians;++j)
				if(y[j])
					pop[i].cromo[k++]=j;
			//printf("%f\n",pop[i].fitness);
			for(j=0;j<n;++j){
				pop[i].y[j] = y[j];
				for(k=0;k<n;++k){
					pop[i].x[j][k]=x[j][k];
				}
			}
		}
	}
}else{
			for(j=0;j<n;++j){
				y[j]=pop->y[j];
				for(k=0;k<n;++k)
					x[j][k]=pop->x[j][k];
			}

			for(j=0;j<medians;++j){
				meds[j]=pop->cromo[j];
				//y[meds[j]]=1;
			}
		
			//designar(costu, medians, n, dist, cap,x, meds);
			lsheur2(costu, n, dist, medians, cap, x, y, max_it);
			pop->fitness = calculate(costu,n,dist,medians,cap,x,y);
			for(j=0,k=0;j<n && k<medians;++j)
				if(y[j])
					pop->cromo[k++]=j;
	
			for(j=0;j<n;++j){
				pop->y[j] = y[j];
				for(k=0;k<n;++k){
					pop->x[j][k]=x[j][k];
				}
			}

}

	free(meds);

}

void hypertroca(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds,double answer){

	int i,j,k,melhor_i,melhor_fr,velho;
	int fr, *L1, *L2;
	double lucro,melhor=9e+10,new_ans;
	int it,max_it=1;

	L1 = (int *)malloc(sizeof(int)*n);
	L2 = (int *)malloc(sizeof(int)*n);

for(it=0;it<max_it;it++){	
	
	atua_l(L1, L2, n, meds, medians, dist);
	for(i=0;i<n;++i){

		if(!y[i]){
			melhor_subs(meds,i,L1,L2,n,medians,dist,&fr,&lucro);

			y[i]=1;
			y[meds[fr]]=0;
			velho=meds[fr];
			meds[fr]=i;

			for(j=0;j<n;++j) for(k=0;k<n;++k) x[j][k]=0;

			designar(costu, medians, n, dist, cap,x, meds);
			new_ans = calculate(costu,n,dist,medians,cap,x,y);
			if(new_ans < melhor){
				melhor = new_ans;
				melhor_i = i;
				melhor_fr = fr;
			}

			meds[fr]=velho;
			y[i]=0;
			y[meds[fr]]=1;
		}
	}
	if(melhor < answer){
		y[melhor_i]=1;
		y[meds[melhor_fr]]=0;
		meds[melhor_fr]=melhor_i;
		for(j=0;j<n;++j) for(k=0;k<n;++k) x[j][k]=0;
		designar(costu, medians, n, dist, cap,x, meds);
	}
}

	free(L1);
	free(L2);
}

void ag_select(int *sels,int pais,ag_pmc *pop,int pops){

	/* Tournament */
	tournament(sels,pais,pop,pops);

	/* SUS */
}

void tournament(int *sels,int pais,ag_pmc *pop,int pops){

	int i,sel1,sel2;

	for(i=0;i<pais;++i){
		if((i%2)==1){
			do{
				sel1 = floor((double)rand()*pops/RAND_MAX);
				sel1=sel1%pops;
				do{
					sel2 = floor((double)rand()*pops/RAND_MAX);
					sel2 = sel2%pops;
				}while(sel2==sel1);
				if(pop[sel1].fitness < pop[sel2].fitness){
					sels[i]=sel1;
				}else{
					sels[i]=sel2;
				}
			}while(sels[i]==sels[i-1]);
		}else{
			sel1 = floor((double)rand()*pops/RAND_MAX);
			sel1 = sel1%pops;
			do{
				sel2 = floor((double)rand()*pops/RAND_MAX);
				sel2 = sel2%pops;
			}while(sel2==sel1);
			if(pop[sel1].fitness < pop[sel2].fitness){
				sels[i]=sel1;
			}else{
				sels[i]=sel2;
			}			
		}
	}

}

void mutation(int *sels,int pais,ag_pmc *pop,int pops,int n,int medians,double **dist, double pm){

	int i,j,k,l,idx=0,menor_pto;
	int achou,rnd,temp;
	int *sorteio=NULL, *mask=NULL;
	double menor;

	sorteio = (int *)malloc(sizeof(int)*(n-medians));
	mask = (int *)malloc(sizeof(int)*medians);
	
	
	for(i=0;i<pais;++i){

		if((double)rand()/RAND_MAX <= pm){

			for(j=0,idx=0;j<n && idx<(n-medians);++j){
				achou=0;
				for(k=0;k<medians;++k){
					if(j==pop[sels[i]].cromo[k]){
						achou=1;
						break;
					}
				}
				if(!achou){
					sorteio[idx]=j;
					++idx;
				}
			}
			for(j=0;j<medians;++j){
				mask[j] = bern(0.2);
			}
	/*	
			for(j=0;j<n-medians;++j){
				rnd = floor((double)rand()*(n-medians)/RAND_MAX);
				rnd = rnd%(n-medians);
				temp = sorteio[j];
				sorteio[j]=sorteio[rnd];
				sorteio[rnd]=temp;
			}
		
			for(j=0;j<medians;++j){
				if(mask[j]){
					pop[sels[i]].cromo[j] = sorteio[j];
				}
			}
		*/

			
			for(j=0;j<medians;++j){
				if(mask[j]){
					menor=9e+10;
					menor_pto = pop[sels[i]].cromo[j];
					for(k=0;k<idx;++k){
						if(dist[pop[sels[i]].cromo[j]][sorteio[k]] < menor){
							achou=0;
							for(l=0;l<j;++l){
								if(pop[sels[i]].cromo[l]==sorteio[k]){
									achou=1;
								}
							}
							if(!achou){
								menor = dist[pop[sels[i]].cromo[j]][sorteio[k]];
								menor_pto = sorteio[k];
							}
						}
					}
//				 printf("troca %d por %d\n",pop[sels[i]].cromo[j],menor_pto);
					pop[sels[i]].cromo[j] = menor_pto;
				}
			}
			
		}
	}

	if(mask != NULL)
		free(mask);
	if(sorteio != NULL)
		free(sorteio);
}

void crossover(int *sels,int pais,ag_pmc *filho,ag_pmc *pop,int pops,int n,int medians,double pc){

	int i,k,l,m,j,o,p;
	int pai1,pai2,igual1,igual2;
	int *troca1,*troca2,*ntroca1,*ntroca2,*mask1,*mask2;

	troca1 = (int *)malloc(sizeof(int)*medians);
	troca2 = (int *)malloc(sizeof(int)*medians);
	ntroca1 = (int *)malloc(sizeof(int)*medians);
	ntroca2 = (int *)malloc(sizeof(int)*medians);
	mask1 = (int *)malloc(sizeof(int)*medians);
	mask2 = (int *)malloc(sizeof(int)*medians);

	for(k=0;k<pais;k+=2){

		pai1=sels[k];
		pai2=sels[k+1];
		l=0;
		m=0;
		o=0;
		p=0;
//printf("pai1: %d,pai2: %d\n",sels[k],sels[k+1]);		
		for(i=0;i<medians;++i){
			igual1=0;
			igual2=0;
			for(j=0;j<medians;++j){
				if(pop[pai1].cromo[i]==pop[pai2].cromo[j]){
					igual1=1;
				}
				if(pop[pai2].cromo[i]==pop[pai1].cromo[j]){
					igual2=1;
				}
			}
			if(!igual1)
				troca1[l++] = i;
			else
				ntroca1[o++] = i;
			if(!igual2)
				troca2[m++] = i;
			else
				ntroca2[p++] = i;

		}
		/*
printf("l=%d m=%d\n",l,m);
if(l!=m){
	for(j=0;j<medians;++j){
		printf("%d %d\n",pop[pai1].cromo[j],pop[pai2].cromo[j]);
	}
	
}
*/
		for(i=0;i<l;++i){
			mask1[i] = bern(0.5);
			mask2[i] = bern(0.5);
		}
		
		
		for(i=0;i<o;++i){
			filho[k].cromo[i]=pop[pai1].cromo[ntroca1[i]];
			filho[k+1].cromo[i]=pop[pai2].cromo[ntroca2[i]];
		}
		for(i=0;i<l;++i){
//			printf("i=%d %d %d %d\n",i,troca1[i],troca2[i],mask1[i]);
			if(mask1[i]){
				filho[k].cromo[o+i]=pop[pai2].cromo[troca2[i]];
			}else{
				filho[k].cromo[o+i]=pop[pai1].cromo[troca1[i]];
			}

			if(mask2[i]){
				filho[k+1].cromo[o+i]=pop[pai1].cromo[troca1[i]];
			}else{
				filho[k+1].cromo[o+i]=pop[pai2].cromo[troca2[i]];
			}			
		}
/*
		for(i=0;i<medians;++i)
			printf("%d ",filho[k].cromo[i]);
		printf("\n");
		for(i=0;i<medians;++i)
			printf("%d ",filho[k+1].cromo[i]);
		printf("\n");
*/
	}

	free(troca1);
	free(troca2);
	free(ntroca1);
	free(ntroca2);
	free(mask1);
	free(mask2);
}

void reproduce(ag_pmc *filho,int pais,ag_pmc *pop,int pops,int medians,int n){

	int *ord_fit,*ord_filho;
	int i,j,k,l,melhor,skip;

	ord_fit = (int *)malloc(sizeof(int)*pops);
	ord_filho = (int *)malloc(sizeof(int)*pais);
	
	for(i=0;i<pops;++i) ord_fit[i]=i;
	for(i=0;i<pais;++i) ord_filho[i]=i;

	/* Ordena do pior para o melhor */
	for(i=0;i<pops-1;++i){
		for(j=i+1;j<pops;++j){
			if(pop[ord_fit[i]].fitness < pop[ord_fit[j]].fitness){
				ord_fit[i]^=ord_fit[j]^=ord_fit[i]^=ord_fit[j];
			}
		}
	}
	
	for(i=0;i<pais-1;++i){
		for(j=i+1;j<pais;++j){
			if(filho[ord_filho[i]].fitness < filho[ord_filho[j]].fitness){
				ord_filho[i]^=ord_filho[j]^=ord_filho[i]^=ord_filho[j];
			}
		}
	}

	/* verifica a partir de qual filho sera copiado */
	melhor=pais;
	for(i=0;i<pais;++i){
		if(filho[ord_filho[i]].fitness < pop[ord_fit[0]].fitness){
			melhor = i;
			break;
		}
	}
	for(i=melhor,j=0;i<pais;++i){
		skip=0;
		for(k=0;k<pops;++k){
			if(pop[k].fitness == filho[ord_filho[i]].fitness){
				skip = 1;
				break;				
			}
		}
		if(skip) continue;
		
		pop[ord_fit[j]].fitness = filho[ord_filho[i]].fitness;
		for(k=0;k<medians;++k){
			pop[ord_fit[j]].cromo[k] = filho[ord_filho[i]].cromo[k];
		}
		for(k=0;k<n;++k){
			pop[ord_fit[j]].y[k] = filho[ord_filho[i]].y[k];
			for(l=0;l<n;++l){
				pop[ord_fit[j]].x[k][l]=filho[ord_filho[i]].x[k][l];
			}
		}
		++j;
	}
	
	free(ord_fit);
	free(ord_filho);
}

int bern(double per){

	double rnd;

	rnd = (double)rand()/RAND_MAX;

	return rnd<=per;
}
