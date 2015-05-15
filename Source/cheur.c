#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "pmedian.h"

#define MAXPILHA 1000

void partition(designa *assig, costumer *costu, double **dist, int *meds, int lb, int ub, int *pj);

void cheur(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y,int max_it){

	int i,j,k,mudou;
	int maior_i,maior_j,temp,menor_med;
	double maior,menor,c;
	int *meds, *old_meds;

	
	meds = (int *)malloc(sizeof(int)*medians);
	old_meds = (int *)malloc(sizeof(int)*medians);
	

	/* Zera os valores de x e y */
	for(i=0;i<n;++i) for(j=0;j<n;++j) x[i][j]=0;
	for(i=0;i<n;++i) y[i]=0;
	
	/* gera solucao incial */

	maior = 0;
	maior_i = 0;
	maior_j = 0;
	
	/* Procura os dois pontos mais distantes */
	for(i=0;i<n;++i){
		for(j=i+1;j<n;++j){
			if(dist[i][j] > maior){
				maior_i=i;
				maior_j=j;
				maior=dist[i][j];
			}
		}
	}
	meds[0] = maior_i;
	meds[1] = maior_j;
	y[meds[0]] = 1;
	y[meds[1]] = 1;

	//i=2;
	/* Para as outras medianas encontre o ponto j tal que prod(dist[j][k]) seja maximo */
	for(i=2;i<medians;++i){
		maior=0;
		for(j=0;j<n;++j){
			c=1;
			for(k=0;k<i;++k){
				c*=dist[j][meds[k]];//*costu[j].demand;
				if(j==meds[k])c=0; 
			}
			if(c>maior){
				maior_i=j;
				maior=c;
			}
		}
		meds[i]=maior_i;
		y[meds[i]]=1;
	}

	/* desginar candidatos */
	designar(costu, medians,n,dist,cap,x,meds);

	for(k=0;k<max_it;++k){
		/* realoca medianas */
		for(i=0;i<n;++i) y[i]=0;
		for(i=0;i<medians;++i) old_meds[i] = meds[i];
		
		re_medians(n,medians,meds,x,y,dist);

		mudou=0;
		for(i=0;i<medians;++i) if(old_meds[i]!=meds[i]) mudou=1;

		if(!mudou) break;
	
	
		/* desginar candidatos de novo */
		for(i=0;i<n;++i) for(j=0;j<n;++j) x[i][j]=0;
		designar(costu, medians,n,dist,cap,x,meds);
	}

	printf("its necessarias passo 1 (const): %d\n",k);
		
	free(meds);
	free(old_meds);
}

/* Realoca as medianas calculando um novo ponto que minimize as distancias dos elementos do grupo */
void re_medians(int n, int medians,int *meds, int **x, int *y, double **dist){
	
	int i,j,k,menor_med;
	double menor,c;
	double media_x,media_y,media,a,b,dista;
	

	/* para cada mediana */
	for(i=0;i<medians;++i){
		menor=9e+10;
		menor_med = -1;
		/* procura o vértice... */
		for(j=0;j<n;++j){
			c = 0;
			/* ...pertencente a mediana atual... */
			if(x[j][meds[i]]){
				/* que minimiza a distancia total do cluster */
				for(k=0;k<n;++k){
					if(x[k][meds[i]]){
						c+=dist[j][k];
					}
				}
				if(c < menor){
					menor_med = j;
					menor=c;
				}
			}
		}
		/* se a mediana foi alterada, atribui o novo valor */
		meds[i]=menor_med;
		y[meds[i]]=1;
	}
}

/* Funcao para designar os vertices as medianas */
void designar(costumer *costu, int medians, int n, double **dist, int cap,int **x, int *meds){

	int temp,i,j,med_menor;
	double menor;
	designa * assig;
	int *caps;
	
	assig = (designa *)malloc(sizeof(designa)*n);
	caps = (int *)malloc(sizeof(int)*medians);

	/* atribui a capacidade de cada mediana */
	for(i=0;i<medians;++i) caps[i]=cap;
	
	/* para cada vértice, verifica qual a mediana mais próxima e relaciona na estrutura assig */
	for(i=0;i<n;++i){
		assig[i].custom = i;
		menor = 9e+10;
		for(j=0;j<medians;++j){
			if(dist[i][meds[j]] < menor){
				menor = dist[i][meds[j]];
				assig[i].median = j;
			}
		}
	}

	/* Ordena os vértices */
	sort(assig,n,dist,costu,meds,medians);

	/* para cada vértice... */
	for(i=0;i<n;++i){
		/* ...se couber na mediana mais proxima atribui a ela e seta seu valor em x...*/
		if(costu[assig[i].custom].demand <= caps[assig[i].median]){
			caps[assig[i].median]-=costu[assig[i].custom].demand;
			x[assig[i].custom][meds[assig[i].median]] = 1;
		/* ...senão encontra a próxima mediana menor que tenha capacidade suficiente para atender sua demanda */
		}else{
			menor = 9e+10;
			med_menor=-1;
			for(j=0;j<medians;++j){
				if(costu[assig[i].custom].demand <= caps[j] && dist[assig[i].custom][meds[j]]<menor){
					menor=dist[assig[i].custom][meds[j]];
					med_menor = j;
				}

			}
			if(med_menor != -1){
				caps[med_menor]-=costu[assig[i].custom].demand;
				x[assig[i].custom][meds[med_menor]] = 1;
			}
		}		
	}
	free(assig);
	free(caps);

}

/* Funcao de ordenacao dos vértices */

void sort(designa *assig,int n,double **dist,costumer *costu,int *meds,int medians){

	int i,j;

	struct bndtype {
		int lb;
		int ub;
	}newbnds;

	struct {
		int top;
		struct bndtype bounds[MAXPILHA];
	}stack;

	stack.top = -1;
	newbnds.lb = 0;
	newbnds.ub = n-1;
	stack.top++;
	stack.bounds[stack.top].lb = newbnds.lb;
	stack.bounds[stack.top].ub = newbnds.ub;	

	while(stack.top >= 0){
		newbnds.lb = stack.bounds[stack.top].lb;
		newbnds.ub = stack.bounds[stack.top].ub;
		stack.top--;
		while(newbnds.ub > newbnds.lb){
			partition(assig,costu, dist, meds, newbnds.lb, newbnds.ub, &j);
			if(j-newbnds.lb > newbnds.ub-j){
				i=newbnds.ub;
				newbnds.ub = j-1;

				stack.top++;
				stack.bounds[stack.top].lb = newbnds.lb;
				stack.bounds[stack.top].ub = newbnds.ub;	
				
				newbnds.lb = j+1;
				newbnds.ub = i;
			}else{
				i = newbnds.lb;
				newbnds.lb = j+1;
				stack.top++;
				stack.bounds[stack.top].lb = newbnds.lb;
				stack.bounds[stack.top].ub = newbnds.ub;	
				newbnds.lb = i;
				newbnds.ub = j-1;
			}
		}
	}
}

void partition(designa *assig, costumer *costu, double **dist, int *meds, int lb, int ub, int *pj){

	int down, temp, up;
	double a;

	a = dist[assig[lb].custom][meds[assig[lb].median]] + ((double)costu[assig[lb].custom].demand)/100000.0;
	up = ub;
	down = lb;
	while(down < up){
		while(((dist[assig[down].custom][meds[assig[down].median]] + ((double)costu[assig[down].custom].demand)/100000.0) <= a) && down < ub)
			down++;
		while((dist[assig[up].custom][meds[assig[up].median]] + ((double)costu[assig[up].custom].demand)/100000.0) > a && up > lb)
			up--;
		if(down < up){
				temp = assig[down].custom;
				assig[down].custom = assig[up].custom;
				assig[up].custom = temp;

				temp = assig[down].median;
				assig[down].median = assig[up].median;
				assig[up].median = temp;
			
		}
	}
	temp = assig[lb].custom;
	assig[lb].custom = assig[up].custom;
	assig[up].custom = temp;

	temp = assig[lb].median;
	assig[lb].median = assig[up].median;
	assig[up].median = temp;
	*pj = up;
}
