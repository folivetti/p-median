#include <stdio.h>
#include <math.h>
#include "pmedian.h"



void cheur(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y,int max_it){

	int i,j,k;
	int maior_i,maior_j,temp,menor_med;
	double maior,menor,c;
	int *meds;
	int *caps;

	designa * assig;
	
	meds = (int *)malloc(sizeof(int)*medians);
	caps = (int *)malloc(sizeof(int)*medians);
	assig = (designa *)malloc(sizeof(designa)*n);
	

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

	i=2;
	/* Para as outras medianas encontre o ponto j tal que prod(dist[j][k]) seja maximo */
	for(;i<medians;++i){
		maior=0;
		for(j=0;j<n;++j){

			c=1;
			for(k=0;k<i;++k){
				c*=dist[j][meds[k]]*costu[j].demand;
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
	desginar(costu, assig, medians,n,dist,caps,cap,x,meds);
	//desginar_rate(costu, assig, medians, n, dist, caps,cap,x,meds);

	for(k=0;k<max_it;++k){
		/* realoca medianas */
		for(i=0;i<n;++i) y[i]=0;
		re_medians(n,medians,meds,x,y,dist);
	
	
		/* desginar candidatos de novo */
		for(i=0;i<n;++i) for(j=0;j<n;++j) x[i][j]=0;
		desginar(costu, assig, medians,n,dist,caps,cap,x,meds);
		//desginar_rate(costu, assig, medians, n, dist, caps,cap,x,meds);
	}
		
	free(meds);
	free(caps);
	free(assig);
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

	/*
	for(i=0;i<medians;++i){
		media_x=0;
		media_y=0;
		media = 0;
		for(j=0;j<n;++j){
			if(x[j][meds[i]]){
				media_x += costu[j].x*costu[j].demand;
				media_y += costu[j].y*costu[j].demand;
				media += costu[j].demand;
			}
		}
		media_x /= media;
		media_y /= media;
		menor_med=0;
		menor = 19e+10;
		for(j=0;j<n;++j){
			a = media_x - costu[j].x;
			b = media_y - costu[j].y;
			dista = sqrt(pow(a,2)+pow(b,2));
			if(dista < menor){
				menor_med = j;
				menor = dista;
			}
		}
		meds[i] = menor_med;
		y[meds[i]]=1;
		//printf("%d %f\t",menor_med,menor);
	}
	*/

}

/* Funcao para designar os vertices as medianas */
void desginar(costumer *costu, designa *assig, int medians, int n, double **dist, int *caps,int cap,int **x, int *meds){

	int temp,i,j,med_menor;
	double menor;

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
			for(j=0;j<medians;++j){
				if(costu[assig[i].custom].demand <= caps[j] && dist[assig[i].custom][meds[j]]<menor){
					menor=dist[assig[i].custom][meds[j]];
					med_menor = j;
				}

			}
			caps[med_menor]-=costu[assig[i].custom].demand;
			x[assig[i].custom][meds[med_menor]] = 1;
		}		
	}

}

/* Funcao de ordenacao dos vértices */
void sort(designa *assig,int n,double **dist,costumer *costu,int *meds,int medians){

	int i,j,k,temp,med_menor[2];
	double menor[2];

	/* Método bubblesort de ordenacao */
	for(i=0;i<n;++i){
		for(j=i+1;j<n;++j){
					
			/* Ordena em ordem crescente de distancia (1)
			 * alterar para <= para o caso (2)
			 * */
			if(dist[assig[j].custom][meds[assig[j].median]] < dist[assig[i].custom][meds[assig[i].median]]){
				temp = assig[i].custom;
				assig[i].custom = assig[j].custom;
				assig[j].custom = temp;

				temp = assig[i].median;
				assig[i].median = assig[j].median;
				assig[j].median = temp;
			}//end 1
		/*	Caso 3/4 - critério de desempate
		 *      se as distancias sao iguais, ordena pela ordem 
		 *	decrescente/crescente de demanda
		 *	(descomente para ativar)
		 */
			
		 
			else if(dist[assig[j].custom][meds[assig[j].median]] == dist[assig[i].custom][meds[assig[i].median]]){
				
				if(costu[j].demand < costu[i].demand){
					temp = assig[i].custom;
					assig[i].custom = assig[j].custom;
					assig[j].custom = temp;

					temp = assig[i].median;
					assig[i].median = assig[j].median;
					assig[j].median = temp;

				}
	
		/* Caso 5:
		 * verifica qual a segunda mediana menor
		 * e ordena pela ordem decrescente
		 * (descomente para ativar, deixando o caso 3/4 comentado
		 */
		/*	

				menor[0] = 9e+10;
				menor[1] = 9e+10;
				med_menor[0] = 0;
				med_menor[1] = 0;
				for(k=0;k<medians;++k){
					if(dist[assig[i].custom][meds[k]]<menor[0] && assig[i].median!=meds[k]){
						menor[0]=dist[assig[i].custom][meds[k]];
						med_menor[0] = k;
					}
					if(dist[assig[j].custom][meds[k]]<menor[1] && assig[j].median!=meds[k]){
						menor[1]=dist[assig[j].custom][meds[k]];
						med_menor[1] = k;
					}

				}
				
				if(dist[assig[j].custom][meds[med_menor[1]]] > dist[assig[i].custom][meds[med_menor[0]]]){
					temp = assig[i].custom;
					assig[i].custom = assig[j].custom;
					assig[j].custom = temp;

					temp = assig[i].median;
					assig[i].median = assig[j].median;
					assig[j].median = temp;
					
				}//end 5
		*/		
			}//end 3/4/5
		}
	}

}

/* Designa os candidatos segundo o método (6), calculando a razao entre L1 e L2 */
void desginar_rate(costumer *costu, designa *assig, int medians, int n, double **dist, int *caps,int cap,int **x, int *meds){

	int temp,i,j,med_menor, *ord, *l1, *l2;
	double menor, *L1, *L2, *r;
	
	L1 = (double *)malloc(sizeof(double)*n);
	L2 = (double *)malloc(sizeof(double)*n);
	r = (double *)malloc(sizeof(double)*n);
	l1 = (int *)malloc(sizeof(int)*n);
	l2 = (int *)malloc(sizeof(int)*n);
	ord = (int *)malloc(sizeof(int)*n);
	
	/* seta as capacidades de cada mediana */
	for(i=0;i<medians;++i) caps[i]=cap;
	
	/* encontra a mediana de menor distancia para cada vertice
	 * e coloca a distancia na array L1 e a mediana na array l1
	 * */
	for(i=0;i<n;++i){
		menor=9e+10;
		for(j=0;j<medians;++j){
			if(dist[i][meds[j]] < menor){
				menor=dist[i][meds[j]];
				med_menor=j;
			}
		}
		L1[i]=menor;
		l1[i]=med_menor;
	}

	/* encontra a mediana de segunda menor distancia
	 * para cada vertice e coloca a dist. na array L2
	 * e a mediana na array l2, e calcula a razao
	 */
	for(i=0;i<n;++i){
		menor=9e+10;
		for(j=0;j<medians;++j){
			if(dist[i][meds[j]] < menor && j!=l1[i]){
				menor=dist[i][meds[j]];
				med_menor=j;
			}
		}
		L2[i]=menor;
		l2[i]=med_menor;
		r[i]=L1[i]/L2[i];
	}

	/* inicializa a ordem de designação como sequencial */
	for(i=0;i<n;++i) ord[i]=i;

	/* bubblesort */
	for(i=0;i<n;++i)
		for(j=i+1;j<n;++j)
			if(r[ord[i]] > r[ord[j]])
				ord[i]^=ord[j]^=ord[i]^=ord[j];				


	/* designa os candidatos na ordem descrita na array ord
	 * caso a mediana nao tenha mais capacidade, atribui-se 
	 * para o segundo menor, e caso essa tb nao possa
	 * atribui-se pro proximo que caiba
	 * */
	for(i=0;i<n;++i){
		if(costu[ord[i]].demand <= caps[l1[ord[i]]]){
			caps[l1[ord[i]]]-=costu[ord[i]].demand;
			x[ord[i]][meds[l1[ord[i]]]]=1;
		}else if(costu[ord[i]].demand <= caps[l2[ord[i]]]){
			caps[l2[ord[i]]]-=costu[ord[i]].demand;
			x[ord[i]][meds[l2[ord[i]]]]=1;
		}else{
			for(j=0;j<medians;++j){
				if(costu[ord[i]].demand <= caps[j]){
					caps[j]-=costu[ord[i]].demand;
					x[ord[i]][meds[j]]=1;
					break;
				}
			}
		}
	}

	free(L1);
	free(L2);
	free(l1);
	free(l2);
	free(r);
	free(ord);

}
