#include <stdio.h>
#include <malloc.h>
#include "pmedian.h"
#include <time.h>

//#define min(x,y) ((x<y)?x:y)



double refine(int **x, int *y, int n, int medians, double answer, double **dist,int cap, int *meds, costumer *costu){

	int *caps, *cli_med;
	int i,j;
	double new_ans;
	
	caps = (int *)malloc(sizeof(int)*medians);
	cli_med = (int *)malloc(sizeof(int)*n);

	for(i=0;i<medians;++i){
		caps[i]=cap;
	}

	for(i=0;i<medians;++i){
		for(j=0;j<n;++j){
			if(x[j][meds[i]]){
				caps[i]-=costu[j].demand;
				cli_med[j]=i;
			}
		}
	}

	/* mesma coisa só que pros clientes */
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			new_ans = answer - dist[i][meds[cli_med[i]]] - dist[j][meds[cli_med[j]]] + \
				dist[i][meds[cli_med[j]]] + dist[j][meds[cli_med[i]]];
			
			if(new_ans < answer && caps[cli_med[i]]+costu[i].demand-costu[j].demand >= 0 \
					&& caps[cli_med[j]]+costu[j].demand-costu[i].demand >= 0){
			
				x[i][meds[cli_med[j]]] = 1;
				x[i][meds[cli_med[i]]] = 0;
				x[j][meds[cli_med[i]]] = 1;
				x[j][meds[cli_med[j]]] = 0;
				caps[cli_med[i]] += costu[i].demand-costu[j].demand;
				caps[cli_med[j]] += costu[j].demand-costu[i].demand;
				cli_med[j]^=cli_med[i]^=cli_med[j]^=cli_med[i];
				answer=new_ans;
			}
		}
	}

	free(caps);
	free(cli_med);

	return answer;
}

int factivel(costumer *costum,int n, int medians, int cap, int **x, int *y){

	int soma,i,j,med=-1;
	
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
		if(soma > cap) return 0;
	}

	soma=0;
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			soma+=x[i][j];
			if(x[i][j] && !y[j]) return 0;
		}
	}
	if(soma != n) return 0;

	soma=0;
	for(i=0;i<n;++i) soma+=y[i];
	if(soma != medians) return 0;
	
	return 1;
}


/* Resende */


void lsheur2(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it){
	
	int it,i,j,l,m,mudou=1,fr,trocas=0;
	int *meds, *L1, *L2;
	double answer, new_ans, lucro;
	int **new_x, *new_y, *new_meds;

	
	answer = calculate(costu, n, dist, medians, cap, x, y);
	
	meds = (int *)malloc(sizeof(int)*medians);
	new_meds = (int *)malloc(sizeof(int)*medians);
	L1 = (int *)malloc(sizeof(int)*n);
	L2 = (int *)malloc(sizeof(int)*n);

	new_y = (int *)malloc(sizeof(int)*n);
	new_x = (int **)malloc(sizeof(int *)*n);
	for(i=0;i<n;++i) new_x[i] = (int *)malloc(sizeof(int)*n);


	for(i=0,j=0;i<n && j < medians;++i){
		if(y[i]) meds[j++]=i;
	}

	for(i=0;i<medians;++i) new_meds[i] = meds[i];
	for(i=0;i<n;++i){
		new_y[i]=y[i];
		for(j=0;j<n;++j){
			new_x[i][j]=x[i][j];
		}
	}

	atua_l(L1, L2, n, meds, medians, dist);
	
	for(it=0;it<max_it && mudou;++it){

		mudou=0;
		/* pra cada vertice */

		for(j=0;j<n;++j){

			if(!y[j]){

				melhor_subs(meds,j,L1,L2,n,medians,dist,&fr,&lucro);
				if(lucro > 0){
					
					new_y[j]=1;
					new_y[new_meds[fr]]=0;
					new_meds[fr]=j;
					for(l=0;l<n;++l) for(m=0;m<n;++m) new_x[l][m]=0;
				
					redesignar(costu, medians,n,dist,cap,new_x,new_meds,j,fr,L1,L2);
					new_ans = calculate(costu, n, dist, medians, cap, new_x, new_y);

					if(new_ans < answer && factivel(costu,n,medians,cap,new_x,new_y)){
						
						y[j]=1;
						y[meds[fr]]=0;
						meds[fr]=j;
						for(l=0;l<n;++l)
							for(m=0;m<n;++m)
								x[l][m]=new_x[l][m];
						answer = new_ans;
						mudou=1;
						++trocas;
						reatua_l(L1, L2, n, meds, medians, dist,j,fr);
					}else{
						new_y[j]=0;
						new_y[meds[fr]]=1;
						new_meds[fr]=meds[fr];
						
					}
					
				}

			}
			
		}
		
	}
	
	if(factivel(costu,n,medians,cap,x,y)){
		answer = refine(x, y, n, medians, answer, dist,cap,meds,costu);
	}
				
	free(meds);
	free(L1);
	free(L2);
	free(new_meds);
	free(new_y);
	for(i=0;i<n;++i) free(new_x[i]);
	free(new_x);
	
}	


double d_min(double a, double b){

	if(a<b) return a;
	else return b;
}

void melhor_subs(int *meds,int fi, int *L1, int *L2, int n, int medians, double **d, int *fr,double *lucro){

	double w=0.0, *v;
	int i;

	v = (double *)malloc(sizeof(double)*medians);

	for(i=0;i<medians;++i) 
		v[i]=0.0;
	for(i=0;i<n;++i){
		if(d[i][fi] < d[i][meds[L1[i]]]){
			w += (d[i][meds[L1[i]]] - d[i][fi]);
		}else{
			v[L1[i]] += d_min(d[i][fi],d[i][meds[L2[i]]]) - d[i][meds[L1[i]]];
		}
	}

	*fr = 0;
	for(i=1;i<medians;++i){
		if(v[i] < v[*fr]){
			*fr=i;
		}
	}
	//printf("w: %f\n",w);
	//printf("v: %f\n\n",v[*fr]);
	*lucro = w - v[*fr];

	free(v);
}

void redesignar(costumer *costu, int medians, int n, double **dist, int cap,int **x, int *meds, int fi, int fr, int *L1, int *L2){

	int i,j,med_menor;
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
		if(L1[i]==fr){
			if(dist[i][fi] < dist[i][meds[L2[i]]]){
				assig[i].median=fr;
			}else{
				assig[i].median=L2[i];
			}
		}else{
			if(dist[i][fi] < dist[i][meds[L1[i]]]){
				assig[i].median=fr;
			}else{
				assig[i].median=L1[i];
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


void atua_l(int *L1, int *L2, int n, int *meds, int medians, double **dist){

	int i,j;

	for(i=0;i<n;++i){
		L1[i] = 0;
		L2[i] = 0;
		for(j=1;j<medians;++j){
			if(dist[i][meds[j]] < dist[i][meds[L1[i]]]){
				L1[i] = j;
			}
		}
		
		if(L2[i] == L1[i]) L2[i]=(L2[i]+1)%n;
		
		
		for(j=0;j<medians;++j){
			if((dist[i][meds[j]] < dist[i][meds[L2[i]]]) && L1[i] != j){
				L2[i] = j;
			}
		}
		
	}
}

void reatua_l(int *L1, int *L2, int n, int *meds, int medians, double **dist,int fi, int fr){

	int i,j,menor_j;
	double menor;

	for(i=0;i<n;++i){
		if(L1[i]==fr){
			if(dist[i][fi] < dist[i][meds[L2[i]]]){
				L1[i]=fr;
			}else{
				L1[i]=L2[i];
				L2[i]=fr;
			}
		}else{
			if(dist[i][fi] < dist[i][meds[L1[i]]]){
				L1[i]=fr;
			}
		}
		
		if(L2[i]==fr){
			if(L1[i] != fr){
				menor=dist[i][meds[fr]];
				menor_j=fr;
				for(j=0;j<medians;++j){
					if(dist[i][meds[j]] < menor && j!=fr && j!=L1[i]){
						menor = dist[i][meds[j]];
						menor_j = j;
					}
				}
				L2[i]=menor_j;
			}else{
				menor=9e+10;
				menor_j=0;
				for(j=0;j<medians;++j){
					if(dist[i][meds[j]] < menor && j!=fr && j!=L1[i]){
						menor = dist[i][meds[j]];
						menor_j = j;
					}
				}
				L2[i]=menor_j;
			}
			
		}else{
			if(L1[i]!=fr && dist[i][fi] < dist[i][meds[L2[i]]]){
				L2[i]=fr;
			}
		}
				
	}
}

