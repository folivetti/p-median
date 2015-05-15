#include <stdio.h>
#include <malloc.h>
#include "pmedian.h"
#include <time.h>

double p_calculate(int vertex, double **dist, int **x, int *freq, double k,int fi);


/* calcula penalidade baseado na frequencia */
double p_calculate(int vertex, double **dist, int **x, int *freq, double k,int fi){

	double answer=0.0;
	int i,j;

	/* para cada vertice...*/
	for(i=0;i<vertex;++i){
		/* ...para os outros vertices.. */
		for(j=0;j<vertex;++j){
			/* ... se for cliente->mediana, soma a distancia*/
			if(x[i][j]){
				answer += (dist[i][j]);
				if(j==fi)
					answer += k*freq[j];
			}
		}
	}

	
	return answer;

}

void ts(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int T, int max_it){

	int it,i,j,l,m;
	double answer,best_ans,lucro,*rank,new_ans,melhor;
	int **best_x,*best_y,*best_meds,*meds;
	int **new_x,*new_y,*new_meds;
	int *sorteio,*L1,*L2,*med_troca, *tabu1,*tabu2;
	int swap,fr,melhor_i,melhor_fr,velho,sai,melhorou;
	double k=0.0;
	int *freq1, *freq2, m_freq=-1,achou;

	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			if(dist[i][j] > k) k = dist[i][j];
	freq1 = (int *)malloc(sizeof(int)*n);
	freq2 = (int *)malloc(sizeof(int)*n);
	for(i=0;i<n;++i){
		freq1[i]=0;
		freq2[i]=0;
	}

	answer = calculate(costu, n, dist, medians, cap, x, y);
	
	meds = (int *)malloc(sizeof(int)*medians);
	best_meds = (int *)malloc(sizeof(int)*medians);
	new_meds = (int *)malloc(sizeof(int)*medians);
	L1 = (int *)malloc(sizeof(int)*n);
	L2 = (int *)malloc(sizeof(int)*n);

	best_y = (int *)malloc(sizeof(int)*n);
	best_x = (int **)malloc(sizeof(int *)*n);
	for(i=0;i<n;++i) best_x[i] = (int *)malloc(sizeof(int)*n);

	new_y = (int *)malloc(sizeof(int)*n);
	new_x = (int **)malloc(sizeof(int *)*n);
	for(i=0;i<n;++i) new_x[i] = (int *)malloc(sizeof(int)*n);


	for(i=0,j=0;i<n && j < medians;++i){
		if(y[i]) meds[j++]=i;
	}

	for(i=0;i<medians;++i) best_meds[i] = meds[i];



	sorteio = (int *)malloc(sizeof(int)*(n-medians));
	med_troca = (int *)malloc(sizeof(int)*(n-medians));
	rank = (double *)malloc(sizeof(double)*(n-medians));
	tabu1 = (int *)malloc(sizeof(int )*n);
	tabu2 = (int *)malloc(sizeof(int )*n);
	
	/* calcula o valor tabu de forma aleatoria entre [2n/15;8n/15]*/
	T = 2*n/15 + rand()*(6*n/15)/RAND_MAX;
	for(i=0;i<n;++i){ tabu1[i]=-T-1;tabu2[i]=-T-1; }

	atua_l(L1, L2, n, meds, medians, dist);

	
	best_ans=answer;
	for(l=0;l<n;++l) for(m=0;m<n;++m){ best_x[l][m]=x[l][m];new_x[l][m]=x[l][m];}
	for(l=0;l<n;++l){ best_y[l] = y[l]; new_y[l] = y[l]; }
	for(l=0;l<medians;++l){ best_meds[l]=meds[l]; new_meds[l]=meds[l];}

	for(it=0;it<max_it;++it){

		/* se ficou 20% das iteracoes sem melhora, atualiza tempo tabu */
		if(melhorou > max_it*20/100){
			T = 2*n/15 + rand()*(6*n/15)/RAND_MAX;
			melhorou=0;
			printf("tabu: %d\n",T);
		}

		/* se for os 85% finais fixa a mediana mais usada na solucao caso ela ja nao esteja la */
		if(it>85*max_it/100 && m_freq==-1){
			m_freq=0;
			for(i=1;i<n;++i){
				if(freq1[i] > freq1[m_freq]){
					m_freq=i;
				}
			}
			achou=0;
			for(i=0;i<medians;++i){
				if(meds[i]==m_freq)
					achou=1;
			}
			if(!achou){
				melhor_subs(meds,m_freq,L1,L2,n,medians,dist,&fr,&lucro);
				y[meds[fr]]=0;
				y[m_freq]=1;
				meds[fr]=m_freq;
				for(l=0;l<n;++l) for(m=0;m<n;++m) new_x[l][m]=0;
				redesignar(costu, medians,n,dist,cap,new_x,new_meds,m_freq,fr,L1,L2);
				answer = calculate(costu, n, dist, medians, cap, x, y);
				reatua_l(L1, L2, n, meds, medians, dist,m_freq,fr);
				if(answer < best_ans && factivel(costu,n,medians,cap,x,y)){
					printf("%d melhorou: %f\n",it,answer);
					best_ans=answer;
					for(l=0;l<n;++l) for(m=0;m<n;++m) best_x[l][m]=x[l][m];
					for(l=0;l<n;++l) best_y[l] = y[l];
					for(l=0;l<medians;++l) best_meds[l]=meds[l];
				}
			}
		}
		melhor=9e+10;
		for(l=0;l<n;++l) for(m=0;m<n;++m){ ;new_x[l][m]=x[l][m];}
		for(l=0;l<n;++l){  new_y[l] = y[l]; }
		for(l=0;l<medians;++l){ new_meds[l]=meds[l];}
		/* para cada cliente fora da solucao  */
		for(i=0;i<n;++i){
			if(!y[i]){
				/* calcula o lucro */
				melhor_subs(meds,i,L1,L2,n,medians,dist,&fr,&lucro);
				if(meds[fr]!=m_freq){
					/*se for positivo faz todo o processo da busca local */
				if(lucro>0){
					new_y[new_meds[fr]]=0;
					new_y[i]=1;
					velho = new_meds[fr];
					new_meds[fr]=i;
					for(l=0;l<n;++l) for(m=0;m<n;++m) new_x[l][m]=0;
					redesignar(costu, medians,n,dist,cap,new_x,new_meds,i,fr,L1,L2);

					/* caso esteja no periodo de diversificacao aplica penalizacao
					 * na resposta
					 * */
					if(it>45*max_it/100 && it < 75*max_it/100){
						new_ans = p_calculate(n, dist, new_x, freq1, k,i);
					}else{
						new_ans = calculate(costu, n, dist, medians, cap, new_x, new_y);
					}
					/* entra no processo de escolha se
					 * a resposta nao for tabu ou 
					 * se a resposta for melhor que a melhor de todas, aspira 
					 * */
					if(new_ans < best_ans || it-T>tabu2[i]){
					if(new_ans < melhor){
						melhor_i = i;
						melhor_fr = fr;
						melhor = new_ans;
					}
					}
					new_y[velho]=1;
					new_y[i]=0;
					new_meds[fr]=velho;
				}
				}
			}
		}
				
		/* executa o melhor nao tabu e atualiza lista de frequencia */
		++freq1[melhor_i];
		++freq2[meds[melhor_fr]];
		y[melhor_i]=1;
		y[meds[melhor_fr]]=0;
		sai = meds[melhor_fr];
		meds[melhor_fr]=melhor_i;
		for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
		redesignar(costu, medians,n,dist,cap,x,meds,melhor_i,melhor_fr,L1,L2);
		if(factivel(costu,n,medians,cap,x,y)) lsheur2(costu, n, dist, medians, cap, x, y, 10);
		for(i=0,j=0;i<n && j < medians;++i){
			if(y[i]) meds[j++]=i;
		}


		answer = calculate(costu, n, dist, medians, cap, x, y);
		reatua_l(L1, L2, n, meds, medians, dist,melhor_i,melhor_fr);

		/* Verifica se e´ melhor do que a melhor encontrada ate agora */
		if(answer < best_ans && factivel(costu,n,medians,cap,x,y)){
			printf("%d melhorou: %f\n",it,answer);
			best_ans=answer;
			for(l=0;l<n;++l) for(m=0;m<n;++m) best_x[l][m]=x[l][m];
			for(l=0;l<n;++l) best_y[l] = y[l];
			for(l=0;l<medians;++l) best_meds[l]=meds[l];
			melhorou=0;
		}else{
			++melhorou;
		}
		
		/* seta movimento como tabu */
		tabu1[sai]=it;
		tabu2[melhor_i]=it;
	}
	for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=best_x[l][m];
	for(l=0;l<n;++l) y[l] = best_y[l];
	for(l=0;l<medians;++l) meds[l]=best_meds[l];
	answer = best_ans;
	answer = refine(x, y, n, medians, answer, dist,cap,meds,costu);

	free(meds);
	free(best_meds);
	free(new_meds);
	free(L1);
	free(L2);
	free(best_y);
	for(i=0;i<n;++i) free(best_x[i]);
	free(best_x);
	free(new_y);
	for(i=0;i<n;++i) free(new_x[i]);
	free(new_x);
	free(sorteio);
	free(med_troca);
	free(rank);
	free(tabu1);
	free(tabu2);

}
