#include <stdio.h>
#include <malloc.h>
#include "pmedian.h"
#include <time.h>

void ADD(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq,int it);
void DELETE(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq);
void SWAP(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq,int it);
int checa_best(double answer, int **x, int *y, int *meds, double *best_answer, int **best_x, int *best_y, int *best_meds, int n, int medians,int it);
double p_calculate(int vertex, double **dist, int **x, int *freq, double k,int fi);

void ts4(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int T, int max_it){

	/* ADD, DELETE, SWAP */
	int **best_x, *best_y, *best_meds;
	double best_answer,answer;
	int *L1, *L2, *meds, it, *tabu, *freq,i,j;
	int melhorou, qt_mel=0;

	meds = (int *)malloc(sizeof(int)*(medians+1));
	best_meds = (int *)malloc(sizeof(int)*medians);
	L1 = (int *)malloc(sizeof(int)*n);
	L2 = (int *)malloc(sizeof(int)*n);
	tabu = (int *)malloc(sizeof(int)*n);
	freq = (int *)malloc(sizeof(int)*n);

	best_y = (int *)malloc(sizeof(int)*n);
	best_x = (int **)malloc(sizeof(int *)*n);
	for(i=0;i<n;++i) best_x[i] = (int *)malloc(sizeof(int)*n);

	for(i=0,j=0;i<n && j < medians;++i){
		if(y[i]) meds[j++]=i;
	}

	for(i=0;i<medians;++i)
		best_meds[i]=meds[i];

	for(i=0;i<n;++i){
		best_y[i]=y[i];
		freq[i]=0;
		for(j=0;j<n;++j){
			best_x[i][j]=x[i][j];
		}
	}


	srand((int)time(NULL));
	/* [n/6,n/4] */
	T = n/15 + rand()*(2*n/15)/RAND_MAX;
	for(i=0;i<n;++i) tabu[i] = -T-1;

	best_answer = calculate(costu, n, dist, medians, cap, x, y);

	atua_l(L1, L2, n, meds, medians, dist);

	for(it=0;it<max_it;++it){

		
		DELETE(costu, n, dist, medians, cap, x, y, meds, T,L1,L2,tabu,freq);
		answer = calculate(costu, n, dist, medians-1, cap, x, y);
		//update_tabu(tabu,n);

		ADD(costu, n, dist, medians-1, cap, x, y, meds, T,L1,L2,tabu,freq,it);
		answer = calculate(costu, n, dist, medians, cap, x, y);
		//update_tabu(tabu,n);

		if(factivel(costu,n,medians,cap,x,y)){ 
			melhorou = checa_best(answer,x,y,meds,&best_answer,best_x,best_y,best_meds,n,medians,it);
			if(!melhorou)
				++qt_mel;
			else
				qt_mel=0;
		}

		SWAP(costu, n, dist, medians, cap, x, y,meds, T,L1,L2,tabu,freq,it);
		answer = calculate(costu, n, dist, medians, cap, x, y);
		//update_tabu(tabu,n);
		
		if(factivel(costu,n,medians,cap,x,y)){
			melhorou = checa_best(answer,x,y,meds,&best_answer,best_x,best_y,best_meds,n,medians,it);
			if(!melhorou)
				++qt_mel;
			else
				qt_mel=0;
		}

		ADD(costu, n, dist, medians, cap, x, y, meds, T,L1,L2,tabu,freq,it);
		answer = calculate(costu, n, dist, medians+1, cap, x, y);
		//update_tabu(tabu,n);
		
		DELETE(costu, n, dist, medians+1, cap, x, y,meds, T,L1,L2,tabu,freq);
		answer = calculate(costu, n, dist, medians, cap, x, y);
		//update_tabu(tabu,n);
	
		if(factivel(costu,n,medians,cap,x,y)){
			melhorou = checa_best(answer,x,y,meds,&best_answer,best_x,best_y,best_meds,n,medians,it);
			if(!melhorou)
				++qt_mel;
			else
				qt_mel=0;
		}
	
		SWAP(costu, n, dist, medians, cap, x, y, meds,T,L1,L2,tabu,freq,it);
		answer = calculate(costu, n, dist, medians, cap, x, y);
		//update_tabu(tabu,n);
		
		if(factivel(costu,n,medians,cap,x,y)){
			melhorou = checa_best(answer,x,y,meds,&best_answer,best_x,best_y,best_meds,n,medians,it);
			if(!melhorou)
				++qt_mel;
			else
				qt_mel=0;
		}
		if(qt_mel > max_it*20/100){
			//T = n/6 + rand()*(n/12)/RAND_MAX;
			T = n/15 + rand()*(2*n/15)/RAND_MAX;
			qt_mel=0;
		}
	}

	for(i=0;i<n;++i){
		y[i]=best_y[i];
		for(j=0;j<n;++j)
			x[i][j]=best_x[i][j];
	}
	for(i=0;i<medians;++i) meds[i] = best_meds[i];

	if(factivel(costu,n,medians,cap,x,y))
		lsheur2(costu, n, dist, medians, cap, x, y, 10);

	free(meds);
	free(best_meds);
	free(best_y);
	for(i=0;i<n;++i) free(best_x[i]);
	free(best_x);
	free(L1);
	free(L2);
	free(tabu);
	free(freq);
	
}

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

int checa_best(double answer, int **x, int *y, int *meds, double *best_answer, int **best_x, int *best_y, int *best_meds, int n, int medians,int it){

	int i,j;

	if(answer < *best_answer){
		printf("%d %f %f\n",it, answer, *best_answer);
		for(i=0;i<n;++i){
			best_y[i] = y[i];
			for(j=0;j<n;++j)
				best_x[i][j]=x[i][j];
		}
		for(i=0;i<medians;++i)
			best_meds[i]=meds[i];
		*best_answer = answer;
		return 1;
	}
	return 0;
}

void ADD(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq,int it){
	
	int i,fi,melhor_fi=0,l,m;
	double w,melhor_lucro=9e+10;
	static int k;

	if(k==0){
		for(l=0;l<n;++l)
			for(m=0;m<n;++m)
				if(dist[l][m] > k)
					k=dist[l][m];
	}

	atua_l(L1, L2, n, meds, medians, dist);
	for(fi=0;fi<n;++fi){
		w=0.0;
		if(!y[fi] && it-T>tabu[fi]){
		/*	
			for(i=0;i<n;++i){
				if(dist[i][fi] < dist[i][meds[L1[i]]]){
					w += (dist[i][meds[L1[i]]] - dist[i][fi]);
					w -= k*freq[fi];
				}
			}
		*/	
			meds[medians] = fi;
			for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
			designar(costu, medians+1,n,dist,cap,x,meds);
			y[fi]=1;
			//w = calculate(costu, n, dist, medians+1, cap, x, y);
			//w += k*freq[fi]*n/medians;
			w = p_calculate(n, dist, x, freq, k,fi);
			y[fi]=0;
			
			if(w < melhor_lucro){
				melhor_lucro = w;
				melhor_fi = fi;
			}
		}
	}

	fi = melhor_fi;
	y[fi] = 1;
	
	//meds = (int *)realloc(meds,sizeof(int)*(medians+1));
	meds[medians]=fi;
	for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
	designar(costu, medians+1,n,dist,cap,x,meds);
	tabu[fi]=it;
 	++freq[fi];

//	for(i=0;i<medians+1;++i) printf("%d\t",meds[i]);

	return;

}

void DELETE(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq){

	int i,fr,melhor_fr=0,l,m,old;
	double v,melhor_lucro=1e+9;

	//for(i=0;i<medians;++i) printf("%d\t",meds[i]);
	atua_l(L1, L2, n, meds, medians, dist);
	for(fr=0;fr<medians;++fr){
		v=0.0;
		//if(!tabu[meds[fr]]){
		/*
			for(i=0;i<n;++i){
				if(x[i][meds[fr]]){
					v += dist[i][meds[L2[i]]] - dist[i][meds[fr]];
				}				
			}
		*/
			y[meds[fr]]=0;
			if(fr!=medians-1){
				old = meds[fr];
				meds[fr]=meds[medians-1];
			}
			for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
			designar(costu, medians-1,n,dist,cap,x,meds);
			y[old]=1;
			if(fr!=medians-1){
				meds[medians-1]=meds[fr];
				meds[fr]=old;
			}
			v = calculate(costu, n, dist, medians-1, cap, x, y);
			
			if(v < melhor_lucro){
				melhor_lucro = v;
				melhor_fr = fr;
			}
		//}
	}

	fr = melhor_fr;
	y[meds[fr]]=0;
	if(fr!=medians-1){
		meds[fr] = meds[medians-1];
	}
//	meds = (int *)realloc(meds,sizeof(int)*(medians-1));
	for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
	designar(costu, medians-1,n,dist,cap,x,meds);

}

void SWAP(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int *meds, int T,int *L1,int *L2,int *tabu,int *freq,int it){

	int fi,fr,l,m,i,velho;
	double lucro,melhor_lucro=9e+10,new_ans;
	int melhor_fi=0,melhor_fr=0;
	static int k;
	int *new_y, *new_meds;

	new_meds = (int *)malloc(sizeof(int)*medians);
	new_y = (int *)malloc(sizeof(int)*n);

	for(i=0;i<medians;++i) new_meds[i] = meds[i];

	for(i=0;i<n;++i) new_y[i]=y[i];
	
	if(k==0){
		for(l=0;l<n;++l)
			for(m=0;m<n;++m)
				if(dist[l][m] > k)
					k=dist[l][m];
	}

	/* = busca local */
	atua_l(L1, L2, n, meds, medians, dist);
	for(fi=0;fi<n;++fi){
		if(!y[fi] && it-T>tabu[fi]){
			melhor_subs(meds,fi,L1,L2,n,medians,dist,&fr,&lucro);
			//lucro -= k*freq[fi]*n/medians;
			//if(lucro>0){
				new_y[new_meds[fr]]=0;
				new_y[fi]=1;
				velho = new_meds[fr];
				new_meds[fr]=fi;
				for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
				redesignar(costu, medians,n,dist,cap,x,new_meds,fi,fr,L1,L2);
				//new_ans = calculate(costu, n, dist, medians, cap,x, new_y);
				//new_ans+=n*k*(freq[i])/medians;
				new_ans = p_calculate(n, dist, x, freq, k,fi);
				if(new_ans < melhor_lucro){
					melhor_fi = fi;
					melhor_fr = fr;
					melhor_lucro = new_ans;
				}
				new_y[velho]=1;
				new_y[fi]=0;
				new_meds[fr]=velho;

			//}
			/*
			if(lucro > melhor_lucro){
				melhor_lucro = lucro;
				melhor_fi = fi;
				melhor_fr = fr;
			}
			*/
		}
	}
	fr = melhor_fr;
	fi = melhor_fi;
	y[meds[fr]]=0;
	y[fi]=1;
	meds[fr]=fi;
	for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
	redesignar(costu, medians,n,dist,cap,x,meds,fi,fr,L1,L2);
	tabu[fi] = it;
	++freq[fi];
	
	free(new_y);
	free(new_meds);

}
/*
void path_relink(){

	for(i=0;i<medians;++i){
		for(j=0;j<medians;++j){
			if(meds1[i]==meds2[j]){
				igual=1;
			}
		}
		if(igual==0){
			swap1[s]=
		}
}
*/
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
	//for(i=0;i<n;++i) tabu[i] = (int *)malloc(sizeof(int)*n);
	
	T = 2*n/15 + rand()*(6*n/15)/RAND_MAX;
	for(i=0;i<n;++i){ tabu1[i]=-T-1;tabu2[i]=-T-1; }

	atua_l(L1, L2, n, meds, medians, dist);

	
	best_ans=answer;
	for(l=0;l<n;++l) for(m=0;m<n;++m){ best_x[l][m]=x[l][m];new_x[l][m]=x[l][m];}
	for(l=0;l<n;++l){ best_y[l] = y[l]; new_y[l] = y[l]; }
	for(l=0;l<medians;++l){ best_meds[l]=meds[l]; new_meds[l]=meds[l];}

	for(it=0;it<max_it;++it){
		if(melhorou > max_it*20/100){
			T = 2*n/15 + rand()*(6*n/15)/RAND_MAX;
			melhorou=0;
			printf("tabu: %d\n",T);
		}
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
		
		for(i=0;i<n;++i){
			if(!y[i]){// && it-T>tabu2[i]){
				melhor_subs(meds,i,L1,L2,n,medians,dist,&fr,&lucro);
				if(meds[fr]!=m_freq){
				//if(lucro>0){// && it-T>tabu1[fr]){
					new_y[new_meds[fr]]=0;
					new_y[i]=1;
					velho = new_meds[fr];
					new_meds[fr]=i;
					for(l=0;l<n;++l) for(m=0;m<n;++m) new_x[l][m]=0;
					redesignar(costu, medians,n,dist,cap,new_x,new_meds,i,fr,L1,L2);
					if(it>45*max_it/100 && it < 75*max_it/100){
						new_ans = p_calculate(n, dist, new_x, freq1, k,i);
					}else{
						new_ans = calculate(costu, n, dist, medians, cap, new_x, new_y);
					}
					//new_ans+=n*k*(freq1[i])/medians;
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
				
		/* executa o melhor nao tabu */
		++freq1[melhor_i];
		++freq2[meds[melhor_fr]];
		y[melhor_i]=1;
		y[meds[melhor_fr]]=0;
		sai = meds[melhor_fr];
		meds[melhor_fr]=melhor_i;
		for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
		//designar(costu, medians,n,dist,cap,x,meds);
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
		
		/* seta movimento como tabu e atualiza lista tabu */
		//update_tabu(tabu1,n);
		//update_tabu(tabu2,n);
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
	//for(i=0;i<n;++i) free(tabu[i]);
	free(tabu1);
	free(tabu2);

}

void ts2(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int T, int max_it){

	int it,i,j,l,m;
	double answer,best_ans,lucro,*rank;
	int **best_x,*best_y,*best_meds,*meds;
	int **new_x,*new_y,*new_meds;
	int *sorteio,*L1,*L2,*med_troca, *tabu;
	int swap,fr,melhor,melhor_i,sai;

	srand((int)time(NULL));
	
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
	tabu = (int *)malloc(sizeof(int )*n);
	//for(i=0;i<n;++i) tabu[i] = (int *)malloc(sizeof(int)*n);
	
	for(i=0;i<n;++i) tabu[i]=0;

	//answer = refine(x, y, n, medians, answer, dist,cap,meds,costu);
	atua_l(L1, L2, n, meds, medians, dist);

	
	best_ans=answer;
	for(l=0;l<n;++l) for(m=0;m<n;++m){ best_x[l][m]=x[l][m];new_x[l][m]=x[l][m];}
	for(l=0;l<n;++l){ best_y[l] = y[l]; new_y[l] = y[l]; }
	for(l=0;l<medians;++l){ best_meds[l]=meds[l]; new_meds[l]=meds[l];}

	for(it=0;it<max_it;++it){
		for(i=0,j=0;i<n && j<n-medians;++i){
			if(!y[i]){
				sorteio[j]=i;
				++j;
			}
		}
		for(i=0;i<n-medians;++i){
			swap = rand()%(n-medians);
			sorteio[i]^=sorteio[swap]^=sorteio[i]^=sorteio[swap];
		}
		/* seleciona os movimentos aleatorios e rankeia */
		for(i=0;i<n-medians;++i){
			swap = rand()%(n-medians);
			sorteio[i]^=sorteio[swap]^=sorteio[i]^=sorteio[swap];
		}
		for(i=0;i<n-medians;++i){
			j = sorteio[i];
			melhor_subs(meds,j,L1,L2,n,medians,dist,&fr,&lucro);
			rank[i]=lucro;
			med_troca[i]=fr;
		}
				
		/* executa o melhor nao tabu */
		melhor=0;melhor_i=-1;
		for(i=0;i<n-medians;++i)
			if(!tabu[sorteio[i]]) {melhor=sorteio[i];melhor_i=i;break;}
		for(;i<n-medians;++i){
			if(rank[i] > rank[melhor_i] && !tabu[sorteio[i]]){
				melhor=sorteio[i];
				melhor_i=i;
			}
		}
		if(melhor_i==-1){ printf("noway\n");continue;}
		y[melhor]=1;
		sai = meds[med_troca[melhor_i]];
		y[meds[med_troca[melhor_i]]]=0;
		meds[med_troca[melhor_i]]=j;
		for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=0;
		redesignar(costu, medians,n,dist,cap,x,meds,melhor,med_troca[melhor_i],L1,L2);
		lsheur2(costu, n, dist, medians, cap, x, y, 10);
		answer = calculate(costu, n, dist, medians, cap, x, y);
		reatua_l(L1, L2, n, meds, medians, dist,melhor,med_troca[melhor_i]);
//printf("%f\n",answer);
		/* Verifica se e´ melhor do que a melhor encontrada ate agora */
		if(answer < best_ans && factivel(costu,n,medians,cap,x,y)){
			printf("melhorou: %f\n",answer);
			best_ans=answer;
			for(l=0;l<n;++l) for(m=0;m<n;++m) best_x[l][m]=x[l][m];
			for(l=0;l<n;++l) best_y[l] = y[l];
			for(l=0;l<medians;++l) best_meds[l]=meds[l];
		}
		
		/* seta movimento como tabu e atualiza lista tabu */
		update_tabu(tabu,n);
		tabu[sai]=T;
	}
	for(l=0;l<n;++l) for(m=0;m<n;++m) x[l][m]=best_x[l][m];
	for(l=0;l<n;++l) y[l] = best_y[l];
	for(l=0;l<medians;++l) meds[l]=best_meds[l];

	free(meds);
	free(best_meds);
	free(L1);
	free(L2);
	free(best_y);
	for(i=0;i<n;++i) free(best_x[i]);
	free(best_x);
	free(sorteio);
	free(med_troca);
	free(rank);
	//for(i=0;i<n;++i) free(tabu[i]);
	free(tabu);

}

void update_tabu(int *tabu,int n){

	int i;
	
	for(i=0;i<n;++i)
		if(tabu[i])
			--tabu[i];
}
