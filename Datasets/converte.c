#include <stdio.h>

int main(){

	FILE *arq;
	int i,j,points,vertex,caps,demand;
	double x,y;
	char arqv[100];

	printf(" 11\n");
	for(i=1;i<12;++i){
		sprintf(arqv,"%d.txt",i);
		arq = fopen(arqv,"r");
		printf(" %d %d\n",i,0);
		fscanf(arq,"%d %d",&points,&vertex);
		printf(" %d %d",points,vertex);
		fscanf(arq,"%lf %lf %d %d",&x,&y,&caps,&demand);
		printf(" %d\n",caps);
		printf(" %d %.0f %.0f %d\n",1,x,y,demand);
		for(j=2;j<=points;++j){
			fscanf(arq,"%lf %lf %d %d",&x,&y,&caps,&demand);
			printf(" %d %.0f %.0f %d\n",j,x,y,demand);
		}
		fclose(arq);
	}
		
		

	return 0;
}
