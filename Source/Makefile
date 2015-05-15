all: pmedian.o cheur.o lsheur.o ts.o ag.o aco.o geo.o
	gcc -o pmedian.exe pmedian.o cheur.o lsheur.o ts.o ag.o aco.o geo.o

aco.o: aco.c
	gcc -Wall -c aco.c
	
geo.o: geo.c
	gcc -Wall -c geo.c
	
ag.o: ag.c
	gcc -Wall -c ag.c
	
ts.o: ts.c
	gcc -Wall -c ts.c
	
lsheur.o: lsheur.c
	gcc -Wall -c lsheur.c
	
cheur.o: cheur.c
	gcc -Wall -c cheur.c

pmedian.o: pmedian.c
	gcc -Wall -c pmedian.c
