d1:	stats.cxx LR.h LR.c
	g++ -g -c LR.c -o LR.o -I. 
	g++ -g -c -o stats.o  stats.cxx -I. 
	g++ -g -o stats LR.o stats.o -lgomp -lgsl -lgslcblas -L/usr/lib

notebook:	stats.cxx LR.h LR.c
	g++ -g -c LR.c -o LR.o -I. -I/home/notebook/Minuit2-5.28.00/inc/
	g++ -g -c -o stats.o  stats.cxx -I. -I/home/notebook/Minuit2-5.28.00/inc/
	g++ -g -o stats LR.o stats.o /home/notebook/Minuit2-5.28.00/src/.libs/libMinuit2.a -lgomp -lgsl -lgslcblas -L/usr/lib


