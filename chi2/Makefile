d1:	stats.cxx LR.h LR.c
	g++ -g -c LR.c -o LR.o -I. -std=c++11
	g++ -g -c -o stats.o  stats.cxx -I. -std=c++11
	g++ -g -o stats LR.o stats.o -lgomp -lnlopt -lgsl -lgslcblas -L/usr/lib -std=c++11

notebook:	stats.cxx LR.h LR.c
	g++ -g -c LR.c -o LR.o -I. -I/home/notebook/Minuit2-5.28.00/inc/ -std=c++11
	g++ -g -c -o stats.o  stats.cxx -I. -I/home/notebook/Minuit2-5.28.00/inc/ -std=c++11
	g++ -g -o stats LR.o stats.o /home/notebook/Minuit2-5.28.00/src/.libs/libMinuit2.a -lnlopt -lgomp -lgsl -lgslcblas -L/usr/lib -std=c++11


