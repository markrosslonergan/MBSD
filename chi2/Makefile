d1:	stats.cxx LR.h LR.c
	g++ -g -std=c++0x -c LR.c -o LR.o -I. -I/mt/home/ballett/etc/include
	g++ -g -std=c++0x -c -o stats.o  stats.cxx -I. -I/mt/home/ballett/etc/include
	g++ -g -std=c++0x -o stats LR.o stats.o -lgomp -lnlopt -lgsl -lgslcblas -L/usr/lib -L/mt/home/ballett/etc/lib

#notebook:	stats.cxx LR.h LR.c
#	g++ -g -c LR.c -o LR.o -I. -I/home/notebook/Minuit2-5.28.00/inc/ -std=c++11
#	g++ -g -c -o stats.o  stats.cxx -I. -I/home/notebook/Minuit2-5.28.00/inc/ -std=c++11
#	g++ -g -o stats LR.o stats.o /home/notebook/Minuit2-5.28.00/src/.libs/libMinuit2.a -lnlopt -lgomp -lgsl -lgslcblas -L/usr/lib -std=c++11


