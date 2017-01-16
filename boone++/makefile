all:
	g++ boone.cxx -o boone.out globals.cxx  `pkg-config --cflags --libs gsl` -std=c++11  -I`globes-config --include` -L`globes-config --libs` -lglobes -lboost_program_options
#-static for movable code
