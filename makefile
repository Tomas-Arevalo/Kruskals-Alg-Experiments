all: randmst

randmst: randmst.cpp
	g++ -std=c++11 -O3 -o randmst randmst.cpp
