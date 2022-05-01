CXX=g++
CC=$(CXX)
CXXFLAGS=-std=c++20
CPPFLAGS=-O3 -Wall

all: main_matrix

main_matrix: main_matrix.o mmio.o

main_matrix.o: main_matrix.cpp Matrix.hpp Matrix_imp.hpp

mmio.o: mmio.c mmio.h

clean:
	$(RM) *.o main_matrix
