CXX=g++
CC=$(CXX)
CXXFLAGS=-Wall -std=c++20
all: main_matrix

main_matrix: main_matrix.o

main_matrix.o: main_matrix.cpp Matrix.hpp Matrix_imp.hpp

clean:
	$(RM) *.o main_matrix

