CXX=g++
CC=$(CXX)
CXXFLAGS=-Wall -std=c++20 -fPIC
CPPFLAGS=-I${HOME}/PACS/pacs-examples/Examples/include
LDFLAGS=-L${HOME}/PACS/pacs-examples/Examples/lib -L. -lpacs -Wl,-rpath={HOME}/PACS/pacs-examples/Examples/lib

all: main_matrix

main_matrix: main_matrix.o mmio.o

main_matrix.o: main_matrix.cpp Matrix.hpp Matrix_imp.hpp

mmio.o: mmio.c mmio.h

clean:
	$(RM) *.o main_matrix
