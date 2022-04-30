PACS_ROOT_DIR?="${HOME}/PACS/pacs-examples/Examples"

CXX=g++
CC=$(CXX)
CXXFLAGS=-Wall -std=c++20 -fPIC
CPPFLAGS=-O3 -I${PACS_ROOT_DIR}/include
LDLIBS=-L${PACS_ROOT_DIR}/lib -lpacs -Wl,-rpath,${PACS_ROOT_DIR}/lib
all: main_matrix

main_matrix: main_matrix.o mmio.o

main_matrix.o: main_matrix.cpp Matrix.hpp Matrix_imp.hpp

mmio.o: mmio.c mmio.h

clean:
	$(RM) *.o main_matrix
