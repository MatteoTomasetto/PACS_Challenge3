# Sparse Matrix class #

This program implements a sparse matrix class called Matrix that, in order to store non-zero values only, can support different formats:
	
	1) Uncompressed format: std::map with indices (i, j) as keys and the corresponding elements A_{i,j} as values;
	
	2) CSR Compressed Sparse Row format;

	3) CSC Compressed Sparse Column format;

Instead `main_matrix.cpp` makes some experiments and computations to check the functionalities in class Matrix; morever, thank to the Chrono utility, it allows to see which format is faster making operations.

Main features of Matrix class:
	
	1) Through a template parameter, you can select the storage oder (ROWMAJOR, COLUMNMAJOR);

	2) The `operator*()` compute matrix-vector multiplication efficiently with respect to the storage type and the format used;

	3) It is possible to compute the norm of the matrix thanks to the function member `norm()`: to select the wanted norm you can use different possible template parameters such as ONE, INF and FROBENIUS;

	4) `reader()` function member read the matrix in `.mtx` format thanks to the Marker Matrix Reader (see http://math.nist.gov/MatrixMarket for details).

The name of the file `.mtx` storing the matrix can be taken in input from command line thanks to GetPot:
Example of execution: `./main_matrix -f filename`
By default: filename = "lnsp_131.mtx".

In this directory, `make` produces the executable which is just called `main_matrix`.

ATTENTION: in order to compile the code with `make`, you should change PACS_ROOT_DIR to link correctly the library used for Chrono utilities.

