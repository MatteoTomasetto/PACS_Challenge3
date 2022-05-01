#include <iostream>
#include <complex>
#include <chrono>
#include "Matrix.hpp"
#include "GetPot"

using timer = std::chrono::high_resolution_clock;

int main(int argc, char **argv)
{

using namespace apsc;
{
	std::cout << "Creating matrix A uncompressed - ROWMAJOR... " << std::endl;

	auto n = 4u;
	auto m = 5u;
	Matrix<double,ROWMAJOR> A(n,m);
	A.ref(0,0) = 1;
	A.ref(3,2) = -1;

	std::cout << "nrows= " << A.nRows << " " << "ncols= " << A.nCols << std::endl;
	for (auto i = 0u; i < n; ++i)
	{
		for (auto j = 0u; j<m; ++j)
		{
			std::cout << A(i,j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v{1,2,3,4,5};
	std::vector<double> res = A*v;
	for (auto i = 0u; i < A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;


	// Matrix-vector multiplication with a matrix (nCols x 1)
	std::cout << "Computing A*v with v a matrix (nCols x 1)..." << std::endl;
	Matrix<double, ROWMAJOR> Vmat(1,5);
	Vmat.ref(0,0) = 1;
	Vmat.ref(0,1) = 2;
	Vmat.ref(0,2) = 3;
	Vmat.ref(0,3) = 4;
	Vmat.ref(0,4) = 5;
	res = A*Vmat;
	for (auto i = 0u; i < A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;


	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;

	// Compressed case
	std::cout << "Compressed case..." << std::endl;
	A.makeCompressed();
	for (auto i = 0u; i < n; ++i)
	{
		for (auto j = 0u; j < m; ++j)
		{
			std::cout << A(i,j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	auto inner = A.getInnerIndexes();
	auto outer = A.getOuterIndexes();
	auto values = A.getValues();
	std::cout << "#non zero elements=" << A.nnz() << std::endl;

	std::cout << " Inner indexes:";
	for (auto i:inner) std::cout << i << " ";
	std::cout << std::endl;
	std::cout << " Outer indexes:";
	for (auto i:outer) std::cout << i <<" ";
	std::cout << std::endl;
	std::cout << " Values:";
	for (auto i:values) std::cout << i << " ";
	std::cout << std::endl;

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	res = A*v;
	for (auto i = 0u; i < A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;

	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;

	// Change matrix
	std::cout << "Changing A..." << std::endl;
	A.decompress();
	A.ref(6,6)=10;
	A.makeCompressed();
	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;
	for (auto i=0u;i<A.nRows;++i)
	{
		for (auto j=0u;j<A.nCols;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;
	inner=A.getInnerIndexes();
	outer=A.getOuterIndexes();
	values=A.getValues();

	std::cout<<" Inner indexes:";
	for (auto i:inner) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Outer indexes:";
	for (auto i:outer) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Values:";
	for (auto i:values) std::cout<<i<<" ";
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v1{1,2,3,4,5,6,7};
	res=A*v1;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
}
	std::cout<<"\n\n***************************************************************\n\n";
{
	std::cout << "Creating matrix A uncompressed - COLUMNMAJOR... " << std::endl;
	auto n=4u;
	auto m=5u;
	Matrix<double,COLUMNMAJOR> A(n,m);
	A.ref(0,0)=1;
	A.ref(3,2)=-1;
	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;
	for (auto i=0u;i<n;++i)
	{
		for (auto j=0u;j<m;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v{1,2,3,4,5};
	std::vector<double> res=A*v;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;

	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	
	
	// Compressed case
	std::cout << "Compressed case..." << std::endl;
	A.makeCompressed();
	for (auto i=0u;i<n;++i)
	{
		for (auto j=0u;j<m;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	auto inner=A.getInnerIndexes();
	auto outer=A.getOuterIndexes();
	auto values=A.getValues();
	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;
	std::cout<<" Inner indexes:";
	for (auto i:inner) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Outer indexes:";
	for (auto i:outer) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Values:";
	for (auto i:values) std::cout<<i<<" ";
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	res=A*v;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;


	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	

	// Change matrix
	std::cout << "Changing A..." << std::endl;
	A.decompress();
	A.ref(6,6)=10;
	A.makeCompressed();
	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;
	for (auto i=0u;i<A.nRows;++i)
	{
		for (auto j=0u;j<A.nCols;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;
	inner=A.getInnerIndexes();
	outer=A.getOuterIndexes();
	values=A.getValues();

	std::cout<<" Inner indexes:";
	for (auto i:inner) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Outer indexes:";
	for (auto i:outer) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Values:";
	for (auto i:values) std::cout<<i<<" ";
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v1{1,2,3,4,5,6,7};
	res=A*v1;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
}
	std::cout<<"\n\n***************************************************************\n\n";
{
	using namespace std::complex_literals;
	std::cout << "Creating matrix A uncompressed - ROWMAJOR with complex values... " << std::endl;
	auto n=4u;
	auto m=5u;
	Matrix<std::complex<double>,ROWMAJOR> A(n,m);
	A.ref(0,0)=1. + 1i;
	A.ref(3,2)=-1.;
	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;

	for (auto i=0u;i<n;++i)
	{
		for (auto j=0u;j<m;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<std::complex<double>> v{1i,2i,3i,4i,5i};
	std::vector<std::complex<double>> res=A*v;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;
	
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;


	// Compressed case
	std::cout << "Compressed case..." << std::endl;
	A.makeCompressed();
	for (auto i=0u;i<n;++i)
	{
		for (auto j=0u;j<m;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	auto inner=A.getInnerIndexes();
	auto outer=A.getOuterIndexes();
	auto values=A.getValues();
	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;

	std::cout<<" Inner indexes:";
	for (auto i:inner) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Outer indexes:";
	for (auto i:outer) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Values:";
	for (auto i:values) std::cout<<i<<" ";
	std::cout<<std::endl;
	

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	res=A*v;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;
	
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;


	// Change matrix
	std::cout << "Changing A..." << std::endl;
	A.decompress();
	A.ref(6,6)=10.;
	A.makeCompressed();
	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;
	for (auto i=0u;i<A.nRows;++i)
	{
		for (auto j=0u;j<A.nCols;++j)
		{
			std::cout<<A(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;
	inner=A.getInnerIndexes();
	outer=A.getOuterIndexes();
	values=A.getValues();

	std::cout<<" Inner indexes:";
	for (auto i:inner) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Outer indexes:";
	for (auto i:outer) std::cout<<i<<" ";
	std::cout<<std::endl;
	std::cout<<" Values:";
	for (auto i:values) std::cout<<i<<" ";
	std::cout<<std::endl;


	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<std::complex<double>> v1{1i,2i,3i,4i,5i,6i,7i};
	res=A*v1;
	for (auto i=0u; i<A.nRows; ++i)
		std::cout << res[i] << " ";
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
}
	std::cout<<"\n\n***************************************************************\n\n";
{
	std::cout << "Reading matrix using Matrix Market - ROWMAJOR case..." << std::endl;
	Matrix<double, ROWMAJOR> A;

	// Read the datafile name from command line
	GetPot cl(argc, argv);
	const char* filename = cl.follow("lnsp_131.mtx", 2, "-f", "--file");

	A.reader(filename);

	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v(A.nCols, 1.0);
	std::vector<double> res;

	auto start = timer::now();

	res=A*v;

	auto end = timer::now();
	
	auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	std::cout << "Elapsed Time: " << elapsed_time << "microseconds" << std::endl;

	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	
	
	// Compressed case
	std::cout << "Compressed case..." << std::endl;
	A.makeCompressed();

	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;

	start = timer::now();

	res=A*v;

	end = timer::now();

	elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	std::cout << "Elapsed Time: " << elapsed_time << "microseconds" << std::endl;
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	
}
	std::cout<<"\n\n***************************************************************\n\n";
{
	std::cout << "Reading matrix using Matrix Market - COLUMNMAJOR case..." << std::endl;
	Matrix<double, COLUMNMAJOR> A;

	// Read the datafile name from command line
	GetPot cl(argc, argv);
	const char* filename = cl.follow("lnsp_131.mtx", 2, "-f", "--file");

	A.reader(filename);

	std::cout<<"nrows= "<<A.nRows<< " " << "ncols= "<<A.nCols<<std::endl;

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;
	std::vector<double> v(A.nCols, 1.0);
	std::vector<double> res;

	auto start = timer::now();

	res=A*v;

	auto end = timer::now();

	auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	std::cout << "Elapsed Time: " << elapsed_time << "microseconds" << std::endl;
	
	
	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	
	
	// Compressed case
	std::cout << "Compressed case..." << std::endl;
	A.makeCompressed();

	std::cout<<"#non zero elements="<<A.nnz()<<std::endl;

	// Matrix-vector multiplication
	std::cout << "Computing A*v..." << std::endl;

	start = timer::now();

	res=A*v;

	end = timer::now();

	elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	std::cout << "Elapsed Time: " << elapsed_time << "microseconds" << std::endl;

	// Norm
	std::cout << "L1 Norm " << A.norm<ONE>() << std::endl;
	std::cout << "Inf norm " << A.norm<INF>() << std::endl;
	std::cout << "Frobenius norm " << A.norm<FROBENIUS>() << std::endl;
	
}
}


