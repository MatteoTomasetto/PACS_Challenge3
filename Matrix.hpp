#ifndef MATRIX_MATRIX_HPP_
#define MATRIX_MATRIX_HPP_

#include <array>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>
#include <exception>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mmio.h"

namespace apsc
{

enum StorageOrder
{
	ROWMAJOR,
	COLUMNMAJOR
};

enum NormType
{
	ONE,
	INF,
	FROBENIUS
};

template <typename Scalar, StorageOrder Order=ROWMAJOR>
class Matrix
{
public:
	using COO = std::array<std::size_t, 2u>; // type to store a pair of indexes

	Matrix()=default; // default constructor

	Matrix(std::size_t nr,std::size_t nc):nRows{nr},nCols{nc}{}; // constructor that takes the number of rows and cols

	void reader(const char*);

	std::size_t size() const {return nRows*nCols;};

	void resize(std::size_t nr,std::size_t nc);

	void reset();

	void makeCompressed();

	void decompress();

	std::size_t nnz() const; // number of non zeros

	bool is_compressed() const {return compressed_;}

	Scalar operator()(std::size_t i,std::size_t j) const; // get the value in position (i,j)

	Scalar operator[](std::size_t i) const; // get the value in position (i,1), useful for matrix (nRows x 1) or (1 x nCols)

	Scalar& ref(std::size_t i,std::size_t j); // setter for the value in position (i,j)

	template<class VectorInputType>
	std::vector<Scalar> operator*(const VectorInputType& v) const; // matrix-vector multiplication

	std::size_t nInner()const {return Order==ROWMAJOR? nRows : nCols;}

	std::size_t nOuter()const {return Order==ROWMAJOR? nCols : nRows;}

	std::size_t nRows{0};

	std::size_t nCols{0};

	auto getInnerIndexes()const {return innerIndexes;}

	auto getOuterIndexes()const {return outerIndexes;}

	auto getValues()const {return values;}

	template<NormType normtype = FROBENIUS >
	double norm() const; 


private:
	// function useful to manage easily ROWMAJOR and COLUMNMAJOR cases in compressed case
	static std::size_t const &  getInner(COO const & c)
	{
		if constexpr (Order==ROWMAJOR)
			return c[0];
		else
			return c[1];
	};

	// function useful to manage easily ROWMAJOR and COLUMNMAJOR cases in compressed case
	static std::size_t const & getOuter(COO const & c)
	{
		if constexpr (Order==ROWMAJOR)
			return c[1];
		else
			return c[0];
	};

	// Comparison to give an order in the uncompressed representation (std::map):
	// if Order == ROWMAJOR then values in the map are ordered by rows;
	// if Order == COLUMNMAJOR then values in the map are ordered by columns.
	static constexpr auto compare = [](COO const &coo1, COO const &coo2)->bool
	{
		return (std::tie(getInner(coo1), getOuter(coo1)) < std::tie(getInner(coo2),getOuter(coo2)));
	};

	// uncompressed representation
	using NCData = std::map<COO, Scalar, decltype(compare)>;
	NCData						ncData{compare};

	// compressed representation
	std::vector<std::size_t>	innerIndexes;
	std::vector<std::size_t>	outerIndexes;
	std::vector<Scalar>			values;

	bool						compressed_{false}; // boolean to understand which representation is used
};

#include "Matrix_imp.hpp"

} // end namespace apsc

#endif /* MATRIX_MATRIX_HPP_ */
