#ifndef MATRIX_MATRIX_IMP_HPP_
#define MATRIX_MATRIX_IMP_HPP_

template<typename Scalar, StorageOrder Order>
inline void
Matrix<Scalar, Order>::reader(const char* filename)
{
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;   
	int i, *I, *J;
	double *val;

	if ((f = fopen(filename, "r")) == NULL)
		throw std::runtime_error("I cannot open the file with the matrix, sorry");

	if (mm_read_banner(f, &matcode) != 0)
	throw std::runtime_error("Could not process Matrix Market banner, sorry");

	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) )
		throw std::runtime_error("Sorry, this application does not support this Matrix Market type");

	/* find out size of sparse matrix .... */

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
		throw std::runtime_error("Sorry, error reading the matrix");

	/* reseve memory for matrices */

	I = (int *) malloc(nz * sizeof(int));
	J = (int *) malloc(nz * sizeof(int));
	val = (double *) malloc(nz * sizeof(double));

	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	nRows = M;
	nCols = N;    

	for (i=0; i<nz; i++)
	{
		fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
  
		I[i]--;	/* adjust from 1-based to 0-based */
		J[i]--;	/* adjust from 1-based to 0-based */

		ncData.insert( { {I[i], J[i]}, val[i]} );
    }
	if (f != stdin) fclose(f);

	return;
}


template<typename Scalar, StorageOrder Order>
inline void
Matrix<Scalar, Order>::resize (std::size_t nr, std::size_t nc)
{
	if (nr<nRows || nc <nCols)
	{
		throw std::runtime_error("I can only expand a Matrix with resize, sorry");
	}
    
	nRows=nr;
	nCols=nc;
	if (compressed_)
	{
		decompress();
		makeCompressed();
	}
}


template<typename Scalar, StorageOrder Order>
inline void
Matrix<Scalar, Order>::reset ()
{
	nRows=0;
	nCols=0;
	ncData.clear();
	innerIndexes.clear();
	outerIndexes.clear();
	values.clear();
	compressed_=false;
}


template<typename Scalar, StorageOrder Order>
inline void
Matrix<Scalar, Order>::makeCompressed ()
{
	if(compressed_) return;

	if constexpr(Order==ROWMAJOR)
	{
		innerIndexes.resize(nRows+1u);

	}else{

		innerIndexes.resize(nCols+1u);
	}
	
	std::fill(innerIndexes.begin(), innerIndexes.end(), 0u);
	outerIndexes.resize(nnz());
	values.resize(nnz());
	std::size_t counter{0u};
    
	for (const auto & [c,val]:ncData)
	{
		auto inner=getInner(c);
		auto outer=getOuter(c);

		++innerIndexes[inner+1u];
		outerIndexes[counter]=outer;
		values[counter]=val;
		++counter;
	}

	// Fix inner
	for (auto i=1u; i<innerIndexes.size(); ++i)
	{
		innerIndexes[i]+=innerIndexes[i-1u];
	}

	ncData.clear();
	compressed_=true;
	}


template<typename Scalar, StorageOrder Order>
inline void
Matrix<Scalar, Order>::decompress ()
{
	if (!compressed_)return;

	for (auto in=0u; in<nInner(); ++in)
	{
		for(auto k=innerIndexes[in]; k<innerIndexes[in+1]; ++k)
		{
			auto out = outerIndexes[k];

			if constexpr(Order == ROWMAJOR)
				ncData.insert({{in,out},values[k]});

			else
				ncData.insert({{out,in},values[k]});
		}
	}

	innerIndexes.clear();
	outerIndexes.clear();
	values.clear();
	compressed_=false;
}


template<typename Scalar, StorageOrder Order>
inline Scalar
Matrix<Scalar, Order>::operator () (std::size_t i, std::size_t j) const
{
	assert(i<nRows && j<nCols);

	if(compressed_)
	{
		// save the indexes to find easily the element exploiting the compressed form
		std::size_t inner = getInner({i,j});
		std::size_t outer = getOuter({i,j});
		std::size_t offset= innerIndexes[inner];
		std::size_t last  = innerIndexes[inner+1u];
		auto res =std::find(outerIndexes.begin()+offset,outerIndexes.begin()+last,outer);

		if(res==outerIndexes.begin()+last)
		{
			return Scalar{0};
		}

		else
		{
			return values[std::distance(outerIndexes.begin(),res)];
		}
	}

	else
	{
		auto res = ncData.find({i,j});
		if( res == ncData.end() )
			return Scalar{0};
		else
		{
			return res->second;
		}
	}
}


template<typename Scalar, StorageOrder Order>
inline Scalar
Matrix<Scalar, Order>::operator [] (std::size_t i) const
{
	if(nRows == 1)
		return this->operator()(0,i);
	if(nCols == 1)
		return this->operator()(i,0);
	else
		throw std::runtime_error("I cannot use operator[] for matrices having nRows!=1 and nCols!=1, sorry");
}


template<typename Scalar, StorageOrder Order>
inline Scalar &
Matrix<Scalar, Order>::ref (std::size_t i, std::size_t j)
{
	if(compressed_)
	{
		assert(i<nRows && j<nCols);

		// save the indexes to find easily the element exploiting the compressed form
		std::size_t inner = getInner({i,j});
		std::size_t outer = getOuter({i,j});
		std::size_t offset= innerIndexes[inner];
		std::size_t last  = innerIndexes[inner+1u];
		auto res =std::find(outerIndexes.begin()+offset,outerIndexes.begin()+last,outer);

		if(res==outerIndexes.begin()+last)
		{
			throw std::runtime_error("Cannot add element in compressed state");
		}

		else
		{
			return values[std::distance(outerIndexes.begin(),res)];
		}
	}

	else
	{
		auto res = ncData.find({i,j});

		if( res == ncData.end() )
		{
			nRows=std::max(nRows,i+1u);
			nCols=std::max(nCols,j+1u);
			auto pos=ncData.insert({{i,j},Scalar{0}}); // if {i,j} not found, add it to the map with value = 0
			return pos.first->second;
		}

		else
		{
			return res->second;
		}
	}
}


template<typename Scalar, StorageOrder Order>
template<class VectorInputType> 
inline std::vector<Scalar> 
Matrix<Scalar, Order>::operator* (const VectorInputType& v) const
{	
	if( !std::is_convertible_v<Scalar, decltype(v[0])> )
		throw std::runtime_error("Matrix and vector elements do not match the type, I'm sorry");

	if(v.size() != nCols)
		throw std::runtime_error("Matrix and vector dimensions do not match, I'm sorry");

	std::vector<Scalar> res(nRows, 0.0); // initialize the solution

	if(compressed_)
	{
		for (auto in=0u; in<nInner(); ++in)
		{
			for(auto k=innerIndexes[in]; k<innerIndexes[in+1]; ++k)
			{
				if constexpr (Order == ROWMAJOR)
				{
					// in --> row index
					// outerIndexes[k] --> col index
					// values[k] --> value
					res[in] += values[k] * v[outerIndexes[k]];
				}

				else
				{
					// in --> col index
					// outerIndexes[k] --> row index
					// values[k] --> value
					res[outerIndexes[k]] += values[k] * v[in];
				}
			}
		}
	}


	else
	{
		for (const auto & [c, val] : ncData)
		{
			// c[0] --> row index
			// c[1] --> col index
			// val --> value
			res[c[0]] += val * v[c[1]];
		}
	}


	return res;
}


template<typename Scalar, StorageOrder Order>
inline std::size_t
Matrix<Scalar, Order>::nnz () const
{
	if(compressed_)
		return outerIndexes.size();
	else
		return ncData.size();
}


template<typename Scalar, StorageOrder Order>
template<NormType normtype>
inline double
Matrix<Scalar, Order>::norm() const
{
	if constexpr(normtype == ONE)
	{
		std::vector<double> ColSums(nCols, 0.0);
	
		if(compressed_)
		{
			for (auto in=0u; in<nInner(); ++in)
			{
				for(auto k=innerIndexes[in]; k<innerIndexes[in+1]; ++k)
				{
					if constexpr(Order == ROWMAJOR)
					{
						// in --> row index
						// outerIndexes[k] --> col index
						// values[k] --> values
						ColSums[outerIndexes[k]] += std::abs(values[k]);
					}

					else
					{
						// in --> col index
						// outerIndexes[k] --> row index
						// values[k] --> values
						ColSums[in] += std::abs(values[k]);
					}
				}
			}
		}

		else
		{
			for (const auto & [c,val] : ncData)
				ColSums[c[1]] += std::abs(val);
		}

		return *std::max_element(ColSums.begin(), ColSums.end());

	}

	else if constexpr (normtype == INF)
	{
		std::vector<double> RowSums(nCols, 0.0);

		if(compressed_)
		{
			for (auto in=0u; in<nInner(); ++in)
			{
				for(auto k=innerIndexes[in]; k<innerIndexes[in+1]; ++k)
				{
					if constexpr(Order == ROWMAJOR)
					{
						// in --> row index
						// outerIndexes[k] --> col index
						// values[k] --> values
						RowSums[in] += std::abs(values[k]);
					}

					else
					{
						// in --> col index
						// outerIndexes[k] --> row index
						// values[k] --> values
						RowSums[outerIndexes[k]] += std::abs(values[k]);
					}
				}
			}
		}

		else
		{
			for (const auto & [c,val] : ncData)
				RowSums[c[0]] += std::abs(val);
		}

		return *std::max_element(RowSums.begin(), RowSums.end());
	}

	else if constexpr(normtype == FROBENIUS)
	{
		double Sum{0.0};
	
		if(compressed_)
		{
			for (std::size_t i = 0; i < values.size(); ++i)
			{
				Sum += std::abs(values[i])*std::abs(values[i]);
			}
		}

		else
		{
			for (const auto & [c,val] : ncData)
				Sum += std::abs(val)*std::abs(val);
		}


		return std::sqrt(Sum);
	}

	else
		throw std::runtime_error("Norm type not implemented yet, sorry");
}





#endif
