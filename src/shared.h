#ifndef DX_SHARED_
#define DX_SHARED_

#include <vector>
#include <stdexcept>
#include <algorithm>

inline double SQR(const double v){ return v*v; };
inline long double SQR(const long double v){ return v*v; };
inline int SQR(const int v){ return v*v; };
inline double CUB(const double v){ return v*v*v; };
inline long double CUB(const long double v){ return v*v*v; };
inline int CUB(const int v){ return v*v*v; };
inline double QRT(const double v){ return v*v*v*v; };
inline long double QRT(const long double v){ return v*v*v*v; };
inline int QRT(const int v){ return v*v*v*v; };


arma::mat mat_init(const arma::mat& orig);

arma::vec vec_init(const arma::vec& orig);


#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = mat_init(omp_orig) )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
initializer( omp_priv = vec_init(omp_orig) )



// very minimal drop ins for non-existent armadillo matrix and vector of long double
// for use in templated functions
// only defines the members that are actually used in elsym

// contiguous column proxy
class ldproxycol
{
	private:
		long double *ptr;
	public:
		int n_elem; 
		ldproxycol(long double *p, const int n_rows)
		{
			ptr = p;
			n_elem = n_rows;
		}
		
		long double* memptr()
		{
			return ptr;
		}


		template <class V2>
		void operator = (V2& x)
		{
			std::copy(x.memptr(), x.memptr() + x.n_elem, ptr);
		}
};

class ldmat
{
protected:
	std::vector<long double> data;
private:	
	
	// non contiguous row proxy
	class ldrow
	{
	private:
		long double *ptr;
		int nelem, offset;
	public:
		ldrow(long double *p, const int n_rows, const int n_cols)
		{
			ptr = p;
			nelem = n_cols;
			offset = n_rows;
		}
		
		void zeros()
		{
			for(int j=0;j<nelem; j++) *(ptr + j*offset) = 0;
		}
		
		void ones()
		{
			for(int j=0;j<nelem; j++) *(ptr + j*offset) = 1;
		}

	};

public:
	int n_cols, n_rows, n_elem;	

	ldmat(const int n_row, const int n_column, const bool init_zero=true)
	{
		if(n_row<=0 || n_column<= 0)
			throw std::invalid_argument("ldmat cannot be sized zero in any dimension");
		
		n_elem = n_row*n_column;
		
		data = std::vector<long double>(n_elem);
		
		n_cols = n_column;
		n_rows = n_row;
		
		if(init_zero)
			std::fill(data.begin(), data.end(), 0);
	}
	

	void zeros()
	{
		std::fill(data.begin(), data.end(), 0);
	}
	
	void ones()
	{
		std::fill(data.begin(), data.end(), 1);
	}

	long double* memptr()
	{
		return &data[0];
	}
	
	long double* colptr(const int j)
	{
		return &data[j&n_rows];
	}
	
	long double& at(const int i, const int j)
	{
		return data[j*n_rows+i];
	}
	
	ldproxycol col(const int j)
	{
		return ldproxycol(&data[j*n_rows], n_rows);
	}
	
	ldrow row(const int i)
	{
		return ldrow(&data[i],n_rows,n_cols);
	}
};

class ldvec : public ldmat
{
public:
	ldvec(const int n_elements) : ldmat(n_elements, 1, true){}
	
	void operator = (ldproxycol x)
	{
		std::copy(x.memptr(), x.memptr() + x.n_elem, data.begin());
	}

	void operator = (ldvec& x)
	{
		std::copy(x.memptr(), x.memptr() + x.n_elem, data.begin());
	}

	long double& operator[](const int i)
	{
		return data[i];
	}
	
};




#endif