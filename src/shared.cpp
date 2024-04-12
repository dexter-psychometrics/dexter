#include <RcppArmadillo.h>
#include "shared.h"

using namespace arma;


mat mat_init(const mat& orig)
{
	return mat(orig.n_rows, orig.n_cols, fill::zeros);
}

vec vec_init(const vec& orig)
{
	return vec(orig.n_elem, fill::zeros);
}