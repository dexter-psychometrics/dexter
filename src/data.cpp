#include <RcppArmadillo.h>
#include <stack>
using namespace Rcpp;


// faster factor creation
// http://gallery.rcpp.org/articles/fast-factor-generation

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x, bool as_int ) 
{
    Vector<RTYPE> levs = sort_unique(x);
    IntegerVector out = match(x, levs);
	if(!as_int)
	{
		out.attr("levels") = as<CharacterVector>(levs);
		out.attr("class") = "factor";
	}
    return out;
}

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x, const Vector<RTYPE>& levs, bool as_int ) 
{
    IntegerVector out = match(x, levs);
	if(!as_int)
	{
		out.attr("levels") = as<CharacterVector>(levs);
		out.attr("class") = "factor";
	}
    return out;
}

// [[Rcpp::export]]
SEXP fast_factor( SEXP x, bool as_int) 
{
    switch( TYPEOF(x) ) {
    case INTSXP: return fast_factor_template<INTSXP>(x, as_int);
    case REALSXP: return fast_factor_template<REALSXP>(x, as_int);
    case STRSXP: return fast_factor_template<STRSXP>(x, as_int);
    }
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP fast_factor_lev( SEXP x, SEXP levs, bool as_int) 
{
    switch( TYPEOF(x) ) {
    case INTSXP: return fast_factor_template<INTSXP>(x, levs, as_int);
    case REALSXP: return fast_factor_template<REALSXP>(x, levs, as_int);
    case STRSXP: return fast_factor_template<STRSXP>(x, levs, as_int);
    }
    return R_NilValue;
}


// [[Rcpp::export]]
std::string ppoint(SEXP x)
{
	return tfm::format("%p", x);
}


// assumed there are no double responses
// to do: can be faster if we unroll the loop
// [[Rcpp::export]]
void fill_resp_matrix(const IntegerVector& person_id, const IntegerVector& item_id, const IntegerVector& item_score, arma::imat& out)
{
	const int n = person_id.length();

#pragma omp parallel for
	for(int i=0; i<n; i++)
	{
		out.at(person_id[i]-1, item_id[i]-1) = item_score[i];
	}
}

// [[Rcpp::export]]
IntegerVector ds_connected_groups(const IntegerMatrix& a)
{
	const int n = a.ncol();	
	IntegerVector group(n);
	int g=0;
	std::stack<int> st;

	for(int j=0;j<n;j++)
	{
		if(group[j]==0)
		{
			st.push(j);
			group[j] = ++g;
			while(!st.empty())
			{
				int s = st.top();
				st.pop();
				for(int i=0;i<n;i++)
				{
					if(a(i,s)>0 && group[i]==0)
					{
						group[i]=g;
						st.push(i);
					}
				}	
			}
		}
	}
	return group;
}

// needs two groups with values 1 and 2 in group_id
// already tried parallel, practically no gain
// [[Rcpp::export]]
std::vector<int> unequal_categories_C(const IntegerVector& group_id, const IntegerVector& item_id, const IntegerVector& item_score, const int nit, const int max_score)
{	
	std::vector<int> scratch((nit+1)*(max_score+1),0);
	
	std::vector<int> out;
	out.reserve(nit);
	
	const int nr = item_id.length();	
	const int isz = max_score+1;
	
	for(int i=0;i<nr;i++)
	{
		if(scratch[item_id[i] * isz + item_score[i]] == group_id[i])
			scratch[item_id[i] * isz + item_score[i]] = 3;
		else if(scratch[item_id[i] * isz + item_score[i]] == 0)
			scratch[item_id[i] * isz + item_score[i]] = 3 - group_id[i]; //2->1, 1->2
	}	
	
	for(int i=1; i<=nit; i++)
		for(int s=0; s<=max_score; s++)
		{
			if(scratch[i*isz+s] == 1 || scratch[i*isz+s] == 2)
			{
				out.push_back(i);
				break;
			}
		}
	out.shrink_to_fit();
	return out;	
}


/* ******************************************************************************************* 
* Booklets, sumscores and designs in a mutate and summarize variant for data.frame columns
*
* all columns must be sorted by person_id
* item_id and booklet_id must be factors
* non-const input columns are changed in place
*
* ********************************************************************************************/



/* ********************************************************************************************
/ the vector booklet_id is always overwritten (if not desirable, use a copy)
/ In case of merge: all original booklet id's are lost (not mapped at all)
/ In case of no merge: returned df map_booklet contains a map of new to old booklet id's
/ 	booklets with different id's that have equal items are maintained as different booklets
/
/
******************************************************************************************** */

template<bool merge>
List make_booklets_tmpl(const IntegerVector& person_id, const IntegerVector& item_id, const IntegerVector& item_score, IntegerVector& booklet_id, IntegerVector& booklet_score)
{
	const int nit = as<CharacterVector>(item_id.attr("levels")).length();
	const int nr = item_id.length();
		
	int nbk = 1;
	int nitbk = 0;
	
	typedef std::pair<std::vector<bool>, int> key;	
	
	struct key_hash
	{
		std::size_t operator()(const key& k) const
		{
			// hash bool vector xor int
			auto hash1 = std::hash<std::vector<bool>>{}(k.first); 
			auto hash2 = std::hash<int>{}(k.second); 
			return hash1 ^ hash2; 
		}	
	};
		
	
	std::vector<bool> bk(nit+1);
	std::unordered_map<key, int, key_hash> booklets;
	
	std::fill(bk.begin(), bk.end(), false);
	
	bool overlap = false;
	int person_start = 0;
	int ss = item_score[0];
	bk[item_id[0]] = true;	
	
	for(int r=1; r<nr; r++) 
	{
		if(person_id[r] != person_id[r-1] || !(booklet_id[r] == booklet_id[r-1] || merge))
		{
			const auto& ret = booklets.insert( std::make_pair(std::make_pair(bk, merge ? 1 : booklet_id[r-1]), nbk));
			
			std::fill(booklet_id.begin()+person_start, booklet_id.begin() + r, ret.first->second);
			std::fill(booklet_score.begin()+person_start, booklet_score.begin() + r, ss);
			
			if (ret.second) // did not already exist
			{
				nbk++;
				for(int i=1; i<=nit; i++) //1 indexed items
					if(bk[i])
						nitbk += 1;
			}
			
			std::fill(bk.begin(), bk.end(), false);
			ss = 0;	
			person_start = r;
		}
		
		if(bk[item_id[r]])
			overlap = true;
		
		bk[item_id[r]] = true;	
		ss += item_score[r];
	}
	if(overlap)
	{
		stop("at least one person has answered at least one item more than once, this is not allowed");
	}
	// finally
	const auto& ret = booklets.insert( std::make_pair(std::make_pair(bk, merge ? 1 : booklet_id[nr-1]), nbk));
	std::fill(booklet_id.begin()+person_start, booklet_id.end(), ret.first->second);
	std::fill(booklet_score.begin()+person_start, booklet_score.end(), ss);
	
	if (ret.second)
	{
		nbk++;
		for(int i=1; i<=nit; i++) //1 indexed items
			if(bk[i])
				nitbk += 1;
	}


	// design is collateral (if not item_position and other design columns necessary)
	IntegerVector dbooklet(nitbk), ditem(nitbk);
	IntegerVector map_booklet_id(nbk-1), map_orig_booklet(nbk-1);	
	
	int indx = 0;
	int bkindx = 0;
	for(auto& iter: booklets )
	{
		const auto& bki =  iter.first.first;
		int bnum = iter.second;
		if(!merge)
		{
			map_orig_booklet[bkindx] = iter.first.second;
			map_booklet_id[bkindx++] = bnum;
		}
		for(int i=1;i<=nit;i++)
		{
			if(bki[i])			
			{				
				dbooklet[indx] = bnum;
				ditem[indx++] = i; //1-indexed item_id from r
			}
		}		
	}
	
	ditem.attr("levels") = as<CharacterVector>(item_id.attr("levels"));
	ditem.attr("class") = "factor";	

	if(Rf_isFactor(booklet_id))
	{
		if(!merge)
		{
			map_orig_booklet.attr("levels") = as<CharacterVector>(booklet_id.attr("levels"));
			map_orig_booklet.attr("class") = "factor";
		}	
		booklet_id.attr("levels") = R_NilValue;
		booklet_id.attr("class") = "integer";	
	}
	//map booklet is only relevant in case there is no merge
	return List::create(
		Named("design") = DataFrame::create( 
			Named("booklet_id") = dbooklet,
			Named("item_id") = ditem),
		Named("map_booklet") = DataFrame::create( 	
			Named("booklet_id") = map_booklet_id,
			Named("org_booklet_id") = map_orig_booklet)
		);

}









template<bool merge>
List make_booklets_summed_tmpl(IntegerVector& person_id, IntegerVector& booklet_id, IntegerVector& item_id, IntegerVector& item_score)
{
	const int nit = as<CharacterVector>(item_id.attr("levels")).length();
	const int nr = item_id.length();
	
	typedef std::pair<std::vector<bool>, int> key;	
	
	struct key_hash
	{
		std::size_t operator()(const key& k) const
		{
			// hash bool vector xor int
			auto hash1 = std::hash<std::vector<bool>>{}(k.first); 
			auto hash2 = std::hash<int>{}(k.second); 
			return hash1 ^ hash2; 
		}	
	};
	
	std::vector<bool> bk(nit+1, false);
	std::unordered_map<key, int, key_hash> booklets;	

	int np = 0;
	int nitbk = 0;
	int nbk = 1;
	bool overlap = false;
	int ss = item_score[0];
	bk[item_id[0]] = true;	

	for(int r=1; r<nr; r++)
	{
		if(person_id[r] != person_id[r-1] || !(booklet_id[r] == booklet_id[r-1] || merge))
		{
			const auto& ret = booklets.insert( std::make_pair(std::make_pair(bk, merge ? 1 : booklet_id[r-1]), nbk));
			
			booklet_id[np] = ret.first->second;
			item_score[np] = ss;
			person_id[np] = person_id[r-1];
			item_id[np] = r; //permutation index, 1-indexed for R
			if (ret.second) // did not already exist
			{
				nbk++;
				for(int i=1; i<=nit; i++) //1 indexed items
					if(bk[i])
						nitbk += 1;
			}
			std::fill(bk.begin(), bk.end(), false);
			ss = 0;	
			np++;
		}
		if(bk[item_id[r]])
			overlap = true;
		
		bk[item_id[r]] = true;	
		ss += item_score[r];
	}
	if(overlap)
	{
		stop("at least one person has answered at least one item more than once, this is not allowed");
	}

	// finally
	const auto& ret = booklets.insert( std::make_pair(std::make_pair(bk, merge ? 1 : booklet_id[nr-1]), nbk));
	booklet_id[np] = ret.first->second;
	item_score[np] = ss;
	person_id[np] = person_id[nr-1];
	item_id[np] = nr;
	
	if (ret.second)
	{
		nbk++;
		for(int i=1; i<=nit; i++) //1 indexed items
			if(bk[i])
				nitbk += 1;
	}


	// design is collateral (if not item_position and other design columns necessary)
	IntegerVector dbooklet(nitbk), ditem(nitbk);
	IntegerVector map_booklet_id(nbk-1), map_orig_booklet(nbk-1);	
	
	int indx = 0;
	int bkindx = 0;

	for(auto& iter: booklets )
	{
		const auto& bki =  iter.first.first;
		int bnum = iter.second;
		if(!merge)
		{
			map_orig_booklet[bkindx] = iter.first.second;
			map_booklet_id[bkindx++] = bnum;
		}
		for(int i=1;i<=nit;i++)
		{
			if(bki[i])			
			{				
				dbooklet[indx] = bnum;
				ditem[indx++] = i; //1-indexed item_id from r
			}
		}		
	}
	
	
	ditem.attr("levels") = as<CharacterVector>(item_id.attr("levels"));
	ditem.attr("class") = "factor";	
	
	item_id.attr("levels") = R_NilValue;
	item_id.attr("class") = "integer";

	if(Rf_isFactor(booklet_id))
	{
		if(!merge)
		{
			map_orig_booklet.attr("levels") = as<CharacterVector>(booklet_id.attr("levels"));
			map_orig_booklet.attr("class") = "factor";
		}	
		booklet_id.attr("levels") = R_NilValue;
		booklet_id.attr("class") = "integer";
	}

	return List::create(
		Named("np") = np+1, 
		Named("design") = DataFrame::create( 
			Named("booklet_id") = dbooklet,
			Named("item_id") = ditem),
		Named("map_booklet") = DataFrame::create( 	
			Named("booklet_id") = map_booklet_id,
			Named("org_booklet_id") = map_orig_booklet)
		);
}









// [[Rcpp::export]]
List make_booklets(const IntegerVector& person_id, const IntegerVector& item_id, const IntegerVector& item_score, IntegerVector& booklet_id, IntegerVector& booklet_score, const bool merged)
{
	if(merged)
		return make_booklets_tmpl<true>(person_id, item_id, item_score, booklet_id, booklet_score);
		
	return make_booklets_tmpl<false>(person_id, item_id, item_score, booklet_id, booklet_score);
}


// [[Rcpp::export]]
List make_booklets_summed(IntegerVector& person_id, IntegerVector& booklet_id, IntegerVector& item_id, IntegerVector& item_score, const bool merged)
{
	// template argument is !merged
	if(merged)
		return make_booklets_summed_tmpl<true>(person_id, booklet_id, item_id, item_score);
		
	return make_booklets_summed_tmpl<false>(person_id, booklet_id, item_id, item_score);
}



/* ******************************************************************************************* 
* Booklets, sumscores and designs in a mutate and summarize variant for integer matrices
*
* rows=persons, columns=items, supplied as IntegerVector
* outputs x and design
*
* ********************************************************************************************/

// parallel will only make this slower because of critical insert in map

// [[Rcpp::export]]
List make_booklets_summed_matrix(const IntegerVector& mtx, const int ncol, const int nrow)
{
	// nbk start with 1 to make booklet_id's nicer in R
	int nbk = 1, nitbk = 0, ss;
	
	typedef std::vector<bool> key;		
	// unordered map probably faster
	key bk(ncol);
	std::map<key,int> booklets;
	std::pair<std::map<key,int>::iterator, bool> ret; 
	
	IntegerVector booklet_id(nrow), booklet_score(nrow);
	
	for(int p=0;p<nrow;p++)
	{
		std::fill(bk.begin(), bk.end(), false);
		ss = 0;
		for(int i=0; i< ncol; i++)
		{
			if(!IntegerVector::is_na(mtx[i*nrow+p]))
			{
				bk[i] = true;
				ss += mtx[i*nrow+p];
			}
		}
		ret = booklets.insert( std::pair<key,int>(bk, nbk));
		if (ret.second) // did not already exist
		{
			nbk++;
			for(int i=0; i<ncol; i++)
				if(bk[i])
					nitbk++;
		}
		booklet_id[p] = ret.first->second;
		booklet_score[p] = ss;
	}
	//finally not necessary
	
	IntegerVector dbooklet(nitbk), ditem(nitbk);
	
	int indx = 0;

	for(std::map<key,int>::iterator iter = booklets.begin(); iter != booklets.end(); ++iter)
	{
		bk =  iter->first;
		for(int i=0;i<ncol;i++)
		{
			if(bk[i])						
			{				
				dbooklet[indx] = iter->second;
				ditem[indx++] = i+1; //1-indexed item_id
			}
		}		
	}
	
	
	return List::create(Named("x") = DataFrame::create(Named("booklet_id") = booklet_id, Named("booklet_score") = booklet_score),
						Named("design") = DataFrame::create(Named("booklet_id") = dbooklet, Named("item_id") = ditem));

}


// equivalent to gather %>% make_booklets


// [[Rcpp::export]]
List make_booklets_matrix(const IntegerVector& mtx, const int ncol, const int nrow)
{
	// nbk start with 1 to make booklet_id's nicer in R
	int nbk = 1, nitbk = 0, ss, i_row=0, person_start=0;
	const int sz = mtx.length();
	
	typedef std::vector<bool> key;		

	key bk(ncol);
	std::unordered_map<key,int> booklets;
	
	// matrix typically contains NA's
	// choice is either two passes or create oversized vectors and shrink later or create smaller vectors and expand
	// from experiments, this is slightly faster than reserve and push_back (probably because std::fill is very fast)
	// arma vec maybe more convenient, also doesn't init to 0
	std::vector<int> booklet_id(sz), item_id(sz), booklet_score(sz), item_score(sz), person_id(sz); 
	
	for(int p=0;p<nrow;p++)
	{
		std::fill(bk.begin(), bk.end(), false);
		ss = 0;
		for(int i=0; i< ncol; i++)
		{
			int indx = i*nrow+p;
			if(!IntegerVector::is_na(mtx[indx]))
			{
				bk[i] = true;
				ss += mtx[indx];
				item_score[i_row] = mtx[indx];
				item_id[i_row++] = i+1; // 1-indexed
			}
		}
		const auto& ret = booklets.insert( std::pair<key,int>(bk, nbk));
		if (ret.second) // did not already exist
		{
			nbk++;
			for(int i=0; i<ncol; i++)
				if(bk[i])
					nitbk++;
		}
		std::fill(booklet_id.begin() + person_start, booklet_id.begin() + i_row, ret.first->second);
		std::fill(booklet_score.begin() + person_start, booklet_score.begin() + i_row, ss);
		std::fill(person_id.begin() + person_start, person_id.begin() + i_row, p+1);
		person_start = i_row;
	}
	//finally not necessary
	
	// final size
	booklet_id.resize(i_row);
	item_id.resize(i_row);
	item_score.resize(i_row);
	booklet_score.resize(i_row);
	person_id.resize(i_row);
	
	
	IntegerVector dbooklet(nitbk), ditem(nitbk);
	
	int indx = 0;

	for(auto& iter: booklets)
	{
		auto& bki = iter.first;
		for(int i=0;i<ncol;i++)
		{
			if(bki[i])						
			{				
				dbooklet[indx] = iter.second;
				ditem[indx++] = i+1; //1-indexed item_id
			}
		}		
	}
	
	
	return List::create(
				Named("x") = DataFrame::create(	Named("person_id") = person_id, Named("booklet_id") = booklet_id, Named("item_id") = item_id, 
												Named("item_score") = item_score, Named("booklet_score") = booklet_score),
				Named("design") = DataFrame::create(Named("booklet_id") = dbooklet, Named("item_id") = ditem));

}






/* *****************************************************************************************
* summarise and mutate for trusted datasets
* only for datasets ordered by person and booklet
* equivalent to:
* group_by(person_id, boooklet_id) %>% mutate(booklet_score=sum(item_score))
* or
* group_by(person_id, boooklet_id) %>% summarise(booklet_score=sum(item_score))
*
* ******************************************************************************************/


// system.time({y=y %>% group_by(person_id,booklet_id) %>% mutate(booklet_score=sum(item_score)) %>% ungroup()})
// system.time({ss2 = mutate_booklet_score(y$person_id, y$booklet_id, y$item_score)})
// all(ss2==y$booklet_score)
// about 2.95 vs .04 sec for pisa
// [[Rcpp::export]]
IntegerVector mutate_booklet_score(const IntegerVector& person_id, const IntegerVector& booklet_id, const IntegerVector& item_score)
{
	const int nr = person_id.length();
	int booklet = booklet_id[0];
	int person = person_id[0];
	int start = 0;
	int ss=0;

	IntegerVector booklet_score(nr);
	
	for(int r=0;r<nr;r++)
	{
		if(person != person_id[r] || booklet != booklet_id[r])
		{
			std::fill(booklet_score.begin() + start, booklet_score.begin() + r, ss);
			start = r;
			ss = 0;
			person = person_id[r];
			booklet = booklet_id[r];
		}
		ss += item_score[r];
	}
	//finally
	std::fill(booklet_score.begin() + start, booklet_score.end(), ss);
	
	return booklet_score;
}




// item_id becomes permutation index
// item_score becomes booklet_score
// returns new number of rows, shrinking must be done in R

// [[Rcpp::export]]
int summarise_booklet_score(IntegerVector& person_id, IntegerVector& booklet_id, IntegerVector& item_id, IntegerVector& item_score)
{
	const int nr = person_id.length();
	int booklet = booklet_id[0];
	int person = person_id[0];
	int ss=0;
	int np = 0;
	
	item_id.attr("levels") = R_NilValue;
	item_id.attr("class") = "integer";


	for(int r=0; r< nr; r++)
	{
		if(person != person_id[r] || booklet != booklet_id[r])
		{
			booklet_id[np] = booklet_id[r-1];
			item_score[np] = ss;
			person_id[np] = person;
			item_id[np] = r; //permutation index, 1-indexed for R
			np++;
			ss=0;
			person = person_id[r];
			booklet = booklet_id[r];
		}
		ss += item_score[r];
	}
	//finally
	booklet_id[np] = booklet_id[nr-1];
	item_score[np] = ss;
	person_id[np] = person;
	item_id[np] = nr;
	
	return np+1;
}



// merge over persons for a person-booklet ordered dataset,  sumscore not updated (not useful for summarised)
// [[Rcpp::export]]
DataFrame merge_booklets(IntegerVector& booklet_id,  const IntegerVector& person_id, const IntegerVector ds_booklet_id, const int maxb)
{
	const int n = booklet_id.length();
	typedef std::vector<bool> key;		
	std::unordered_map<key, int> bkcombi;
	
	key books(maxb+1, false); // plus one for 1 indexed factors in R
	
	books[booklet_id[0]] = true; // do the first iter
	int start = 0, nb=1;
	
	for(int i=1; i<n; i++)
	{
		if(person_id[i] != person_id[i-1])
		{
			const auto& ret = bkcombi.insert(std::make_pair(books,nb)); 
			if(ret.second)
				nb++;
			std::fill(booklet_id.begin() + start, booklet_id.begin() + i, ret.first->second); 
			std::fill(books.begin(), books.end(), false);
			start = i;
		}
		books[booklet_id[i]] = true;
	}
	// finalize
	std::fill(booklet_id.begin() + start, booklet_id.end(), bkcombi.insert(std::make_pair(books,nb)).first->second);
	
	
	// need first pass to determine size to avoid shrink
	int sz = 0;
	for(std::unordered_map<key, int>::iterator iter = bkcombi.begin(); iter != bkcombi.end(); ++iter)
	{
		const auto& k = iter->first;
		for(int i=1; i<=maxb; i++)
			if(k[i])
				sz++;
	}
	
	IntegerVector bk_new(sz), bk_old(sz);
	int indx = 0;

	for(std::unordered_map<key, int>::iterator iter = bkcombi.begin(); iter != bkcombi.end(); ++iter)
	{
		const auto& k = iter->first;
		int v = iter->second;
		for(int i=1; i<=maxb; i++)
			if(k[i])
			{
				bk_new[indx] = v;
				bk_old[indx] = i;
				indx++;
			}
	}

	bk_old.attr("class") = "factor";
	bk_old.attr("levels") = booklet_id.attr("levels");
	
	return DataFrame::create(Named("booklet_id") = bk_new, Named("old_booklet_id") = bk_old);
} 





// get the design for a trusted dataset (e.g. from database but mutilated by an unsafe predicate), no side effects

// system.time({ds=distinct(y,booklet_id,item_id) %>% arrange(booklet_id,item_id)})
// system.time({ds2=get_design_C(y$booklet_id,y$item_id)})
// all(ds$booklet_id == ds2$booklet_id)
// all(ds$item_id == ds2$item_id)
// todo: check timing for a ridiculously large item*booklet combination like stex
// about .23 vs .015 for pisa
// [[Rcpp::export]]
DataFrame get_design_C(const IntegerVector& booklet_id, const IntegerVector& item_id)
{
	// hell will break loose if booklet is not a factor
	const int nb = as<CharacterVector>(booklet_id.attr("levels")).length();
	const int nit = as<CharacterVector>(item_id.attr("levels")).length();
	const int nn = nit * nb;
	const int n = booklet_id.length();
	int nbi = 0;
	
	std::vector<bool> seen(nn, false); // bool vector is not thread safe, do not use openmp
	
	for(int i=0;i<n;i++)
		seen[(booklet_id[i]-1) * nit + (item_id[i]-1)] = true;
	
	for(int i=0; i<nn; i++)
		if(seen[i])
			nbi++;
	
	IntegerVector dbk(nbi), dit(nbi);
	
	int idx = 0;
	for(int b=0; b<nb; b++)
		for(int i=0; i< nit; i++)
			if(seen[b*nit+i])
			{
				dbk[idx] = b+1;
				dit[idx++] = i+1;
			}

	dbk.attr("levels") = as<CharacterVector>(booklet_id.attr("levels"));
	dbk.attr("class") = "factor";
	dit.attr("levels") = as<CharacterVector>(item_id.attr("levels"));
	dit.attr("class") = "factor";
	
	return DataFrame::create(Named("booklet_id") = dbk, Named("item_id") = dit);

}
// [[Rcpp::export]]
List polytomize_C(IntegerVector& booklet_id, IntegerVector& person_id, IntegerVector& item_prop, IntegerVector& item_score, IntegerVector& booklet_score, const int nlev, const int nb)
{
	const int n = booklet_id.length();
	std::vector<int> pscore(nlev+1, 0); // 1-indexed 
	std::vector<bool> pobs(nlev+1, false);
	
	std::vector<bool> bseen(nb+1, false), design(nb*nlev);
	
	int indx = 0;
	pscore[item_prop[0]] = item_score[0];
	pobs[item_prop[0]] = true;
	
	for(int i=1; i<n; i++)
	{
		if(person_id[i] != person_id[i-1] || booklet_id[i] != booklet_id[i-1])
		{
			for(int p=1; p<=nlev; p++) if(pobs[p])
			{
				booklet_score[indx] = booklet_score[i-1];
				person_id[indx] = person_id[i-1];
				booklet_id[indx] = booklet_id[i-1];
				item_score[indx] = pscore[p];
				item_prop[indx] = p;		
				indx++;
			}
			if(!bseen[booklet_id[i]])
			{
				for(int p=1; p<=nlev; p++) if(pobs[p])
				{
					design[(booklet_id[i]-1) * nlev + p - 1] = true;
				}
				bseen[booklet_id[i]] = true;
			}
			
			std::fill(pscore.begin(),pscore.end(),0);
			std::fill(pobs.begin(),pobs.end(),false);		
		}
		pscore[item_prop[i]] += item_score[i];
		pobs[item_prop[i]] = true;
	}
	// finally
	for(int p=1; p<=nlev; p++) if(pobs[p])
	{
		booklet_score[indx] = booklet_score[n-1];
		person_id[indx] = person_id[n-1];
		booklet_id[indx] = booklet_id[n-1];
		item_score[indx] = pscore[p];
		item_prop[indx] = p;		
		indx++;
	}
	if(!bseen[booklet_id[n-1]])
		for(int p=1; p<=nlev; p++) if(pobs[p])
			design[(booklet_id[n-1]-1) * nlev + p - 1] = true;
	// end finally
	
	int len_ds = 0;
	for(int i=0; i<nb*nlev; i++)
		if(design[i])
			len_ds++;
			
	IntegerVector ds_booklet_id(len_ds), ds_item_id(len_ds);
	
	int ds_indx = 0;
	for(int b=0; b<nb; b++)
		for(int i=0; i<nlev; i++)
			if(design[b*nlev+i])
			{
				ds_booklet_id[ds_indx] = b+1;
				ds_item_id[ds_indx] = i+1;
				ds_indx++;
			}
	
	ds_booklet_id.attr("class") = "factor";
	ds_booklet_id.attr("levels") = as<CharacterVector>(booklet_id.attr("levels"));
	ds_item_id.attr("class") = "factor";
	ds_item_id.attr("levels") = as<CharacterVector>(item_prop.attr("levels"));
	
	return List::create(Named("design") = DataFrame::create(Named("booklet_id") = ds_booklet_id, Named("item_id") = ds_item_id),
						Named("n") = indx);
}






// http://gallery.rcpp.org/articles/hierarchical-risk-parity/
// https://computing.llnl.gov/tutorials/openMP/

// utility function to check if data is safe to use with above functions

// [[Rcpp::export]]
bool is_person_booklet_sorted(const IntegerVector& booklet_id, const IntegerVector& person_id)
{
	const int n = booklet_id.length();
	std::atomic<bool> sorted(true);
	// prefer parallel time worst case over fast return if not sorted
#pragma omp parallel for
	for(int i=1; i<n; i++)
	{
		if((booklet_id[i] < booklet_id[i-1] && person_id[i] == person_id[i-1]) ||
		   (person_id[i] < person_id[i-1]))
		   sorted = false; 
	}
	return sorted;
}


// 0.52 for size 1e6 * 100 (if loop is necessary)
// item_id must be shrunk to relevant items
// item_id must be ordered
// [[Rcpp::export]]
bool parms_is_superset_matrix(const IntegerMatrix& x, const IntegerVector& item_id, const IntegerVector& item_score, const int maxs)
{
	const int nit = x.ncol();
	const int np = x.nrow();
	const int ns = item_score.length();
	
	if(ns == (maxs * nit + nit)) // we include zero category, if all slots are filled return true
		return true;

	std::vector<bool> scores(maxs * nit + nit, false);
	
	// item scores must be smaller than maxs (checked in R)
	for(int i=0; i< ns; i++)
		scores[(item_id[i]-1) * maxs + item_score[i]] = true;
	
	
	std::atomic<bool> is_super(true);
	
#pragma omp parallel for
	for(int i=0; i< nit; i++)
	{
		for(int j=0; j<np; j++)
			if(x(j,i) > 0 && !scores[i*maxs + x(j,i)])	//NA integer is smaller than zero so this works, zero category is guaranteed to be in parms so can be ignored
				is_super = false;
	}
	
	return is_super;
}




/* ************************************************************************************************************************
* Basic sufficient statistics
*
* ssIS = x %>% 
*    count(.data$item_id, .data$item_score, name='sufI') %>%
*    arrange(.data$item_id, .data$item_score)
*
*
*  plt = x %>% 
*  		group_by(.data$booklet_id, .data$booklet_score, .data$item_id) %>% 
*   	summarise(meanScore=mean(.data$item_score), N=n()) %>% 
*   	ungroup()
*
* *************************************************************************************************************************/



// [[Rcpp::export]]
List suf_stats_nrm(const IntegerVector& booklet_id, const IntegerVector& booklet_score, const IntegerVector& item_id, const IntegerVector& item_score, const int nit, const int max_score)
{
	typedef std::tuple<int, int, int> key;
	typedef std::pair<int, int> val;
	
	const int n = booklet_id.length();	
	
	// prepare map for plt
	struct int3_hash
	{
		std::size_t operator()(const key& k) const
		{
			// booklet, booklet_score, item
			return std::hash<int>()(std::get<0>(k) * 128 + std::get<1>(k)  + std::get<2>(k) * 8192);
		}	
	};

	std::unordered_map<key, val, int3_hash> plt;

	// reserve space for ssIS
	const int nscore = max_score + 1;	
	const int ssIS_len = nit * nscore;
	std::vector<int> ssI(ssIS_len, 0), ssi_item, ssi_item_score;
	ssi_item.reserve(ssIS_len);
	ssi_item_score.reserve(ssIS_len);
	
	// main loop
	for(int i=0; i<n; i++)
	{
		// plt
		auto& s = plt[std::forward_as_tuple(booklet_id[i], booklet_score[i], item_id[i])]; 
		std::get<0>(s) += item_score[i];
		std::get<1>(s)++;
		// ssI
		ssI[(item_id[i]-1) * nscore + item_score[i]]++;
	}
	
	// convert map to data.frame for plot
	const int sz = plt.size();
	IntegerVector plt_booklet_id(sz), plt_item_id(sz), plt_booklet_score(sz), plt_n(sz);
	NumericVector plt_mean(sz);
	
	int i=0;
	for(std::unordered_map<key, val, int3_hash>::iterator iter = plt.begin(); iter != plt.end(); ++iter)
	{
		auto& k = iter->first;
		auto& s = iter->second;
		
		plt_booklet_id[i] = std::get<0>(k);
		plt_booklet_score[i] = std::get<1>(k);
		plt_item_id[i] = std::get<2>(k);
		plt_n[i] = std::get<1>(s);
		plt_mean[i++] = ((double)(std::get<0>(s)))/std::get<1>(s);		
	}
	
	plt_booklet_id.attr("levels") = booklet_id.attr("levels");	
	plt_booklet_id.attr("class") = "factor";
	
	plt_item_id.attr("levels") = item_id.attr("levels");	
	plt_item_id.attr("class") = "factor";
	
	// shrink ssIS to remove non-existing scores
	int indx = 0;
	for(int itm=0; itm<nit; itm++)
	{
		for(int s=0; s<=max_score; s++)
		{
			if(ssI[itm * nscore + s] > 0)
			{
				ssI[indx++] = ssI[itm * nscore + s];
				ssi_item.push_back(itm+1); // 1-indexed
				ssi_item_score.push_back(s);
			}
		}
	}
	
	ssi_item.shrink_to_fit();
	ssi_item_score.shrink_to_fit();
	ssI.resize(indx);
	
	return List::create(
		Named("plt") = DataFrame::create(Named("booklet_id") = plt_booklet_id, Named("booklet_score") = plt_booklet_score, 
										 Named("item_id") = plt_item_id, Named("meanScore") = plt_mean, Named("N") = plt_n),
		Named("ssIS") = DataFrame::create(Named("item_id") = ssi_item, Named("item_score") = ssi_item_score, Named("sufI") = ssI));	

}

// interaction model

// [[Rcpp::export]]
List suf_stats_im(const IntegerVector& booklet_score, const IntegerVector& item_id, const IntegerVector& item_score, const int nit, const int max_score)
{
	const int n = item_id.length();	

	// reserve space for ssIS
	const int nscore = max_score + 1;	
	const int ssIS_len = nit * nscore;
	
	std::vector<int> ss_i(ssIS_len, 0), ss_c(ssIS_len, 0), ssi_item, ssi_item_score;
	ssi_item.reserve(ssIS_len);
	ssi_item_score.reserve(ssIS_len);

	
	// space for plt
	const int max_testscore = nit * max_score;
	const int plt_len = nit * (max_testscore + 1);
	std::vector<int> plt_n(plt_len, 0), plt_sum(plt_len, 0);	
	std::vector<int> plt_item_id, plt_booklet_score;
	std::vector<double> plt_mean;
	plt_item_id.reserve(plt_len);
	plt_booklet_score.reserve(plt_len);
	plt_mean.reserve(plt_len);
		
	// main loop
	for(int i=0; i<n; i++)
	{
		plt_n[booklet_score[i] * nit + item_id[i]-1]++;
		plt_sum[booklet_score[i] * nit + item_id[i]-1] += item_score[i];
		ss_i[(item_id[i]-1) * nscore + item_score[i]]++;
		ss_c[(item_id[i]-1) * nscore + item_score[i]] += booklet_score[i];
	}
	
	// compute plt
	int j = 0;
	for(int sc=0; sc<=max_testscore; sc++)	
	{
		for(int itm = 0; itm<nit; itm++)
		{
			int indx = sc * nit + itm;
			if(plt_n[indx] > 0)
			{
				plt_item_id.push_back(itm+1);
				plt_booklet_score.push_back(sc);
				plt_mean.push_back( ((double)plt_sum[indx]) / plt_n[indx]);
				plt_n[j++] = plt_n[indx];
			}
		}
	}
	plt_item_id.shrink_to_fit();
	plt_mean.shrink_to_fit();
	plt_booklet_score.shrink_to_fit();
	plt_n.resize(plt_item_id.size());

	// shrink ssIS to remove non-existing scores
	int indx = 0;
	for(int itm=0; itm<nit; itm++)
	{
		for(int s=0; s<=max_score; s++)
		{			
			if(ss_i[itm * nscore + s] > 0)
			{
				ss_c[indx] = ss_c[itm * nscore + s] * s;
				ss_i[indx++] = ss_i[itm * nscore + s];
				ssi_item.push_back(itm+1); // 1-indexed
				ssi_item_score.push_back(s);
			}
		}
	}
	
	ssi_item.shrink_to_fit();
	ssi_item_score.shrink_to_fit();
	ss_i.resize(indx);
	ss_c.resize(indx);
	
	return List::create(
		Named("plt") = DataFrame::create(Named("booklet_score") = plt_booklet_score, Named("item_id") = plt_item_id, 
										 Named("meanScore") = plt_mean, Named("N") = plt_n),
		Named("ssIS") = DataFrame::create(Named("item_id") = ssi_item, Named("item_score") = ssi_item_score, 
										  Named("sufI") = ss_i, Named("sufC_") = ss_c));	

}

// scoretab for a single booklet from sumscore vector, a small but noticable speed increase for large vectors

/*
* scrs %>%
*    count(.data$score) %>%
*    right_join(tibble(score=0L:mx),by='score') %>%
*    mutate(n=coalesce(n,0L)) %>%
*    arrange(.data$score) %>%
*    pull(.data$n))
*/


// [[Rcpp::export]]
IntegerVector score_tab_single(const IntegerVector& scores, const int max_score)
{
	const int n = scores.length();
	IntegerVector out(max_score+1); // filled with zero by default
	
	for(int i=0; i<n; i++)
		out[scores[i]]++;
 
	return out;
}


// about 10-30 times faster than R (small bonus for rasch) (to~do: check for severly incomplete designs)
// [[Rcpp::export]]
DataFrame tia_C(const IntegerVector& booklet_id, const IntegerVector& booklet_score, const IntegerVector& item_id, const IntegerVector& item_score, const int nb, const int nit, 
				const IntegerVector& frst_item, const IntegerVector& ds_booklet_id, const IntegerVector& ds_item_id ) 
{
	
	const int n = booklet_id.length();

	// booklet or item things 1 indexed, item_booklet 0 indexed
	std::vector<int> btally(nb+1), bsum(nb+1), bisum(nit*nb), imax(nit+1, 0);
	std::vector<unsigned long long> bsum2(nb+1), bisum2(nit*nb),  bic(nit*nb);
	
	for(int i=0; i<n; i++)
	{
		if(frst_item[booklet_id[i]] == item_id[i])
		{
			bsum[booklet_id[i]] += booklet_score[i];
			bsum2[booklet_id[i]] += booklet_score[i] * booklet_score[i];
			btally[booklet_id[i]]++;
		}
		if(item_score[i] > 0)
		{
			int indx = (booklet_id[i]-1) * nit + item_id[i] - 1;
			int s = item_score[i];
			if(s == 1)
			{
				bisum[indx]++;
				bisum2[indx]++;
				bic[indx] += booklet_score[i];
			}
			else
			{						
				bisum[indx] += s;
				bisum2[indx] += s*s;
				bic[indx] += s * booklet_score[i];					
			}
			if(s > imax[item_id[i]])
				imax[item_id[i]] = s;
		}
	}
	
	const int len_ds = ds_booklet_id.length();
	NumericVector tia_mean(len_ds), tia_sd(len_ds), tia_rit(len_ds), tia_rir(len_ds);
	IntegerVector tia_max(len_ds), tia_n(len_ds);
	
	for(int i=0; i<len_ds; i++)
	{
		tia_max[i] = imax[ds_item_id[i]];
		tia_n[i] = btally[ds_booklet_id[i]];
		
		int indx = (ds_booklet_id[i]-1) * nit + ds_item_id[i] - 1;
		double N = (double)(tia_n[i]);
		double bmean = bsum[ds_booklet_id[i]]/N;
		
		tia_mean[i] = bisum[indx]/N;
		if(tia_n[i] * tia_max[i] != bisum[indx]) // get around float comparison
		{
			tia_sd[i] = std::sqrt(bisum2[indx]/N - tia_mean[i] * tia_mean[i]);
				
			tia_rit[i] = (bic[indx]/N - tia_mean[i]*bmean) /
							(tia_sd[i] * std::sqrt(bsum2[ds_booklet_id[i]]/N - bmean*bmean));
				
			bmean -= tia_mean[i];			
			tia_rir[i] = ((bic[indx] - bisum2[indx])/N - tia_mean[i]*bmean) /
						  (tia_sd[i] * std::sqrt(bisum2[indx]/N + bsum2[ds_booklet_id[i]]/N -2*bic[indx]/N  - bmean*bmean));
		} else
		{
			tia_sd[i] = 0;
			tia_rit[i] = NA_REAL;
			tia_rir[i] = NA_REAL;
		}
	}
	
	
	return DataFrame::create(Named("booklet_id") = ds_booklet_id, Named("item_id") = ds_item_id, Named("meanScore") = tia_mean, 
							 Named("maxScore") = tia_max, Named("sdScore") = tia_sd, Named("rit") = tia_rit, Named("rir") = tia_rir, 
							 Named("n") = tia_n);
}


