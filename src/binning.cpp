
#include <stack>
#include <RcppArmadillo.h>
#include "shared.h"

using namespace arma;

// [[Rcpp::export]]
arma::ivec weighted_binning(const arma::ivec& weights, const int nbins)
{
	const int nw = weights.n_elem;
	const double opt_size = accu(weights)/((double)nbins);
	const ivec weights_remainder = reverse(cumsum(reverse(weights))); 
	
	double cost, current_cost, opt_split2;
	int current_break, w;
	
	// set up hash tables
	struct key {
		int start; 
		int bins;
		bool operator==(const key& other) const
		{
			return (start==other.start && bins == other.bins);
		}
	};
	
	struct key_hash
	{
	  std::size_t operator()(const key& k) const
	  {
		return std::hash<int>{}(k.start*20 + k.bins);
	  }	
	};
	
	std::stack<key> stk;
	key p;
		
	std::unordered_map<key, int, key_hash> parent;
	std::unordered_map<key, bool, key_hash> branched;
    std::unordered_map<key, double, key_hash> memo;
	
	// heuristically determine appropriate min and max size
	w=0;
	double dev = 0;	
	int wbin=0;
	for(int i=0,bin=1; i<nw;i++)
	{
		w += weights[i];
		wbin += weights[i];
		if(w>=opt_size*bin && bin<nbins)
		{
			dev = std::max(dev,std::abs(wbin-opt_size));
			wbin=0;
			bin++;
		}
	}
	dev = std::max(dev,std::abs(wbin-opt_size));
	
	const int min_size = std::floor(opt_size-dev), max_size = std::ceil(opt_size+dev);

	
	stk.push(key{0,nbins});
	
	
	while(!stk.empty())
	{
		p = stk.top();		

		// edge case
		if(nw - p.start <= p.bins)
		{
			memo[p] = std::numeric_limits<double>::infinity();
			stk.pop();
			continue;
		}		
		
		if(p.bins == 2) // one break left, split in two equals
		{
			opt_split2 = weights_remainder[p.start]/2.0;
			w = 0;  
			for(int i=p.start; i<nw; i++)
			{					
				if(w + weights[i] > opt_split2)
				{
					if(std::abs(w-opt_split2) < std::abs(w+weights[i]-opt_split2))
					{
						i -= 1; // previous is best
					} else
					{
						w += weights[i];
					}

					parent[p] = i;  
					memo[p] = SQR(w-opt_size) + SQR(weights_remainder[p.start]-w-opt_size);
					break;
				}
				w += weights[i];     
			}
			stk.pop();	
			continue;
		}
	
		if(branched[p]) 
		{
			// compute
			current_break = 0; w=0;
			current_cost = std::numeric_limits<double>::infinity();
			
			for(int i=p.start;i<nw; i++)
			{			  
				w += weights[i];
				if( w >= min_size)
				{
					if(w > max_size) break;
					
					cost = SQR(w-opt_size) + memo[key{i+1, p.bins-1}];
					if(cost < current_cost)
					{
						current_break = i;
						current_cost = cost;
					}
				}
			}
			parent[p]  = current_break;    
			memo[p] = current_cost;
				
			stk.pop();	
		}
		else
		{
			// branch
			w=0;
			for(int i=p.start;i<nw; i++)
			{
				w += weights[i];
				if( w >= min_size)
				{
					if(w > max_size) break;
					if (branched.find(key{i+1,p.bins-1}) == branched.end())	
					{
						stk.push(key{i+1,p.bins-1});
						branched[key{i+1,p.bins-1}] = false;
					}
						
				}
			}
			// we can do this one next time we encounter it
			branched[p] = true;
		}			
	
	}
	
	// create output vector from parent

	ivec res(nw);
	res.fill(nbins);
     
    int a = 0, b = parent[key{0,nbins}];
    for(int bin=1; bin<=nbins; bin++)
    {
		for(int i=a; i<=b; i++)
		{
			res[i] = bin;
		}
		a = b+1;
		b = parent[key{a,nbins-bin}];		
    }

    return res;
}
  


