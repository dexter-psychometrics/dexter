#ifndef DEXTER_PROG_
#define DEXTER_PROG_

#include <RcppArmadillo.h>

struct progress
{
	int step, step_before_sub, sub_nsteps, sub_length, nsteps, w, sub_step,p,l;
	std::string fmt;
	
	progress(const int nsteps_, const arma::ivec& settings)
	{
		step_before_sub = settings(0);
		sub_nsteps = nsteps_;
		sub_length = settings(2);
		nsteps = settings(3);
		w = settings(4);
		sub_step=0;
		p=0;
		l=0;
		step=0;
		if(w>0)
			fmt = std::string("\r|%-") + std::to_string(w) + std::string("s| %3i%%");	
	}
	void draw_perc()
	{
		if(w>0)
		{
			int old = l+p;
			step = std::min(step, nsteps);
			double part = ((double)step)/nsteps;
			p = (int)(100*part);
			l = (int)(w*part);
			if(old != l+p)
			{
				Rprintf(fmt.c_str(), std::string(l,'=').c_str(), p);
			}
		}		
	}
	void tick(const int nticks = 1)
	{
		sub_step = sub_step + nticks;
		double sub_prog = std::min(sub_step/((double)sub_nsteps),1.0);
		step = step_before_sub + (int)(1e-6 +  sub_prog * sub_length);
		draw_perc();
	}
};	
	
#endif