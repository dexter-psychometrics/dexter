#include <atomic>
#include <string>
#include <RcppArmadillo.h>
#include "progress.h"

// see https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

progress::progress(const int nsteps_, const arma::ivec& settings)
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
	
void progress::draw_perc()
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
	
void progress::tick(const int nticks)
{
	sub_step = sub_step + nticks;
	double sub_prog = std::min(sub_step/((double)sub_nsteps),1.0);
	step = step_before_sub + (int)(1e-6 +  sub_prog * sub_length);
	draw_perc();
}

progress_prl::progress_prl(const int nsteps_, const arma::ivec& settings)
		: progress(nsteps_ , settings)
{
	atm_tick = 0;
	atm_interrupted = false;
}
	
bool progress_prl::interrupted()
{
	return atm_interrupted.load();
}
	
// NEVER call checkInterrupt outside of the main thread
void progress_prl::checkInterrupt() 
{
	if(R_ToplevelExec(chkIntFn, NULL) == FALSE) // the if is necessary
	{
		atm_interrupted = true; 
	}
}
		
void progress_prl::tick(const bool main_thread, const int nticks)
{
	atm_tick += nticks;
	if(main_thread)
	{
		checkInterrupt();
		if(!interrupted()) // attempting to print anything when interrupted causes a crash
		{
			progress::tick(atm_tick.load());		
			atm_tick = 0;
		}
	}		
}	

