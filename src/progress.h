#ifndef DEXTER_PROG_
#define DEXTER_PROG_
#include <atomic>
#include <string>
#include <RcppArmadillo.h>


class progress
{
private:
	int step, step_before_sub, sub_nsteps, sub_length, nsteps, w, sub_step,p,l;
	std::string fmt;
public:	
	progress(const int nsteps_, const arma::ivec& settings);	
	void draw_perc();	
	void tick(const int nticks = 1);

};	

class progress_prl : progress
{
private:
	std::atomic<int> atm_tick;
	std::atomic<bool> atm_interrupted;
public:
	
	progress_prl(const int nsteps_, const arma::ivec& settings);
	
	bool interrupted();
	
	// NEVER call checkInterrupt outside of the main thread
	void checkInterrupt(); 
	// also calls checkinterrupt if main_thread==true	
	void tick(const bool main_thread, const int nticks = 1);
};
	
#endif