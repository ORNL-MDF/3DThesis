//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright 2019 UT-Battelle, LLC
*
* All Rights Reserved
*
* Authors: Benjamin Stump <stumpbc@ornl.gov>, Alex Plotkowski, James Ferguson, Kevin Sisco
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
*	 this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of 3DThesis nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#include <omp.h>
#include <cmath>

#include <iostream>

#include "Run.h"
#include "DataStructs.h"
#include "Point.h"
#include "Calc.h"
#include "Util.h"
#include "Melt.h"
#include "PINT.h"
#include "Out.h"

void Run::Simulate(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {
	Util::EstimateEndTime(sim, segv);
	if (sim.param.mode == 0) { 
		Run::Mode_0(ptv, segv, sim, seg_num); 
	}
	else if (sim.param.mode == 1) { 
		Run::Mode_1(ptv, segv, sim, seg_num); 
	}
	else if (sim.param.mode == 2) {
		if (sim.param.use_PINT) { Run::Mode_2_PINT(ptv, segv, sim, seg_num); }
		else { Run::Mode_2(ptv, segv, sim, seg_num); }
	}
	else if (sim.param.mode == 3) {
		if (sim.param.use_PINT) { Run::Mode_3_PINT(ptv, segv, sim, seg_num); }
		else { Run::Mode_3(ptv, segv, sim, seg_num); }
	}

	// CUSTOM AI CODE
	/*std::vector<Point> ptv_mode0 = ptv;
	Run::Mode_0(ptv_mode0, segv, sim, seg_num);
	Run::Mode_3(ptv, segv, sim, seg_num);*/

	// MATT CODE
	//Run::Mode_2(ptv, segv, sim, seg_num);

	return;
}

void Run::Mode_0(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {

	int itert = seg_num.size(), liq_num = 0;
	double t = sim.util.scanEndTime;

	// CUSTOM AI CODE
	/*t = segv[segv.size() - 2].seg_time;
	itert = floor(t / sim.param.dt);*/

	//Pre-calculate integration loop information
	std::vector<std::vector<int_seg>> isegv_par;
	std::vector<int_seg> isegv;

	while (true) {
		Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, 1);

		int p_tot = 0;
		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(dynamic,1+sim.param.pnum/sim.setting.thnum/128)
		for (int p = 0; p < sim.param.pnum; p++) {
		//for (int p = sim.param.pnum-1; p >=0; p-=10) { // CUSTOM AI CODE
			ptv[p].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0);
			ptv[p].set_output_flag(1);
			#pragma omp atomic
			p_tot++;
			if (!omp_get_thread_num()) {Out::Point_Progress(sim, p_tot);}
		}
		t += sim.param.dt;						//Increment time
		itert++;								//Increment time step counter
		Out::Progress(sim, itert);
		if (true) { break; }					//Check if simulation is finished
	}

	
	return;
}

void Run::Mode_1(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {
	int itert = 0, liq_num = 0;
	double t = segv[0].seg_time;

	//Integration information
	std::vector<std::vector<int_seg>> isegv_par;
	std::vector<int_seg> isegv;

	#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
	for (int p = 0; p < sim.param.pnum; p++) { ptv[p].set_output_flag(1); }

	while (true) {
		Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, sim.setting.thnum);

		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
		for (int p = 0; p < sim.param.pnum; p++) { ptv[p].set_Tlast();}

		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(dynamic,1+sim.param.pnum/sim.setting.thnum/64)
		for (int p = 0; p < sim.param.pnum; p++) {
			if (ptv[p].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq){
				if (ptv[p].get_Tlast() >= sim.mat.T_liq) {ptv[p].CalcGo_2(t, segv, sim);}
			}
			else {
				#pragma omp atomic
				liq_num++;
			}
		}
		
		Out::Progress(sim, itert);
		if (itert && (itert % sim.param.out_freq == 0)) {
			Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); //Output data
		}

		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_num)) { 
			if (itert / sim.param.out_freq) { Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); }
			break; 
		}

		t += sim.param.dt;							//Increment time
		itert++;									//Increment time step counter
		liq_num = 0;
		isegv.clear();
	}
	
	return;
}

void Run::Mode_2(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time
	std::vector<omp_lock_t> lock(sim.param.pnum);
	Util::SetLocks(lock, sim);

	// Vector of Liquid Point numbers
	std::vector<int> liq_pts;

	// Vector of Points to Reset Each Timestep
	std::vector<int> reset_pts;

	//Integration information
	std::vector<std::vector<int_seg>> isegv_par;
	std::vector<int_seg> isegv;

	int itert = 0, liq_num = 0;
	double t = segv[0].seg_time;

	while (true) {
		Out::Progress(sim, itert);

		//Calculate Integration information
		Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, sim.setting.thnum);
		
		//Set T_calc_flag at all points to indicate that they have not yet been calculated
		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
		for (int r = 0; r < reset_pts.size(); r++) {
			ptv[reset_pts[r]].set_Tlast();
			ptv[reset_pts[r]].init_T_calc_flag();
			ptv[reset_pts[r]].set_s_flag(0);
		}
		reset_pts.clear();
		
		//See which points from previous time step have solidified and do appropriate calculations
		std::vector<int> last_liq_pts = liq_pts;
		reset_pts = liq_pts;
		liq_pts.clear();

		//Check liquid points to see if they have solidified
		#pragma omp parallel num_threads(sim.setting.thnum)
		{
			std::vector<int> th_liq_pts;
			std::vector<int> th_reset_pts;
			#pragma omp for schedule(dynamic,1+last_liq_pts.size()/sim.setting.thnum/64)
			for (int it = 0; it < last_liq_pts.size(); it++) {
				if (ptv[last_liq_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) { ptv[last_liq_pts[it]].CalcGo_2(t, segv, sim); }
				else {th_liq_pts.push_back(last_liq_pts[it]);}
			}
			#pragma omp critical
			liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
		}

		//Start search from points that are known to be liquid
		std::vector<int> test_pts = liq_pts;

		// Trace beam path between time steps and add relevant points to the test vector if:
		std::vector<int> bm_tr_pts;
		if (itert && t < sim.util.scanEndTime + sim.param.dt) {Melt::beam_trace(bm_tr_pts, ptv, segv, sim, seg_num, itert, itert - 1);}
		for (int it = 0; it < bm_tr_pts.size(); it++) {
			reset_pts.push_back(bm_tr_pts[it]);
			if (ptv[bm_tr_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
				test_pts.push_back(bm_tr_pts[it]);
				liq_pts.push_back(bm_tr_pts[it]);
			}
		}

		//Iterative loop expanding from previously identified points to find melt pool boundary
		while (true) {		
			if (!test_pts.size()) { break; }
			Melt::neighbor_check(test_pts, liq_pts, reset_pts, ptv, lock, t, isegv, sim, 0);
		}

		//Out::Progress_Iter(sim, itert); //Output progress
		if (itert && (itert % sim.param.out_freq == 0)) {
			Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); //Output data
		}
		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_pts.size())) {
			if (itert / sim.param.out_freq) { Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); }
			break; 
		}	
		t += sim.param.dt;		//Increment time
		itert++;				//Increment time step counter
		isegv.clear();
	}
	
	return;
}

void Run::Mode_2_PINT(std::vector<Point>& ptv_master, std::vector<path_seg>& segv, Simdat& sim_master, std::vector<int>& seg_num) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time (last is for thread specific locking)
	std::vector<omp_lock_t> lock_master(sim_master.param.pnum);
	Util::SetLocks(lock_master, sim_master);
	int num_done = 0;
	int itert_tot = 0;

	std::vector <int> th_run(sim_master.setting.thnum, 1);
	int last_thread = sim_master.setting.thnum - 1;

	#pragma omp parallel for num_threads(sim_master.setting.thnum) schedule(dynamic) shared(ptv_master,lock_master,num_done,th_run,last_thread)
	for (int thread = 0; thread < sim_master.setting.thnum; thread++) {
		Simdat sim = sim_master;
		//Sets locks so only 1 thread can access a point at the same time (last is for thread specific locking)
		std::vector<omp_lock_t> lock(sim.param.pnum+1);
		Util::SetLocks(lock, sim);
		//Vector of ints whose index corresponds to the ptv_master and it's value corresponds to the index of the same point in ptv
		std::vector<int> god_pts(sim_master.param.pnum, -1);
		std::deque<Point> ptv;
		// Vector of Liquid Point numbers
		std::vector<int> liq_pts;
		// Vector of Points to Reset Each Timestep
		std::vector<int> reset_pts;
		//Integration information
		std::vector<std::vector<int_seg>> isegv_par;
		std::vector<int_seg> isegv;
		///////////////////////////////////////
		int itert = 0, liq_num = 0, itert_thread=0, itert_end=0;
		double t = segv[0].seg_time;
		double speed_pow = 0.75;
		PINT::calcSpeedPow(speed_pow, sim);
		PINT::calc_iterts(itert, itert_end, sim, thread, speed_pow);
		t = itert * sim.param.dt;
		///////////////////////////////////////
		int num_free = 0;
		while (true) {
			if (thread == last_thread) { num_free = num_done;}
			//Finds the time segment to start searching
			Util::GetStartSeg(sim, seg_num, itert);

			//Pre-calculate integration loop information
			isegv.clear();
			Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, 1 + num_free);

			//Set T_calc_flag at all points to indicate that they have not yet been calculated
			#pragma omp parallel for num_threads(1+num_free) schedule(static) if(num_free)
			for (int r = 0; r < reset_pts.size(); r++) {
				ptv[god_pts[reset_pts[r]]].init_T_calc_flag();
				ptv[god_pts[reset_pts[r]]].set_s_flag(0);
			}
			reset_pts.clear();

			//See which points from previous time step have solidified and do appropriate calculations
			std::vector<int> last_liq_pts = liq_pts;
			reset_pts = last_liq_pts;
			liq_pts.clear();

			#pragma omp parallel num_threads(1+num_free) if(num_free)
			{
				std::vector<int> th_liq_pts;
				#pragma omp for schedule(dynamic, 1+last_liq_pts.size() / (1+num_free) / 64)
				for (int it = 0; it < last_liq_pts.size(); it++) {
					if (ptv[god_pts[last_liq_pts[it]]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) {
						ptv[god_pts[last_liq_pts[it]]].CalcGo_2(t, segv, sim);
					}
					else { th_liq_pts.push_back(last_liq_pts[it]); }
				}
				omp_set_lock(&(lock[sim.param.pnum]));
				liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
				omp_unset_lock(&(lock[sim.param.pnum]));
			}

			//Start search from points that are known to be liquid
			std::vector<int> test_pts = liq_pts;

			std::vector<int> bm_tr_pts;
			// Trace beam path between time steps and add relevant points to the test vector if:
			if (itert && !itert_thread){ PINT::beam_trace(bm_tr_pts, god_pts, ptv, lock, segv, sim, seg_num, itert, 0); }
			else if (itert && t < sim.util.scanEndTime + sim.param.dt) {PINT::beam_trace(bm_tr_pts, god_pts, ptv, lock, segv, sim, seg_num, itert, itert - 1);}
			
			for (int it = 0; it < bm_tr_pts.size(); it++) {
				reset_pts.push_back(bm_tr_pts[it]);
				if (ptv[god_pts[bm_tr_pts[it]]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
					test_pts.push_back(bm_tr_pts[it]);

				}
			}

			//Iterative loop expanding from previously identified points to find melt pool boundary
			while (true) {
				if (!test_pts.size()) { break; }
				PINT::neighbor_check(test_pts, liq_pts, reset_pts, god_pts, ptv, lock, t, isegv, sim, num_free, 0);
			}

			liq_num = liq_pts.size();
			if (itert && (itert % sim.param.out_freq == 0)) { Out::Write_csv_PINT(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); } //Output data
			//Check if simulation is finished
			if (thread + 1 == sim.setting.thnum && Util::sim_finish(t, sim, liq_num)) { break; }
			else if (thread + 1 < sim.setting.thnum && itert == itert_end) { break; }

			t += sim.param.dt;	//Increment time
			itert++;			//Increment time step counter
			itert_thread++;
		}
		//Out::Write_csv(ptv, sim, Util::ZeroPadNumber(thread), sim.setting.out_mode);
		PINT::GodToPtv(ptv_master, god_pts, ptv, sim, lock_master);
		#pragma omp critical
		{	
			th_run[thread] = 0;
			for (int i = th_run.size() - 1; i >= 0; i--) {
				if (th_run[i]) { 
					last_thread = i; 
					break; 
				}
			}
			itert_tot += itert_thread;
			num_done++;
			Out::Progress(sim, itert_tot);
		}	
	}
	return;
}

void Run::Mode_3(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time
	std::vector<omp_lock_t> lock(sim.param.pnum);
	Util::SetLocks(lock, sim);

	// Vector of Liquid Point numbers
	std::vector<int> liq_pts;
	std::vector<int> depths(sim.param.xnum*sim.param.ynum, 0);

	// Vector of Points to Reset Each Timestep
	std::vector<int> reset_pts;

	//Integration information
	std::vector<std::vector<int_seg>> isegv_par;
	std::vector<int_seg> isegv;
	std::vector<int_seg> isegv_last;

	int itert = 0, liq_num = 0;
	double t = segv[0].seg_time;
	
	//CUSTOM AI CODE
	/*t = segv[segv.size() - 2].seg_time;
	itert = floor(t / sim.param.dt);*/

	while (true) {
		Out::Progress(sim, itert);

		//Calculate Integration information
		Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, sim.setting.thnum);

		//Set T_calc_flag at all points to indicate that they have not yet been calculated
		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
		for (int r = 0; r < reset_pts.size(); r++) {
			ptv[reset_pts[r]].set_Tlast();
			ptv[reset_pts[r]].init_T_calc_flag();
			ptv[reset_pts[r]].set_s_flag(0);
		}
		reset_pts.clear();

		//See which points from previous time step have solidified and do appropriate calculations
		std::vector<int> last_liq_pts = liq_pts;
		reset_pts = liq_pts;
		liq_pts.clear();

		//Check liquid points to see if they have solidified
		#pragma omp parallel num_threads(sim.setting.thnum)
		{
			std::vector<int> th_liq_pts;
			std::vector<int> th_reset_pts;
			#pragma omp for schedule(dynamic,1+last_liq_pts.size()/sim.setting.thnum/64)
			for (int it = 0; it < last_liq_pts.size(); it++) {
				if (ptv[last_liq_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) { 
					ptv[last_liq_pts[it]].CalcGo_2(t, segv, sim);
					///////SOLIDIFY REST OF COLUMN/////////	
					int i = ptv[last_liq_pts[it]].get_i();
					int j = ptv[last_liq_pts[it]].get_j();
					int dnum = i * sim.param.ynum + j;
					int depth = depths[dnum];
					for (int d = depth; d > 0; d--) {
						int pnum = Util::ijk_to_p(i, j, sim.param.znum - 1 - d, sim);
						ptv[pnum].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0);
						if (d != depth) {ptv[pnum].Temp_Calc_Pre_Path(t-sim.param.dt, isegv_last, sim, -1, 0);}
						ptv[pnum].CalcGo_2(t, segv, sim);
						th_reset_pts.push_back(pnum);
					}
					depths[dnum] = 0;
				}
				else {th_liq_pts.push_back(last_liq_pts[it]);}
			}
			#pragma omp critical
			{
				liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
				reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
			}
		}

		//Start search from points that are known to be liquid
		std::vector<int> test_pts = last_liq_pts;

		// Trace beam path between time steps and add relevant points to the test vector
		std::vector<int> bm_tr_pts;
		if (itert && t < sim.util.scanEndTime + sim.param.dt) { Melt::beam_trace(bm_tr_pts, ptv, segv, sim, seg_num, itert, itert - 1); }
		for (int it = 0; it < bm_tr_pts.size(); it++) {
			reset_pts.push_back(bm_tr_pts[it]);
			if (ptv[bm_tr_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
				test_pts.push_back(bm_tr_pts[it]);
				liq_pts.push_back(bm_tr_pts[it]);
			}
		}

		//Iterative loop expanding from previously identified points to find melt pool boundary
		while (true) {
			if (!test_pts.size()) { break; }
			Melt::neighbor_check(test_pts, liq_pts, reset_pts, ptv, lock, t, isegv, sim, 1);
		}

		//Check the depths of the meltpool, solidify if needed
		Melt::calc_depth(depths, liq_pts, reset_pts, ptv, t, isegv, isegv_last, segv, sim);

		//Out::Progress_Iter(sim, itert); //Output progress
		if (itert && (itert % sim.param.out_freq == 0)) {
			Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); //Output data
		}
		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_pts.size())) {
			if (itert / sim.param.out_freq) {Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode);}
			break;
		}
		t += sim.param.dt;			//Increment time
		itert++;					//Increment time step counter
		isegv_last = isegv;
		isegv.clear();
	}
	return;
}

void Run::Mode_3_PINT(std::vector<Point>& ptv_master, std::vector<path_seg>& segv, Simdat& sim_master, std::vector<int>& seg_num) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time (last is for thread specific locking)
	std::vector<omp_lock_t> lock_master(sim_master.param.pnum);
	Util::SetLocks(lock_master, sim_master);
	int num_done = 0;
	int itert_tot = 0;

	std::vector <int> th_run(sim_master.setting.thnum, 1);
	int last_thread = sim_master.setting.thnum - 1;

	#pragma omp parallel for num_threads(sim_master.setting.thnum) schedule(dynamic) shared(ptv_master,lock_master,num_done,th_run,last_thread)
	for (int thread = 0; thread < sim_master.setting.thnum; thread++) {
		Simdat sim = sim_master;
		//Sets locks so only 1 thread can access a point at the same time (last is for thread specific locking)
		std::vector<omp_lock_t> lock(sim.param.pnum+1);
		Util::SetLocks(lock, sim);
		//Vector of ints whose index corresponds to the ptv_master and it's value corresponds to the index of the same point in ptv
		std::vector<int> god_pts(sim_master.param.pnum, -1);
		std::deque<Point> ptv;
		// Vector of Liquid Point numbers
		std::vector<int> liq_pts;
		std::vector<int> depths(sim.param.xnum*sim.param.ynum, 0);
		// Vector of Points to Reset Each Timestep
		std::vector<int> reset_pts;
		//Integration information
		std::vector<std::vector<int_seg>> isegv_par;
		std::vector<int_seg> isegv;
		std::vector<int_seg> isegv_last;
		///////////////////////////////////////
		int itert = 0, liq_num = 0, itert_thread = 0, itert_end = 0;
		double t = segv[0].seg_time;
		double speed_pow = 0.75;
		PINT::calcSpeedPow(speed_pow, sim);
		PINT::calc_iterts(itert, itert_end, sim, thread, speed_pow);
		t = itert * sim.param.dt;
		///////////////////////////////////////
		int num_free = 0;

		while (true) {
			if (thread == last_thread) { num_free = num_done; }

			//Calculate Integration information
			isegv.clear();
			Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, 1 + num_free);

			//Set T_calc_flag at all points to indicate that they have not yet been calculated
			#pragma omp parallel for num_threads(1+num_free) schedule(static) if(num_free)
			for (int r = 0; r < reset_pts.size(); r++) {
				ptv[god_pts[reset_pts[r]]].set_Tlast();
				ptv[god_pts[reset_pts[r]]].init_T_calc_flag();
				ptv[god_pts[reset_pts[r]]].set_s_flag(0);
			}
			reset_pts.clear();

			//See which points from previous time step have solidified and do appropriate calculations
			std::vector<int> last_liq_pts = liq_pts;
			reset_pts = liq_pts;
			liq_pts.clear();

			//Check liquid points to see if they have solidified
			#pragma omp parallel num_threads(1+num_free) if(num_free)
			{
				std::vector<int> th_liq_pts;
				std::vector<int> th_reset_pts;
				#pragma omp for schedule(dynamic,1+ last_liq_pts.size() / (1+num_free) / 64)
				for (int it = 0; it < last_liq_pts.size(); it++) {
					if (ptv[god_pts[last_liq_pts[it]]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) {
						ptv[god_pts[last_liq_pts[it]]].CalcGo_2(t, segv, sim);
						///////SOLIDIFY REST OF COLUMN/////////	
						int i = ptv[god_pts[last_liq_pts[it]]].get_i();
						int j = ptv[god_pts[last_liq_pts[it]]].get_j();
						int dnum = i * sim.param.ynum + j;
						int depth = depths[dnum];
						for (int d = depth; d > 0; d--) {
							int pnum = Util::ijk_to_p(i, j, sim.param.znum - 1 - d, sim);
							PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
							ptv[god_pts[pnum]].set_output_flag(1);
							ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0);
							if (d != depth) { ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t - sim.param.dt, isegv_last, sim, -1, 0); }
							ptv[god_pts[pnum]].CalcGo_2(t, segv, sim);
							th_reset_pts.push_back(pnum);
						}
						depths[dnum] = 0;
					}
					else { th_liq_pts.push_back(last_liq_pts[it]); }
				}
				omp_set_lock(&(lock[sim.param.pnum]));
				{
					liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
					reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
				}
				omp_unset_lock(&(lock[sim.param.pnum]));
			}

			//Start search from points that are known to be liquid
			std::vector<int> test_pts = last_liq_pts;

			// Trace beam path between time steps and add relevant points to the test vector
			std::vector<int> bm_tr_pts;
			if (itert && !itert_thread) { PINT::beam_trace(bm_tr_pts, god_pts, ptv,lock, segv, sim, seg_num, itert, 0); }
			else if (itert && t < sim.util.scanEndTime + sim.param.dt) { PINT::beam_trace(bm_tr_pts, god_pts, ptv,lock, segv, sim, seg_num, itert, itert - 1); }
			//if (itert && t < sim.util.scanEndTime + sim.param.dt) { Melt::beam_trace(bm_tr_pts, ptv, segv, sim, seg_num, itert, itert - 1); }
			for (int it = 0; it < bm_tr_pts.size(); it++) {
				reset_pts.push_back(bm_tr_pts[it]);
				if (ptv[god_pts[bm_tr_pts[it]]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
					test_pts.push_back(bm_tr_pts[it]);
					liq_pts.push_back(bm_tr_pts[it]);
				}
			}

			//Iterative loop expanding from previously identified points to find melt pool boundary
			while (true) {
				if (!test_pts.size()) { break; }
				PINT::neighbor_check(test_pts, liq_pts, reset_pts, god_pts, ptv, lock, t, isegv, sim, num_free, 1);
			}

			//Check the depths of the meltpool, solidify if needed
			PINT::calc_depth(depths, liq_pts, reset_pts, god_pts, ptv,lock, t, isegv, isegv_last, segv, sim, num_free);

			if (itert && (itert % sim.param.out_freq == 0)) { Out::Write_csv_PINT(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); } //Output data
			//Check if simulation is finished
			liq_num = liq_pts.size();
			if (thread + 1 == sim.setting.thnum && Util::sim_finish(t, sim, liq_num)) { break; }
			else if (thread + 1 < sim.setting.thnum && itert == itert_end) { break; }

			t += sim.param.dt;			//Increment time
			itert++;					//Increment time step counter
			itert_thread++;
			isegv_last = isegv;
			isegv.clear();
		}
		//Out::Write_csv(ptv, sim, Util::ZeroPadNumber(thread), sim.setting.out_mode);
		PINT::GodToPtv(ptv_master, god_pts, ptv, sim, lock_master);
		#pragma omp critical
		{
			th_run[thread] = 0;
			for (int i = th_run.size() - 1; i >= 0; i--) {
				if (th_run[i]) {
					last_thread = i;
					break;
				}
			}
			itert_tot += itert_thread;
			num_done++;
			Out::Progress(sim, itert_tot);
		}
	}
	return;
}


//MATT CUSTOM
//void Run::Mode_2(std::vector<Point>& ptv, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num) {
//	sim.mat.T_liq = 1616;
//	sim.mat.T_sol = 1446;
//	
//	omp_set_nested(1);
//
//	//Sets locks so only 1 thread can access a master point at the same time
//	std::vector<omp_lock_t> lock(sim.param.pnum);
//	Util::SetLocks(lock, sim);
//
//	// Vector of Liquid Point numbers
//	std::vector<int> liq_pts;
//
//	// Vector of Points to Reset Each Timestep
//	std::vector<int> reset_pts;
//
//	//Integration information
//	std::vector<std::vector<int_seg>> isegv_par;
//	std::vector<int_seg> isegv;
//
//	int itert = 0, liq_num = 0;
//	double t = segv[0].seg_time;
//
//	while (true) {
//		Out::Progress(sim, itert);
//
//		//Calculate Integration information
//		Calc::Integrate(isegv, isegv_par, segv, sim, seg_num, itert, t, 0, sim.setting.thnum);
//
//		//Set T_calc_flag at all points to indicate that they have not yet been calculated
//		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
//		for (int r = 0; r < reset_pts.size(); r++) {
//			ptv[reset_pts[r]].set_Tlast();
//			ptv[reset_pts[r]].init_T_calc_flag();
//			ptv[reset_pts[r]].set_s_flag(0);
//		}
//		reset_pts.clear();
//
//		//See which points from previous time step have solidified and do appropriate calculations
//		std::vector<int> last_liq_pts = liq_pts;
//		reset_pts = liq_pts;
//		liq_pts.clear();
//
//		//Check liquid points to see if they have solidified
//		#pragma omp parallel num_threads(sim.setting.thnum)
//		{
//			std::vector<int> th_liq_pts;
//			std::vector<int> th_reset_pts;
//			#pragma omp for schedule(dynamic,last_liq_pts.size()/sim.setting.thnum/64)
//			for (int it = 0; it < last_liq_pts.size(); it++) {
//				if (ptv[last_liq_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq && ptv[last_liq_pts[it]].get_Tlast()>=sim.mat.T_liq ) { 
//					ptv[last_liq_pts[it]].set_t_last_liq(ptv[last_liq_pts[it]].Calc_Time(t, segv, sim, sim.mat.T_liq));
//					ptv[last_liq_pts[it]].set_t_last_sol(0.0);
//				}
//				if (ptv[last_liq_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_sol) {
//					if (ptv[last_liq_pts[it]].get_t_last_sol() == 0.0) {
//						ptv[last_liq_pts[it]].set_t_last_sol(ptv[last_liq_pts[it]].Calc_Time(t, segv, sim, sim.mat.T_sol));
//					}
//				}
//				else { th_liq_pts.push_back(last_liq_pts[it]); }
//			}
//			#pragma omp critical
//			liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
//		}
//
//		//Start search from points that are known to be liquid
//		std::vector<int> test_pts = liq_pts;
//
//		// Trace beam path between time steps and add relevant points to the test vector if:
//		std::vector<int> bm_tr_pts;
//		if (itert && t < sim.util.scanEndTime + sim.param.dt) { Melt::beam_trace(bm_tr_pts, ptv, segv, sim, seg_num, itert, itert - 1); }
//		for (int it = 0; it < bm_tr_pts.size(); it++) {
//			reset_pts.push_back(bm_tr_pts[it]);
//			double T = ptv[bm_tr_pts[it]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0);
//			double z = ptv[bm_tr_pts[it]].get_z();
//			if (T >= sim.mat.T_sol && (T <= sim.mat.T_liq + 100.0 || z==0.0) ) {
//				test_pts.push_back(bm_tr_pts[it]);
//				liq_pts.push_back(bm_tr_pts[it]);
//			}
//		}
//
//		//Iterative loop expanding from previously identified points to find melt pool boundary
//		while (true) {
//			if (!test_pts.size()) { break; }
//			Melt::neighbor_check(test_pts, liq_pts, reset_pts, ptv, lock, t, isegv, sim, 0);
//		}
//
//		//Out::Progress_Iter(sim, itert); //Output progress
//		if (itert && (itert % sim.param.out_freq == 0)) {
//			Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); //Output data
//		}
//		//Check if simulation is finished
//		if (Util::sim_finish(t, sim, liq_pts.size())) {
//			if (itert / sim.param.out_freq) { Out::Write_csv(ptv, sim, Util::ZeroPadNumber(itert), sim.setting.out_mode); }
//			break;
//		}
//		t += sim.param.dt;		//Increment time
//		itert++;				//Increment time step counter
//		isegv.clear();
//	}
//
//	return;
//}