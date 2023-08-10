//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright 2019 UT-Battelle, LLC
*
* All Rights Reserved
*
* Authors: Benjamin Stump <stumpbc@ornl.gov> and Alex Plotkowski
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
#include <fstream>

#include "Run.h"
#include "DataStructs.h"
#include "Calc.h"
#include "Util.h"
#include "Melt.h"
#include "Out.h"

void Run::Simulate(Grid& grid, const Simdat& sim) {
	
	if (sim.param.mode == "Snapshots") { Run::Snapshots(grid, sim); }
	else if (sim.param.mode == "Solidification") { Run::Solidify(grid, sim); }
	else { std::cout << "ERROR: Unrecognized mode: " << sim.param.mode << "\n"; }

	// Mode - Solidification
	//// Domain: Domain.txt, Points.txt
	//// Tracking: None, Volume, Perimeter
	//// TimeStep 
	//// OutputFrequency
	//// SecondarySolidification

	// Mode - Temperture History 
	//// Domain: Domain.txt, Points.txt
	//// Timestep: 

	// Mode - Snapshots (only T)
	//// Domain: Domain.txt, Points.txt
	//// Times: (time or scan%)

	// Be able to read unlimited beams

	return;
}

void Run::Snapshots(Grid& grid, const Simdat& sim) {

	if (sim.param.tracking == "None" || sim.domain.customPoints) {
		Run::Snapshots_NoTracking(grid, sim);
	}
	else if (sim.param.tracking == "Volume") {
		Run::Snapshots_Volume(grid, sim);
	}
	else if (sim.param.tracking == "Surface") {
		Run:Snapshots_Volume(grid, sim);
	}
	else {
		std::cout << "ERROR: Unrecognized snapshots tracking: " << sim.param.tracking << "\n";
	}

	return;
}

void Run::Snapshots_NoTracking(Grid& grid, const Simdat& sim) {

	//Pre-calculate integration loop information
	vector<Nodes> isegv_par;
	Nodes nodes;

	// Make sure every point outputs temperature
	#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
	for (int p = 0; p < sim.domain.pnum; p++) {
		grid.set_output_flag(true, p);
	}

	// For all times, set approximate iteration and integrate
	for (int i = 0; i < sim.param.SnapshotTimes.size(); i++) {
		
		// Set time and iteration
		double t = sim.param.SnapshotTimes[i];
		int itert = i;

		// Output time
		std::cout << "Time: " << t << "\n";
		
		// Get quadrature nodes
		Calc::Integrate_Serial(nodes, sim, t, false);
		
		// For all points, calculate temperature
		int p_tot = 0;
		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(dynamic,1+sim.domain.pnum/sim.settings.thnum/128)
		for (int p = 0; p < sim.domain.pnum; p++) {
			grid.Calc_T(t, nodes, sim, true, p);
			#pragma omp atomic
			p_tot++;
			if (!omp_get_thread_num()) { Out::Point_Progress(sim, p_tot); }
		}
		std::cout << "\n";

		// Output 
		grid.Output(sim, "Snapshot." + Util::ZeroPadNumber(itert, 2));
		
		// Clear quadrature nodes
		Util::ClearNodes(nodes);
	}

	return;
}

void Run::Snapshots_Volume(Grid& grid, const Simdat& sim) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time
	vector<omp_lock_t> lock(sim.domain.pnum);
	Util::SetLocks(lock, sim);

	//Pre-calculate integration loop information
	vector<Nodes> isegv_par;
	Nodes nodes;
	
	// Vector of Liquid Point numbers
	vector<int> liq_pts;

	// Vector of Points to Reset Each Timestep
	vector<int> reset_pts;

	// For all times, reset points, set approximate iteration, and get quadrature nodes
	for (int i = 0; i < sim.param.SnapshotTimes.size(); i++) {
		
		// Reset all points to allow for temperature calculation
		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
		for (int r = 0; r < reset_pts.size(); r++) {
			const int p = reset_pts[r];
			grid.set_output_flag(false, p);
			grid.set_T_calc_flag(false, p);
		}
		reset_pts.clear();
		
		// Set time and "iteration"
		const double t = sim.param.SnapshotTimes[i];
		const int itert = i;

		// Output time
		std::cout << "Time: " << t << "\n";

		// Get quadrature nodes
		Calc::Integrate_Serial(nodes, sim, t, false);

		// Trace path from now until start of scan to seed initial points
		vector<int> bm_tr_pts;
		if (itert) { Melt::beam_trace(bm_tr_pts, grid, sim, 0, t); }
		
		// Calculate temperature of seed points and see what is liquid
		vector<int> test_pts = liq_pts;
		for (int it = 0; it < bm_tr_pts.size(); it++) {
			const int p = bm_tr_pts[it];
			reset_pts.push_back(p);
			if (grid.Calc_T(t, nodes, sim, true, p) >= sim.material.T_liq) {
				test_pts.push_back(p);
				liq_pts.push_back(p);
			}
		}

		//Iterative loop expanding from previously identified points to find melt pool boundary
		while (true) {
			if (!test_pts.size()) { break; }
			Melt::neighbor_check(test_pts, liq_pts, reset_pts, grid, lock, nodes, sim, t, false);
		}
		
		// Output results
		grid.Output(sim, "Snapshot." + Util::ZeroPadNumber(i, 2));

		// Clear integration segments
		Util::ClearNodes(nodes);
	}

}

void Run::Solidify(Grid& grid, const Simdat& sim){
	if (sim.param.tracking == "None" || sim.domain.customPoints) {
		Run::Solidify_NoTracking(grid, sim);
	}
	else if (sim.param.tracking == "Volume") {
		Run::Solidify_Volume(grid, sim);
	}
	else if (sim.param.tracking == "Surface") {
		Run::Solidify_Surface(grid, sim);
	}
	else {
		std::cout << "ERROR: Unrecognized solidfication tracking: " << sim.param.tracking << "\n";
	}
	return;
}

void Run::Solidify_NoTracking(Grid& grid, const Simdat& sim) {
	int itert = 0;
	int liq_num = 0;
	double t = 0.0;

	//Integration information
	vector<Nodes> isegv_par;
	Nodes nodes;

	#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
	for (int p = 0; p < sim.domain.pnum; p++) { 
		grid.set_T(sim.material.T_init, p);
	}

	while (true) {
		Out::Progress(sim, itert);
		
		Calc::Integrate_Parallel(nodes, sim, t, false);

		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
		for (int p = 0; p < sim.domain.pnum; p++) { 
			grid.set_T_last(p);
		}

		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(dynamic,1+sim.domain.pnum/sim.settings.thnum/64)
		for (int p = 0; p < sim.domain.pnum; p++) {
			if (grid.Calc_T(t, nodes, sim, true, p) < sim.material.T_liq){
				if (grid.get_T_last(p) >= sim.material.T_liq) {
					grid.Solidify(t, sim, p);
				}
			}
			else {
				#pragma omp atomic
				liq_num++;
			}
		}

		//Output data
		if (itert && (itert % sim.param.out_freq == 0)) {
			grid.Output(sim, Util::ZeroPadNumber(itert));
		}

		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_num)) {
			if (sim.param.out_freq != INT_MAX) {
				itert = ((itert / sim.param.out_freq) + 1) * sim.param.out_freq;
				grid.Output(sim, Util::ZeroPadNumber(itert));
			}
			break;
		}

		t += sim.param.dt;							//Increment time
		itert++;									//Increment time step counter
		liq_num = 0;
		Util::ClearNodes(nodes);
	}
	
	return;
}

void Run::Solidify_Volume(Grid& grid, const Simdat& sim) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time
	vector<omp_lock_t> lock(sim.domain.pnum);
	Util::SetLocks(lock, sim);

	// Vector of Liquid Point numbers
	vector<int> liq_pts;

	// Vector of Points to Reset Each Timestep
	vector<int> reset_pts;
	 
	// Integration information
	vector<Nodes> isegv_par;
	Nodes nodes;

	// Set initial time, iteration, and number of liquid points
	double t = 0.0;
	int itert = 0;
	int liq_num = 0;
	
	while (true) {
		Out::Progress(sim, itert);

		//Calculate Integration information
		Calc::Integrate_Parallel(nodes, sim, t, false);

		//Set T_calc_flag at all points to indicate that they have not yet been calculated
		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
		for (int r = 0; r < reset_pts.size(); r++) {
			const int p = reset_pts[r];
			grid.set_T_last(p);
			grid.set_T_calc_flag(false, p);
		}
		reset_pts.clear();

		//See which points from previous time step have solidified and do appropriate calculations
		vector<int> last_liq_pts = liq_pts;
		reset_pts = liq_pts;
		liq_pts.clear();

		//Check liquid points to see if they have solidified
		#pragma omp parallel num_threads(sim.settings.thnum)
		{
			vector<int> th_liq_pts;
			vector<int> th_reset_pts;
			#pragma omp for schedule(dynamic,1+last_liq_pts.size()/sim.settings.thnum/64)
			for (int it = 0; it < last_liq_pts.size(); it++) {
				const int p = last_liq_pts[it];
				if (grid.Calc_T(t, nodes, sim, true, p) < sim.material.T_liq) {
					grid.Solidify(t, sim, p);
				}
				else { 
					th_liq_pts.push_back(p); 
				}
			}
			#pragma omp critical
			liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
		}

		//Start search from points that are known to be liquid
		vector<int> test_pts = liq_pts;

		// Trace beam path between time steps and add relevant points to the test vector if:
		vector<int> bm_tr_pts;
		if (itert && t < sim.util.allScansEndTime + sim.param.dt) { 
			Melt::beam_trace(bm_tr_pts, grid, sim, t - sim.param.dt, t);
		}
		for (int it = 0; it < bm_tr_pts.size(); it++) {
			const int p = bm_tr_pts[it];
			reset_pts.push_back(p);
			if (grid.Calc_T(t, nodes, sim, true, p) >= sim.material.T_liq) {
				test_pts.push_back(p);
				liq_pts.push_back(p);
			}
		}

		//Iterative loop expanding from previously identified points to find melt pool boundary
		while (true) {
			if (!test_pts.size()) { break; }
			Melt::neighbor_check(test_pts, liq_pts, reset_pts, grid, lock, nodes, sim, t, false);
		}

		//Output data
		if (itert && (itert % sim.param.out_freq == 0)) {
			grid.Output(sim, Util::ZeroPadNumber(itert));
		}

		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_pts.size())) {
			if (sim.param.out_freq != INT_MAX) {
				itert = ((itert / sim.param.out_freq) + 1) * sim.param.out_freq;
				grid.Output(sim, Util::ZeroPadNumber(itert));
			}
			break;
		}

		t += sim.param.dt;		//Increment time
		itert++;				//Increment time step counter
		Util::ClearNodes(nodes);
	}

	return;
}

void Run::Solidify_Surface(Grid& grid, const Simdat& sim) {
	omp_set_nested(1);

	//Sets locks so only 1 thread can access a master point at the same time
	vector<omp_lock_t> lock(sim.domain.pnum);
	Util::SetLocks(lock, sim);

	// Vector of Liquid Point numbers
	vector<int> liq_pts;
	vector<int> depths(sim.domain.xnum*sim.domain.ynum, 0);
	vector<double> depths_max(sim.domain.xnum * sim.domain.ynum, 0);

	// Vector of Points to Reset Each Timestep
	vector<int> reset_pts;

	//Integration information
	Nodes nodes;
	Nodes nodes_last;

	int itert = 0, liq_num = 0;
	double t = 0.0;

	while (true) {
		Out::Progress(sim, itert);

		//Calculate Integration information
		Calc::Integrate_Parallel(nodes, sim, t, false);

		//Set T_calc_flag at all points to indicate that they have not yet been calculated
		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
		for (int r = 0; r < reset_pts.size(); r++) {
			const int p = reset_pts[r];
			if (grid.get_T(p) >= sim.material.T_liq) { grid.add_RDF_tm(t - sim.param.dt, p); }
			grid.set_T_last(p);
			grid.set_T_calc_flag(0, p);
		}
		reset_pts.clear();

		//See which points from previous time step have solidified and do appropriate calculations
		vector<int> last_liq_pts = liq_pts;
		reset_pts = liq_pts;
		liq_pts.clear();

		//Check liquid points on surface to see if they have solidified
		#pragma omp parallel num_threads(sim.settings.thnum)
		{
			vector<int> th_liq_pts;
			vector<int> th_reset_pts;
			#pragma omp for schedule(dynamic,1+last_liq_pts.size()/sim.settings.thnum/64)
			for (int it = 0; it < last_liq_pts.size(); it++) {
				const int p = last_liq_pts[it];
				if (grid.Calc_T(t, nodes, sim, true , p) < sim.material.T_liq){
					grid.Solidify(t, sim, p);
					// SOLIDIFY REST OF COLUMN //	
					const int i = grid.get_i(p);
					const int j = grid.get_j(p);
					const int dnum = i * sim.domain.ynum + j;
					const int depth = depths[dnum];
					for (int d = depth; d > 0; d--) {
						const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - d, sim);
						grid.Calc_T(t, nodes, sim, true, p_temp);
						if (d != depth) {
							const double T_last = grid.Calc_T(t - sim.param.dt, nodes_last, sim, false, p_temp);
							grid.set_T_last(T_last, p_temp);
						}
						grid.Solidify(t, sim, p_temp);
						th_reset_pts.push_back(p_temp);
					}
					depths[dnum] = 0;
				}
				else {th_liq_pts.push_back(p);}
			}
			#pragma omp critical
			{
				liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
				reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
			}
		}

		//Start search from points that are known to be liquid
		vector<int> test_pts = last_liq_pts;

		// Trace beam path between time steps and add relevant points to the test vector
		vector<int> bm_tr_pts;
		if (itert && t < sim.util.allScansEndTime + sim.param.dt) { 
			Melt::beam_trace(bm_tr_pts, grid, sim, t-sim.param.dt, t);
		}
		for (int it = 0; it < bm_tr_pts.size(); it++) {
			const int p = bm_tr_pts[it];
			reset_pts.push_back(p);
			if (grid.Calc_T(t, nodes, sim, true, p) >= sim.material.T_liq) {
				test_pts.push_back(p);
				liq_pts.push_back(p);
			}
		}

		//Iterative loop expanding from previously identified points to find melt pool boundary
		while (true) {
			if (!test_pts.size()) { break; }
			Melt::neighbor_check(test_pts, liq_pts, reset_pts, grid, lock, nodes, sim, t, true);
		}

		//Check the depths of the meltpool, solidify if needed
		Melt::calc_depth(depths, liq_pts, reset_pts, grid, nodes, nodes_last, sim, t);

		// Calculate the maximum depth under points
		Melt::calc_depth_max(depths, depths_max, liq_pts, grid, nodes, sim);

		//Output data
		if (itert && (itert % sim.param.out_freq == 0)) {
			grid.Output(sim, "Solidification."+Util::ZeroPadNumber(itert));
		}

		//Check if simulation is finished
		if (Util::sim_finish(t, sim, liq_pts.size())) {
			if (sim.param.out_freq != INT_MAX) {
				itert = ((itert / sim.param.out_freq) + 1) * sim.param.out_freq;
				grid.Output(sim, "Solidification."+Util::ZeroPadNumber(itert));
			}
			break;
		}

		t += sim.param.dt;			//Increment time
		itert++;					//Increment time step counter
		nodes_last = nodes;
		Util::ClearNodes(nodes);
	}
	return;
}