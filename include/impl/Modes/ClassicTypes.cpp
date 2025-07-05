/****************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC                                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of 3dThesis. 3dThesis is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include "impl/Modes/ClassicTypes.hpp"

#include "impl/Structs/DataStructs.hpp"
#include "impl/Calc/Calc.hpp"
#include "impl/Calc/Util.hpp"
#include "impl/Calc/Melt.hpp"
#include "impl/IO/Out.hpp"

namespace Thesis::impl{
	
	void Modes::Simulate(Grid& grid, const Simdat& sim) 
	{
		if (sim.param.mode == "Snapshots") { 
			Modes::Snapshots(grid, sim); 
		}
		else if (sim.param.mode == "Solidification") { 
			Modes::Solidify(grid, sim); 
		}
		else { 
			std::cout << "ERROR: Unrecognized mode: " << sim.param.mode << "\n"; 
		}
		return;
	}

	void Modes::Snapshots(Grid& grid, const Simdat& sim) 
	{
		if (sim.param.tracking == "None" || sim.domain.customPoints) {
			Modes::Snapshots_NoTracking(grid, sim);
		}
		else if (sim.param.tracking == "Volume") {
			Modes::Snapshots_Volume(grid, sim);
		}
		else if (sim.param.tracking == "Surface") {
			Modes::Snapshots_Volume(grid, sim);
		}
		else if (sim.param.tracking == "Geometry") {
			Modes::Snapshots_GeometryBounds(grid, sim);
		}
		else {
			std::cout << "ERROR: Unrecognized snapshots tracking: " << sim.param.tracking << "\n";
		}
		return;
	}

	void Modes::Snapshots_NoTracking(Grid& grid, const Simdat& sim) 
	{
		// Print Stuff
		int print_prog_last = 0;

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		// Nodes& nodes_last = quad.prev_nodes;
		// vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

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
			Calc::Integrate_Serial(nodes, start_seg, sim, t, false);
			
			// For all points, calculate temperature
			int p_tot = 0;
			#pragma omp parallel for num_threads(sim.settings.thnum) schedule(dynamic,1+sim.domain.pnum/sim.settings.thnum/128)
			for (int p = 0; p < sim.domain.pnum; p++) {
				grid.Calc_T(t, nodes, sim, true, p);
				#pragma omp atomic
				p_tot++;
				if (!omp_get_thread_num()) { Out::Point_Progress(print_prog_last, sim, p_tot); }
			}
			std::cout << "\n";

			// Output 
			grid.Output(sim, "Snapshot." + Util::ZeroPadNumber(itert, 2));
			
			// Clear quadrature nodes
			Util::ClearNodes(nodes);
		}

		return;
	}

	void Modes::Snapshots_Volume(Grid& grid, const Simdat& sim) {

		//Sets locks so only 1 thread can access a master point at the same time
		vector<omp_lock_t> lock(sim.domain.pnum);
		Util::SetLocks(lock, sim);

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		// Nodes& nodes_last = quad.prev_nodes;
		// vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

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
			Calc::Integrate_Serial(nodes, start_seg, sim, t, false);

			// Trace path from now until start of scan to seed initial points
			vector<int> bm_tr_pts;
			Melt::beam_trace(bm_tr_pts, grid, sim, 0, t);
			
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

	void Modes::Snapshots_GeometryBounds(Grid& grid, const Simdat& sim) {

		const vector<string> col_names = {
			"Length Raw (m)",
			"Width Raw (m)",
			"Origin X Raw (m)",
			"Origin Y Raw (m)",
			"Length Rotated (m)",
			"Width Rotated (m)",
			"Origin X Rotated (m)",
			"Origin Y Rotated (m)",
			"Percent Melted",
			"Depth (m)"
		};
		const double zres = sim.domain.zres;

		enum Columns {
			x = 0,
			y = 1,
			z = 2,
			T = 3,
			length_raw = 0,
			width_raw = 1,
			origin_x_raw = 2,
			origin_y_raw = 3,
			length_rotated = 4,
			width_rotated = 5,
			origin_x_rotated = 6,
			origin_y_rotated = 7,
			percent_melted = 8,
			depth = 9
		};

		vector<array<double, 10>> snaps;
		vector<double> times;
		vector<double> angles;
		vector<double> x_col;
		vector<double> y_col;
		double nan = std::numeric_limits<double>::quiet_NaN();

		//Sets locks so only 1 thread can access a master point at the same time
		vector<omp_lock_t> lock(sim.domain.pnum);
		Util::SetLocks(lock, sim);

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		// Nodes& nodes_last = quad.prev_nodes;
		// vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;
		
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

			// Find angle and update time and angle vectors
			double angle;
			const vector<path_seg>& path = sim.paths[0];
			int j = 0;
			for (const path_seg& seg : path){
				if (seg.seg_time >= t){
					break;
				}
				j++;
			}

			double x_0 = path[j-1].sx;
			double y_0 = path[j-1].sy;
			double x_1 = path[j].sx;
			double y_1 = path[j].sy;
			double dx = x_1 - x_0;
			double dy = y_1 - y_0;
			double dt = path[j].seg_time - path[j-1].seg_time;

			double beam_x = x_0 + (dx/dt)*(t - path[j-1].seg_time);
			double beam_y = y_0 + (dy/dt)*(t - path[j-1].seg_time);

			if (path[j].smode == 1){
				angle = 0.0;
				beam_x = x_1;
				beam_y = y_1;
			}
			else{
				angle = atan2(dy, dx);
				beam_x = x_0 + (dx/dt)*(t - path[j-1].seg_time);
				beam_y = y_0 + (dy/dt)*(t - path[j-1].seg_time);
			}
			angles.push_back(angle);
			times.push_back(t);
			x_col.push_back(beam_x);
			y_col.push_back(beam_y);

			// Output time
			std::cout << "Time: " << t << "\n";

			// Get quadrature nodes
			Calc::Integrate_Serial(nodes, start_seg, sim, t, false);

			// Trace path from now until start of scan to seed initial points
			vector<int> bm_tr_pts;
			if (itert) { Melt::beam_trace(bm_tr_pts, grid, sim, 0, t); }
			
			// Calculate temperature of seed points and see what is liquid
			vector<int> test_pts = liq_pts;
			liq_pts.clear();
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

			// Update snapshot data
			vector<vector<double>> pool;
			for (int p : liq_pts){
				pool.push_back({grid.get_x(p), grid.get_y(p), grid.get_z(p), grid.get_T(p)});
			}
			if (pool.size() == 0){snaps.push_back({0,0,nan,nan,0,0,nan,nan,0,0});}
			else{
				array<double, 10> row;
				vector<vector<double>> df_rot = Util::rotateField(pool, angle, x, y);
				array<double, 4> stats = Util::getLengthWidthOrigin(pool, zres, x, y);
				for (int k = 0; k < 4; k++){
					row[k] = stats[k];
				}
				stats = Util::getLengthWidthOrigin(df_rot, zres, x, y);
				for (int k = 0; k < 4; k++){
					row[k + 4] = stats[k];
				}
				double length = row[length_rotated];
				double width = row[width_rotated];
				double depth = Util::getMax(pool, z) - Util::getMin(pool, z);
				if (depth == 0){depth = 0.5 * zres;}
				row[8] = Util::getPerBoxMelted(pool, length, width, zres);
				row[9] = depth;
				snaps.push_back(row);
			}		
			// Clear integration segments
			Util::ClearNodes(nodes);
		}

		// Output data file
		std::ofstream datafile;
		datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		string out_file = sim.files.dataDir + "/" + "snapshot_data.csv";
		try {
			datafile.open(out_file.c_str());
			datafile << "Time (s),Scan Angle (rad)";
			for (int i = 0; i < col_names.size(); i++){
				datafile << "," + col_names[i];
			}
			datafile << ",Beam X,Beam Y";
			datafile << "\n";
			for (int i = 0; i < snaps.size(); i++){
				array<double, 10> row = snaps[i];
				datafile << times[i];
				datafile << ",";
				datafile << angles[i];
				for (int j = 0; j < row.size(); j++){
					datafile << ",";
					datafile << row[j];
				}
				datafile << "," << x_col[i] << "," << y_col[i];
				datafile << "\n";
			}
		}
		catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
		datafile.close();
		return;
	}

	void Modes::Solidify(Grid& grid, const Simdat& sim){
		if (sim.param.tracking == "None" || sim.domain.customPoints) {
			Modes::Solidify_NoTracking(grid, sim);
		}
		else if (sim.param.tracking == "Volume") {
			Modes::Solidify_Volume(grid, sim);
		}
		else if (sim.param.tracking == "Surface") {
			Modes::Solidify_Surface(grid, sim);
		}
		else if (sim.param.tracking == "Stork"){
			Modes::Stork(grid, sim);
		}
		else {
			std::cout << "ERROR: Unrecognized solidfication tracking: " << sim.param.tracking << "\n";
		}
		return;
	}

	void Modes::Solidify_NoTracking(Grid& grid, const Simdat& sim) {
		
		// Print Stuff
		int print_prog_last = 0;
		
		int itert = 0;
		int liq_num = 0;
		double t = 0.0;

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		// Nodes& nodes_last = quad.prev_nodes;
		vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

		#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
		for (int p = 0; p < sim.domain.pnum; p++) { 
			grid.set_T(sim.material.T_init, p);
		}

		while (true) {
			Out::Progress(print_prog_last, sim, itert);
			
			Calc::Integrate_Parallel(nodes, parNodes, start_seg, sim, t, false);

			#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
			for (int p = 0; p < sim.domain.pnum; p++) { 
				grid.set_T_last(p);
			}

			#pragma omp parallel for num_threads(sim.settings.thnum) schedule(dynamic,1+sim.domain.pnum/sim.settings.thnum/64)
			for (int p = 0; p < sim.domain.pnum; p++) {
				if (grid.Calc_T(t, nodes, sim, true, p) < sim.material.T_liq){
					if (grid.get_T_last(p) >= sim.material.T_liq) {
						grid.Solidify(start_seg, t, sim, p);
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

	void Modes::Solidify_Volume(Grid& grid, const Simdat& sim) {

		// Print Stuff
		int print_prog_last = 0;

		//Sets locks so only 1 thread can access a master point at the same time
		vector<omp_lock_t> lock(sim.domain.pnum);
		Util::SetLocks(lock, sim);

		// Vector of Liquid Point numbers
		vector<int> liq_pts;

		// Vector of Points to Reset Each Timestep
		vector<int> reset_pts;
		
		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		// Nodes& nodes_last = quad.prev_nodes;
		vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

		// Set initial time, iteration, and number of liquid points
		double t = 0.0;
		int itert = 0;
		
		while (true) {
			Out::Progress(print_prog_last, sim, itert);

			//Calculate Integration information
			Calc::Integrate_Parallel(nodes, parNodes, start_seg, sim, t, false);

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
						grid.Solidify(start_seg, t, sim, p);
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

	void Modes::Solidify_Surface(Grid& grid, const Simdat& sim) {

		// Print Stuff
		int print_prog_last = 0;

		//Sets locks so only 1 thread can access a master point at the same time
		vector<omp_lock_t> lock(sim.domain.pnum);
		Util::SetLocks(lock, sim);

		// Vector of Liquid Point numbers
		vector<int> liq_pts;
		vector<int> depths(sim.domain.xnum*sim.domain.ynum, 0);
		vector<double> depths_max(sim.domain.xnum * sim.domain.ynum, 0);

		// Vector of Points to Reset Each Timestep
		vector<int> reset_pts;

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		Nodes& nodes_last = quad.prev_nodes;
		vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

		int itert = 0;
		double t = 0.0;

		while (true) {
			Out::Progress(print_prog_last, sim, itert);

			//Calculate Integration information
			Calc::Integrate_Parallel(nodes, parNodes, start_seg, sim, t, false);

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
						grid.Solidify(start_seg, t, sim, p);
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
							grid.Solidify(start_seg, t, sim, p_temp);
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
			Melt::calc_depth(depths, liq_pts, reset_pts, grid, start_seg, nodes, nodes_last, sim, t);

			// Calculate the maximum depth under points
			Melt::calc_depth_max(depths, depths_max, liq_pts, grid, nodes, sim);

			// Calculate meltpool dimensions
			if (!sim.mpi)
				Melt::calc_mp_info(depths, grid, sim, t);
			else if (sim.print)
				std::cout << "Cannot export MeltPool Statistics with MPI domain decomposition." << std::endl;

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

	void Modes::Stork(Grid& grid, const Simdat& sim) {

		// Print Stuff
		int print_prog_last = 0;

		//Sets locks so only 1 thread can access a master point at the same time
		vector<omp_lock_t> lock(sim.domain.pnum);
		Util::SetLocks(lock, sim);

		// Vector of Liquid Point numbers
		vector<int> liq_pts;
		vector<int> depths(sim.domain.xnum*sim.domain.ynum, 0);
		vector<double> depths_max(sim.domain.xnum * sim.domain.ynum, 0);

		// Vector of Points to Reset Each Timestep
		vector<int> reset_pts;

		// Integration information
		QuadDat quad;
		Nodes& nodes = quad.cur_nodes;
		Nodes& nodes_last = quad.prev_nodes;
		vector<Nodes>& parNodes = quad.par_nodes;
		vector<int>& start_seg = quad.start_seg;

		// Grid of liquid info (cells)
		const int c_xnum = (sim.domain.xnum-1);
		const int c_ynum = (sim.domain.ynum-1);
		const int c_znum = (sim.domain.znum-1);
		const int c_pnum = c_xnum*c_ynum*c_znum;
		vector<uint8_t> c_alpha(c_pnum, static_cast<uint8_t>(0));
		vector<uint8_t> c_beta(c_pnum, static_cast<uint8_t>(0));

		// Info for "vertices"
		vector<double> T_alpha(sim.domain.pnum, sim.material.T_init);
		vector<double> T_beta(sim.domain.pnum, sim.material.T_init);
		vector<uint8_t> T_calc_alpha(sim.domain.pnum, static_cast<uint8_t>(1));
		vector<uint8_t> T_calc_beta(sim.domain.pnum, static_cast<uint8_t>(1));
		vector<uint8_t> isLiq(sim.domain.pnum, static_cast<uint8_t>(0));

		// Get vectors to store data into
		vector<uint32_t>& idxs = grid.get_RRDF_idxs();
		vector<double>& ts = grid.get_RRDF_ts();
		vector<double>& Ts = grid.get_RRDF_Ts();

		// Vector to store relevant C's in
		vector<int> c_relevant;
		c_relevant.reserve(c_pnum);

		int itert = 0;
		double t = 0.0;

		while (true) {
			
			Out::Progress(print_prog_last, sim, itert);

			// Determine which c is current and previous
			vector<uint8_t>& c_prev = ((itert%2) ? c_alpha : c_beta);
			vector<uint8_t>& c_cur = ((itert%2) ? c_beta : c_alpha);
			vector<double>& T_prev = ((itert%2) ? T_alpha : T_beta);
			vector<double>& T_cur = ((itert%2) ? T_beta : T_alpha);
			vector<uint8_t>& T_calc_prev = ((itert%2) ? T_calc_alpha : T_calc_beta);
			vector<uint8_t>& T_calc_cur = ((itert%2) ? T_calc_beta : T_calc_alpha);

			// Reset current T_calc flags
			std::fill(c_cur.begin(), c_cur.end(), static_cast<uint8_t>(0));
			std::fill(T_calc_cur.begin(), T_calc_cur.end(), static_cast<uint8_t>(0));
			std::fill(isLiq.begin(),isLiq.end(),static_cast<uint8_t>(0));
			c_relevant.clear();

			//Set T_calc_flag at all points to indicate that they have not yet been calculated
			#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
			for (int r = 0; r < reset_pts.size(); r++) {
				const int p = reset_pts[r];
				grid.set_T_last(p);
				grid.set_T_calc_flag(0, p);
			}
			reset_pts.clear();

			// Set surface liquid points from previous time to have their Temperatures checked
			vector<int> last_liq_pts = liq_pts;
			liq_pts.clear();

			// Calculate Integration information
			Calc::Integrate_Parallel(nodes, parNodes, start_seg, sim, t, false);

			// Check liquid points on surface to see if they have solidified
			#pragma omp parallel num_threads(sim.settings.thnum)
			{
				vector<int> th_liq_pts;
				vector<int> th_reset_pts;
				#pragma omp for schedule(static)
				for (int it = 0; it < last_liq_pts.size(); it++) {
					const int p = last_liq_pts[it];
					// Calculate Temperature
					T_cur[p] = grid.Calc_T(t, nodes, sim, true, p);
					T_calc_cur[p] = 1;
					// Push back points to reset
					th_reset_pts.push_back(p);
					// If it is still molten, add it to liquid points
					if (T_cur[p] >= sim.material.T_liq)
					{
						// Add to surface liquid points to being checking
						th_liq_pts.push_back(p);	
					}
					// If the previously molten surface point has solidified, then every point below it must have as well
					else 
					{
						// Get the <i,j> and depth
						const int i = grid.get_i(p);
						const int j = grid.get_j(p);
						const int dnum = i * sim.domain.ynum + j;
						// Set depth to zero
						depths[dnum] = 0;
					}
				}
				// Concatenate all liquid points together
				#pragma omp critical
				{
					reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
					liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
				}
			}

			// Start search from points that are known to be liquid
			vector<int> test_pts = liq_pts;

			// Trace beam path between time steps and add relevant points to the test vector
			vector<int> bm_tr_pts;
			if (itert && t < sim.util.allScansEndTime + sim.param.dt) { 
				Melt::beam_trace(bm_tr_pts, grid, sim, t-sim.param.dt, t);
			}

			// For all beam trace points, if they are molten, add them to the total liquid points and points to start checking surface neighbors
			for (int it = 0; it < bm_tr_pts.size(); it++) {
				const int p = bm_tr_pts[it];
				reset_pts.push_back(p);
				// Calculate Temperature
				if (T_calc_cur[p]==0){
					T_cur[p] = grid.Calc_T(t, nodes, sim, true , p);
					T_calc_cur[p] = 1;
					if (T_cur[p] >= sim.material.T_liq) {
						test_pts.push_back(p);
						liq_pts.push_back(p);
					}
				}
			}

			// Iterative loop expanding from previously identified points to find melt pool boundary
			while (true) {
				if (!test_pts.size()) { break; }
				vector<int> test_tmp;
				//identify neighbors of liquid points on surface
				#pragma omp parallel num_threads(sim.settings.thnum)
				{
					vector<int> th_liq_pts;
					vector<int> th_test_tmp;
					vector<int> th_reset_pts;
					//Find neighbors in test_pts for checking
					#pragma omp for schedule(static)
					for (int it = 0; it < test_pts.size(); it++) {
						//Get i, j, k location of current point, then construct array of neighbors
						const int p = test_pts[it];
						const int i = grid.get_i(p);
						const int j = grid.get_j(p);
						const int k = grid.get_k(p);
						// For neighbors
						for (int di=0;di<3;di++){
							for (int dj=0;dj<3;dj++){
								const int ni = i + di - 1;
								const int nj = j + dj - 1;
								const bool i_ob = (ni<0 || ni>sim.domain.xnum-1);
								const bool j_ob = (nj<0 || nj>sim.domain.ynum-1);
								if (i_ob || j_ob) { continue;}
								const int p_temp = Util::ijk_to_p(ni, nj, k, sim); 
								bool Tflag = false;
								omp_set_lock(&(lock[p_temp]));
								if (!grid.get_T_calc_flag(p_temp)) {
									grid.set_T_calc_flag(true, p_temp);
									Tflag = true;
								}
								omp_unset_lock(&(lock[p_temp]));
								if (Tflag) {
									T_cur[p_temp] = grid.Calc_T(t, nodes, sim, true , p_temp);
									T_calc_cur[p_temp] = 1;
									if (T_cur[p_temp] >= sim.material.T_liq) {
										th_liq_pts.push_back(p_temp);
										th_test_tmp.push_back(p_temp);
									}
									th_reset_pts.push_back(p_temp);
								}
							}
						}
					}
					#pragma omp critical
					{
						liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
						test_tmp.insert(test_tmp.end(), th_test_tmp.begin(), th_test_tmp.end());
						reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
					}
				}
				test_pts.clear();
				test_pts = test_tmp;
			}
			
			// Find the depths at each <i,j> of the meltpool
			#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static) //dynamic,1+last_liq_pts.size()/sim.settings.thnum/64)
			for (int it = 0; it < liq_pts.size(); it++) {
				// Get <i,j> of the point
				const int p = liq_pts[it];
				const int i = grid.get_i(p);
				const int j = grid.get_j(p);
				// Get starting depth
				const int dnum = i * sim.domain.ynum + j;
				int depth = depths[dnum];
				bool wentDeeper = false;
				while (true) {
					// Increment depth
					depth++;
					// If at bottom of domain, break and display message
					if (depth == sim.domain.znum) { 
						throw std::runtime_error("Error -> 3DThesis -> Maximum Depth Exceeded.");
					}
					// Find point number
					const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth, sim);
					// If the point is solid, go to next loop
					T_cur[p_temp] = grid.Calc_T(t, nodes, sim, false, p_temp);
					T_calc_cur[p_temp] = 1;
					if (T_cur[p_temp] < sim.material.T_liq) {break;}
					// Otherwise, set a flag saying we sucessfully went deeper
					else { wentDeeper = true; }
				}
				while (true) {
					// Decrement depth
					depth--;
					// If at the bottom or we already went deeper, break
					if (depth == 0 || wentDeeper == true) { break; }
					// Find point number
					const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth, sim);
					// If point is liquid, break
					T_cur[p_temp] = grid.Calc_T(t, nodes, sim, false, p_temp);
					T_calc_cur[p_temp] = 1;
					if (T_cur[p_temp] >= sim.material.T_liq) { break; }
				}
				// If at max depth, decrement once
				if (depth == sim.domain.znum - 1) { depth--; }
				// Store depth
				depths[dnum] = depth;
				// For all points in the depth, set to be liquid
				for (int d = 0; d <= depth; d++) {
					const int k = sim.domain.znum - 1 - d;
					const int p_temp = Util::ijk_to_p(i, j, k, sim);
					isLiq[p_temp] = true;
				}
			}

			if (last_liq_pts.size() || liq_pts.size()){
				// Loop over and compress to C
				#pragma omp parallel num_threads(sim.settings.thnum)
				{
					vector<int> th_c_relevant;
					th_c_relevant.reserve(c_pnum);
					#pragma omp for schedule(static)//dynamic,16384)//static)
					for (int c = 0; c < c_pnum; c++){
						// Get <i,j,k> of the min point (which is <ci,cj,ck>)
						const int ci = (c/c_ynum)/(c_znum);
						const int cj = (c/c_znum)%(c_ynum);
						const int ck = (c%c_znum);
						// Loop over to cell vertices find number of liquid points
						for (int di=0;di<2;di++){
							for (int dj=0;dj<2;dj++){
								for (int dk=0;dk<2;dk++){
									const int i = ci + di;
									const int j = cj + dj;
									const int k = ck + dk;
									// Get point number
									const int p = Util::ijk_to_p(i, j, k, sim);
									// Add to counter
									c_cur[c] += isLiq[p];
								}
							}
						}
						// If we should be saving temperatures
						if (c_cur[c]!=c_prev[c] || c_prev[c]%8 || c_cur[c]%8){
							th_c_relevant.push_back(c);
						}
					}
					#pragma omp critical
					{
						c_relevant.insert(c_relevant.end(), th_c_relevant.begin(), th_c_relevant.end());
					}
				}
			}

			if (c_relevant.size()){
				// Put data in a struct, calculate as needed
				#pragma omp parallel num_threads(sim.settings.thnum)
				{
					vector<uint32_t> th_idxs;
					vector<double> th_ts;
					vector<double> th_Ts;
					#pragma omp for schedule(static)
					for (int i=0; i<c_relevant.size();i++){
						// get <i,j,k> of cell
						const int c = c_relevant[i];
						const int ci = (c/c_ynum)/(c_znum);
						const int cj = (c/c_znum)%(c_ynum);
						const int ck = (c%c_znum);
						// Store index data
						th_idxs.push_back(ci);
						th_idxs.push_back(cj);
						th_idxs.push_back(ck);
						// Store time data
						th_ts.push_back(t - sim.param.dt);
						th_ts.push_back(t);
						// Loop over cell vertices to calculate temperatures
						int n=0;
						vector<double> temp_Ts(16);
						for (int di=0;di<2;di++){
							for (int dj=0;dj<2;dj++){
								for (int dk=0;dk<2;dk++){
									const int i = ci + di;
									const int j = cj + dj;
									const int k = ck + dk;
									// Get point number
									const int p = Util::ijk_to_p(i, j, k, sim);
									// "Atomically" determine temperatures to be calculated
									bool T_calc_p_prev = false;
									bool T_calc_p_cur = false;
									omp_set_lock(&(lock[p]));
									{
										if (T_calc_prev[p] == 0) {
											T_calc_prev[p] = 1;
											T_calc_p_prev = true;
										}
										if (T_calc_cur[p] == 0) {
											T_calc_cur[p] = 1;
											T_calc_p_cur = true;
										}
										// Calculate temperatures
										if (T_calc_p_prev){T_prev[p] = grid.Calc_T(t - sim.param.dt, nodes_last, sim, false, p);}
										if (T_calc_p_cur){T_cur[p] = grid.Calc_T(t, nodes, sim, false, p);}
									}
									omp_unset_lock(&(lock[p]));
									// Set temperatures
									temp_Ts[n] = T_prev[p];
									temp_Ts[n+8] = T_cur[p];
									// Increment neighbor number
									n++;
								}
							}
						}
						th_Ts.insert(th_Ts.end(), temp_Ts.begin(), temp_Ts.end());
					}
					// Concatenate data
					#pragma omp critical
					{
						idxs.insert(idxs.end(), th_idxs.begin(), th_idxs.end());
						ts.insert(ts.end(), th_ts.begin(), th_ts.end());
						Ts.insert(Ts.end(), th_Ts.begin(), th_Ts.end());
					}
				}
			}

			//Check if simulation is finished
			if (Util::sim_finish(t, sim, liq_pts.size())) {	break; }

			t += sim.param.dt;			//Increment time
			itert++;					//Increment time step counter
			nodes_last = nodes;
			Util::ClearNodes(nodes);
		}

		return;
	}
}

