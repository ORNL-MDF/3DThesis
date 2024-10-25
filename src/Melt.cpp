//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) - Transformer Resilience and Advanced Components (TRAC) Program.

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

#include "Melt.h"

#include <vector>
#include <omp.h>
#include <cmath>

using std::max;

void Melt::beam_trace(vector<int>& test_pts, Grid& grid, const Simdat& sim, const double t_start, const double t_end) {
	// For each path
	for (const vector<path_seg>& path : sim.paths) {
		
		// Find segment associated with starting time
		int seg_start = 0;
		while (path[seg_start].seg_time < t_start && (seg_start + 1) != path.size()) { seg_start++; }

		// Find segment associated with ending time
		int seg_end = 0;
		while (path[seg_end].seg_time < t_end && (seg_end + 1) != path.size()) { seg_end++; }

		// Variables for integer grid numbers
		int x_grid_num = 0, y_grid_num = 0, z_grid_num = 0;

		// Variables for if the domain the grid number is out of bounds
		int x_flat = 0, y_flat = 0, z_flat = 0;

		// For all path segments between starting and ending segments
		for (int seg = max(seg_start - 1, 0); seg <= seg_end; seg++) {
			
			// If {xnum,ynum,znum}>1, find the grid numbers for the segments {x,y,z}
			if (sim.domain.xnum - 1) { x_grid_num = (int)((path[seg].sx - sim.domain.xmin) / sim.domain.xres); x_flat = 1; }
			if (sim.domain.ynum - 1) { y_grid_num = (int)((path[seg].sy - sim.domain.ymin) / sim.domain.yres); y_flat = 1; }
			if (sim.domain.znum - 1) { z_grid_num = (int)((path[seg].sz - sim.domain.zmin) / sim.domain.zres); z_flat = 1; }

			// If out of bounds, set the the end and make flat
			if (x_grid_num < 0) { x_grid_num = 0; x_flat = 0; }
			else if (x_grid_num >= (sim.domain.xnum - 1)) { x_grid_num = sim.domain.xnum - 1; x_flat = 0; }
			if (y_grid_num < 0) { y_grid_num = 0; y_flat = 0; }
			else if (y_grid_num >= (sim.domain.ynum - 1)) { y_grid_num = sim.domain.ynum - 1; y_flat = 0; }
			if (z_grid_num < 0) { z_grid_num = 0; z_flat = 0; }
			else if (z_grid_num >= (sim.domain.znum - 1)) { z_grid_num = sim.domain.znum - 1; z_flat = 0; }

			// For each direction, add all points in that square
			for (int a = 0; a <= z_flat; a++) {
				for (int b = 0; b <= y_flat; b++) {
					for (int c = 0; c <= x_flat; c++) {
						const int p = (z_grid_num + a) + sim.domain.znum * (y_grid_num + b) + sim.domain.znum * sim.domain.ynum * (x_grid_num + c);
						if (!grid.get_T_calc_flag(p))
						{
							test_pts.push_back(p);
							grid.set_T_calc_flag(true, p);
						}
					}
				}
			}
		}

		// If the "current" segment is a line, also add the current point
		if (path[seg_end].smode == 0 && t_end<path.back().seg_time) {
			int_seg current_beam = Util::GetBeamLoc(t_end, seg_end, path, sim);
			if (sim.domain.xnum - 1) { x_grid_num = (int)((current_beam.xb - sim.domain.xmin) / sim.domain.xres); x_flat = 1; }
			if (sim.domain.ynum - 1) { y_grid_num = (int)((current_beam.yb - sim.domain.ymin) / sim.domain.yres); y_flat = 1; }
			if (sim.domain.znum - 1) { z_grid_num = (int)((current_beam.zb - sim.domain.zmin) / sim.domain.zres); z_flat = 1; }

			if (x_grid_num < 0) { x_grid_num = 0; x_flat = 0; }
			else if (x_grid_num >= (sim.domain.xnum - 1)) { x_grid_num = sim.domain.xnum - 1; x_flat = 0; }
			if (y_grid_num < 0) { y_grid_num = 0; y_flat = 0; }
			else if (y_grid_num >= (sim.domain.ynum - 1)) { y_grid_num = sim.domain.ynum - 1; y_flat = 0; }
			if (z_grid_num < 0) { z_grid_num = 0; z_flat = 0; }
			else if (z_grid_num >= (sim.domain.znum - 1)) { z_grid_num = sim.domain.znum - 1; z_flat = 0; }

			for (int a = 0; a <= z_flat; a++) {
				for (int b = 0; b <= y_flat; b++) {
					for (int c = 0; c <= x_flat; c++) {
						const int p = (z_grid_num + a) + sim.domain.znum * (y_grid_num + b) + sim.domain.znum * sim.domain.ynum * (x_grid_num + c);
						if (!grid.get_T_calc_flag(p)) {
							test_pts.push_back(p);
							grid.set_T_calc_flag(true, p);
						}
					}
				}
			}
		}
	}
	return;
}

void Melt::neighbor_check(vector<int>& test_pts, vector<int>& liq_pts, vector<int>& reset_pts, Grid& grid, vector<omp_lock_t>& lock, const Nodes& nodes, const Simdat& sim, const double t, const bool surface_only) {
	vector<int> test_tmp;
	//identify neighbors of liquid points ONLY ON SURFACE
	#pragma omp parallel num_threads(sim.settings.thnum)
	{
		vector<int> th_liq_pts;
		vector<int> th_test_tmp;
		vector<int> th_reset_pts;
		//Find neighbors in test_pts for checking
		#pragma omp for schedule(dynamic,1+test_pts.size()/sim.settings.thnum/64)
		for (int it = 0; it < test_pts.size(); it++) {
			//Get i, j, k location of current point, then construct array of neighbors
			const int p = test_pts[it];
			const int i = grid.get_i(p);
			const int j = grid.get_j(p);
			const int k = grid.get_k(p);

			if (surface_only){
				for (int di=-1;di<=1;di++){
					for (int dj=-1;dj<=1;dj++){
						const int ni = i + di;
						const int nj = j + dj;
						if (ni<0 || (ni+1)>=sim.domain.xnum || nj<0 || (nj+1)>=sim.domain.ynum){
							continue;
						}
						const int p_temp = Util::ijk_to_p(ni, nj, k, sim); 

						bool Tflag = false;
						omp_set_lock(&(lock[p_temp]));
						if (!grid.get_T_calc_flag(p_temp)) {
							grid.set_T_calc_flag(true, p_temp);
							Tflag = true;
						}
						omp_unset_lock(&(lock[p_temp]));
						if (Tflag) {
							if (grid.Calc_T(t, nodes, sim, true, p_temp) >= sim.material.T_liq) {
								th_liq_pts.push_back(p_temp);
								th_test_tmp.push_back(p_temp);
							}
							th_reset_pts.push_back(p_temp);
						}
					}
				}
			}
			else{
				for (int di=-1;di<=1;di++){
					for (int dj=-1;dj<=1;dj++){
						for (int dk=-1;dk<1;dk++){
							const int ni = i + di;
							const int nj = j + dj;
							const int nk = k + dk;
							if (ni<0 || (ni+1)>=sim.domain.xnum || nj<0 || (nj+1)>=sim.domain.ynum || nk<0){
								continue;
							}
							const int p_temp = Util::ijk_to_p(ni, nj, nk, sim); 

							bool Tflag = false;
							omp_set_lock(&(lock[p_temp]));
							if (!grid.get_T_calc_flag(p_temp)) {
								grid.set_T_calc_flag(true, p_temp);
								Tflag = true;
							}
							omp_unset_lock(&(lock[p_temp]));
							if (Tflag) {
								if (grid.Calc_T(t, nodes, sim, true, p_temp) >= sim.material.T_liq) {
									th_liq_pts.push_back(p_temp);
									th_test_tmp.push_back(p_temp);
								}
								th_reset_pts.push_back(p_temp);
							}
						}
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

void Melt::calc_depth(vector<int>& depths, vector<int>& liq_pts, vector<int>& reset_pts, Grid& grid, const Nodes& nodes, const Nodes& isegv_last, const Simdat& sim, const double t) {
	if (sim.domain.znum == 1) {return;} 
	#pragma omp parallel num_threads(sim.settings.thnum)
	{
		vector<int> th_reset_pts;
		#pragma omp for schedule(dynamic)
		for (int it = 0; it < liq_pts.size(); it++) {
			const int p = liq_pts[it];
			const int i = grid.get_i(p);
			const int j = grid.get_j(p);
			int dnum = i * sim.domain.ynum + j;
			int depth = depths[dnum];
			int flag = 1;
			while (true) {
				depth++;
				if (depth == sim.domain.znum) { break; }
				const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth, sim);
				th_reset_pts.push_back(p_temp);
				if (grid.Calc_T(t, nodes, sim, true, p_temp) < sim.material.T_liq) {break;}
				else { flag = 0; }
			}
			while (true) {
				depth--;
				if (depth == 0 || flag == 0) { break; }
				const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth, sim);
				th_reset_pts.push_back(p_temp);
				if (grid.Calc_T(t, nodes, sim, true, p_temp) < sim.material.T_liq) {
					if (depth != depths[dnum]) { 
						const double T_last = grid.Calc_T(t - sim.param.dt, isegv_last, sim, false, p_temp);
						grid.set_T_last(T_last, p_temp);
					}
					grid.Solidify(t, sim, p_temp);
				}
				else { 
					break; 
				}
			}
			if (depth == sim.domain.znum) { depth--; }
			depths[dnum] = depth;
		}
		#pragma omp critical
		{
			reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
		}
	}
	return;
}

void Melt::calc_depth_max(vector<int>& depths, vector<double>& depth_max, vector<int>& liq_pts, Grid& grid, const Nodes& nodes, const Simdat& sim) {
	if (sim.domain.znum == 1 || sim.output.depth == 0) { return; }
	#pragma omp parallel for num_threads(sim.settings.thnum) schedule(static)
	for (int it = 0; it < liq_pts.size(); it++) {
		const int p = liq_pts[it];
		const int i = grid.get_i(p);
		const int j = grid.get_j(p);
		int dnum = i * sim.domain.ynum + j;
		int depth_liq = depths[dnum];

		double depth;
		if (depth_liq == sim.domain.znum - 1) {
			depth = depth_liq;
		}
		else {
			int depth_sol = depth_liq + 1;

			const int p_liq = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth_liq, sim);
			const int p_sol = Util::ijk_to_p(i, j, sim.domain.znum - 1 - depth_sol, sim);
			const double T_liq = grid.get_T(p_liq);
			const double T_sol = grid.get_T(p_sol);

			depth = depth_liq + (sim.material.T_liq - T_liq) / (T_sol - T_liq);
		}

		// If depth is greater
		if (depth > depth_max[dnum]) {
			// Store maximum depth
			depth_max[dnum] = depth;
			// Set that points and all points above it to have that max depth
			for (int d = depth_liq; d >= 0; d--) {
				const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - d, sim);
				grid.set_depth(depth, p_temp);
			}
		}

	}
	return;
}

void Melt::calc_mp_info(const vector<int>& depths, Grid& grid, const Simdat& sim, const double t){
	// If not outputting, don't do
	if (sim.output.mp_stats == 0) { return;}
	
	// This enables each segment to have its own local pools
	static list<vector<int>> multi_local_liq_pools;
	static list<int> multi_segs;
	static list<vector<double>> multi_trigs;

	// This enables the starting search path segment to be quickly initialized each time. 
	static int seg(0);

	// Set reference just to first path
	const vector<path_seg>& path = sim.paths[0];

	// Keep incrementing up if t is greater than the end of the path segment but also below the end time of the scan
	while ((t > path[seg].seg_time) && (seg + 1 < path.size())) { 
		// Increment segment
		seg++; 
		// Add test points (if power is on)
		if (path[seg].sqmod>0){
			// What segment information should be used?
			int seg_pt;
			if (path[seg].smode==0){seg_pt=seg-1;}
			else{seg_pt=seg;}
			const int x_grid_num = (int)((path[seg_pt].sx - sim.domain.xmin) / sim.domain.xres);
			const int y_grid_num = (int)((path[seg_pt].sy - sim.domain.ymin) / sim.domain.yres);
			const int z_grid_num = (int)((path[seg_pt].sz - sim.domain.zmin) / sim.domain.zres);

			// Vector of test points
			vector<int> test_pts;
			for (int dx=0;dx<=1;dx++){
				for (int dy=0;dy<=1;dy++){
					const int p = (z_grid_num) + sim.domain.znum * (y_grid_num + dy) + sim.domain.znum * sim.domain.ynum * (x_grid_num + dx);
					test_pts.push_back(p);
				}
			}

			// Add to final stuff
			multi_local_liq_pools.push_back(test_pts);
			multi_segs.push_back(seg);	

			// Find angle and update time and angle vectors
			vector<double> trigs(2);
			double angle;
			const double x_0 = path[seg-1].sx;
			const double y_0 = path[seg-1].sy;
			const double x_1 = path[seg].sx;
			const double y_1 = path[seg].sy;
			const double dx = x_1 - x_0;
			const double dy = y_1 - y_0;
			if (path[seg].smode == 1){angle = 0.0;}
			else{angle = atan2(dy, dx);}
			trigs[0] = sin(angle);
			trigs[1] = cos(angle);
			multi_trigs.push_back(trigs);
		}
	}

	// If we are on the same seg (and line melt), add current position too
	if (multi_segs.size()){
		if (multi_segs.back()==seg && path[seg].smode==0){
			vector<int> local_test_pts = multi_local_liq_pools.back();
			int_seg current_beam = Util::GetBeamLoc(t, seg, path, sim);
			const int x_grid_num = (int)((current_beam.xb - sim.domain.xmin) / sim.domain.xres);
			const int y_grid_num = (int)((current_beam.yb - sim.domain.ymin) / sim.domain.yres);
			const int z_grid_num = (int)((current_beam.zb - sim.domain.zmin) / sim.domain.zres);

			// Vector of test points
			for (int dx=0;dx<=1;dx++){
				for (int dy=0;dy<=1;dy++){
					const int p = (z_grid_num) + sim.domain.znum * (y_grid_num + dy) + sim.domain.znum * sim.domain.ynum * (x_grid_num + dx);
					local_test_pts.push_back(p);
				}
			}

			multi_local_liq_pools.back() = local_test_pts;
		}
	}
	else{
		return;
	}
	
	// For all remaining segments
	auto seg_it = multi_segs.rbegin();
	auto liq_it = multi_local_liq_pools.rbegin();
	auto trig_it = multi_trigs.rbegin();
	while (seg_it != multi_segs.rend()){
		
		// Get local values (front to back)
		const int& local_seg = *seg_it;	
		vector<int>& local_test_pts = *liq_it;
		const vector<double>& local_trigs = *trig_it;
		const double& sin_angle = local_trigs[0];
		const double& cos_angle = local_trigs[1];
		
		// Iterative loop from current point to find local, liquid points
		vector<int> local_liq_pts;
		while (!local_test_pts.empty()){
			local_neighbor_check(local_test_pts, local_liq_pts, grid, sim);
		}

		// If no liquid points
		if (local_liq_pts.empty()){
			// Erase data using base iterators
			seg_it = list<int>::reverse_iterator(multi_segs.erase((++seg_it).base()));
			liq_it = list<vector<int>>::reverse_iterator(multi_local_liq_pools.erase((++liq_it).base()));
			trig_it = list<vector<double>>::reverse_iterator(multi_trigs.erase((++trig_it).base()));
		}
		else{
			// Set test points next iteration to current pool
			local_test_pts = local_liq_pts;

			// Calculate rotated x,y of liquid points
			double minRotX = std::numeric_limits<double>::max();
			double maxRotX= std::numeric_limits<double>::lowest();
			double minRotY = std::numeric_limits<double>::max();
			double maxRotY = std::numeric_limits<double>::lowest();
			for (const int& liq_pt:local_liq_pts){
				const double x = grid.get_x(liq_pt);
				const double y = grid.get_y(liq_pt);

				const double x_rot = x*cos_angle + y*sin_angle;
				minRotX = std::min(x_rot, minRotX);
				maxRotX = std::max(x_rot, maxRotX);

				const double y_rot = -x*sin_angle + y*cos_angle;
				minRotY = std::min(y_rot, minRotY);
				maxRotY = std::max(y_rot, maxRotY);
			}

			// Calculate width and depth
			const double width = maxRotY - minRotY;
			const double length = maxRotX - minRotX;
			
			// Get depth
			const double depth = sim.domain.xres * (*std::max_element(depths.begin(), depths.end()));

			// Now add to all relevant points
			for (const int& liq_pt:local_liq_pts){
				// Get i,j and <ij>
				const int i = grid.get_i(liq_pt);
				const int j = grid.get_j(liq_pt);
				int dnum = i * sim.domain.ynum + j;
				// For points in depth
				for (int d=0;d<=depths[dnum];d++){
					const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - d, sim);
					grid.set_mpWidth(width, p_temp);
					grid.set_mpLength(length, p_temp);
					grid.set_mpDepth(depth, p_temp);
				}
			}

			// Update iterators
			++seg_it;
			++liq_it;
			++trig_it;
		}		
	}
}

void Melt::local_neighbor_check(vector<int>& test_pts, vector<int>& liq_pts, Grid& grid, const Simdat& sim) {
	
	vector<int> test_tmp;
	
	//Find neighbors in test_pts for checking
	for (int it = 0; it < test_pts.size(); it++) {
		//Get i, j, k location of current point, then construct array of neighbors
		const int p = test_pts[it];
		const int i = grid.get_i(p);
		const int j = grid.get_j(p);
		const int k = grid.get_k(p);

		for (int di=-1;di<=1;di++){
			for (int dj=-1;dj<=1;dj++){
				const int ni = i + di;
				const int nj = j + dj;
				if (ni<0 || (ni+1)>=sim.domain.xnum || nj<0 || (nj+1)>=sim.domain.ynum){
					continue;
				}
				const int p_temp = Util::ijk_to_p(ni, nj, k, sim); 
				if (grid.get_T_calc_flag(p_temp) && grid.get_T(p_temp)>=sim.material.T_liq) {
					grid.set_T_calc_flag(false, p_temp);
					test_tmp.push_back(p_temp);
					liq_pts.push_back(p_temp);
				}
			}
		}
	}

	test_pts.clear();
	test_pts = test_tmp;
}