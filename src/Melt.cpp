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

			const int n = 1;
			std::vector <int> ijkminmax;
			ijkminmax.push_back(n - i);
			ijkminmax.push_back(i + 1 + n - sim.domain.xnum);
			ijkminmax.push_back(n - j);
			ijkminmax.push_back(j + 1 + n - sim.domain.ynum);
			ijkminmax.push_back(n - k);
			ijkminmax.push_back(k + 1 + n - sim.domain.znum);

			for (int temp = 0; temp < ijkminmax.size(); temp++) {
				if (ijkminmax[temp] < 0) { ijkminmax[temp] = 0; }
			}

			vector<int> nbs;
			int p_temp;
			bool Tflag = false;
			for (int di = -n + ijkminmax[0]; di <= (n - ijkminmax[1]); di++) {
				for (int dj = -n + ijkminmax[2]; dj <= (n - ijkminmax[3]); dj++) {
					for (int dk = -n + ijkminmax[4]; dk <= (n - ijkminmax[5]); dk++) {
						if (surface_only) { 
							p_temp = Util::ijk_to_p(i + di, j + dj, k, sim); 
						}
						else { 
							p_temp = Util::ijk_to_p(i + di, j + dj, k + dk, sim); 
						}
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
						Tflag = false;
						if (surface_only) { break; }
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

void Melt::calc_mp_info(const vector<int>& depths, const vector<int>& liq_pts, Grid& grid, const Simdat& sim, const double t){
	// If not outputting, don't do
	if (sim.output.mp_stats == 0) { return;}
	
	// This enables the starting search path segment to be quickly initialized each time. 
	static int seg(1);

	// Set reference just to first path
	const vector<path_seg>& path = sim.paths[0];

	// Keep incrementing up if t is greater than the end of the path segment but also below the end time of the scan
	while ((t > path[seg].seg_time) && (seg + 1 < path.size())) { seg++; }
	// Keep incrementing down if t is less than end of previous path segment
	while ((t < path[seg - 1].seg_time) && (seg - 1 > 0)) { seg--; }

	// Find angle and update time and angle vectors
	double angle;
	double x_0 = path[seg-1].sx;
	double y_0 = path[seg-1].sy;
	double x_1 = path[seg].sx;
	double y_1 = path[seg].sy;
	double dx = x_1 - x_0;
	double dy = y_1 - y_0;
	if (path[seg].smode == 1){angle = 0.0;}
	else{angle = atan2(dy, dx);}
	
	// Calculate depths and widths
	if (liq_pts.size() != 0 && path[seg].sqmod>0){

		// Calculate rotated x,y of liquid points
		double minRotX = std::numeric_limits<double>::max();
		double maxRotX= std::numeric_limits<double>::lowest();
		double minRotY = std::numeric_limits<double>::max();
		double maxRotY = std::numeric_limits<double>::lowest();
		for (const int& liq_pt:liq_pts){
			const double x = grid.get_x(liq_pt);
			const double y = grid.get_y(liq_pt);

			const double x_rot = x*cos(angle) + y*sin(angle);
			minRotX = std::min(x_rot, minRotX);
	        maxRotX = std::max(x_rot, maxRotX);

			const double y_rot = -x*sin(angle) + y*cos(angle);
			minRotY = std::min(y_rot, minRotY);
	        maxRotY = std::max(y_rot, maxRotY);
		}

		// Calculate width and depth
		const double width = maxRotY - minRotY;
		const double length = maxRotX - minRotX;
		
		// Get depth
		const double depth = sim.domain.xres * (*std::max_element(depths.begin(), depths.end()));

		// Now add to all relevant points
		for (const int& liq_pt:liq_pts){
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
	}
}

// void findEigenValuesAndVectors(double xx, double xy, double yy, double& lambda1, double& lambda2, std::pair<double, double>& vec1, std::pair<double, double>& vec2) {
//     // Calculate eigenvalues
//     double trace = xx + yy;
//     double determinant = xx * yy - xy * xy;
//     double sqrtTerm = std::sqrt(trace * trace - 4 * determinant);

//     lambda1 = (trace + sqrtTerm) / 2.0;
//     lambda2 = (trace - sqrtTerm) / 2.0;

//     // Calculate eigenvectors
//     if (xy != 0) {
//         // For lambda1
//         vec1 = {lambda1 - yy, xy};
//         // For lambda2
//         vec2 = {lambda2 - yy, xy};
//     } else {
//         // Special case when b == 0 (diagonal matrix)
//         vec1 = {1, 0}; // Eigenvector along x-axis
//         vec2 = {0, 1}; // Eigenvector along y-axis
//     }

//     // Normalize the eigenvectors
//     double norm1 = std::sqrt(vec1.first * vec1.first + vec1.second * vec1.second);
//     vec1.first /= norm1;
//     vec1.second /= norm1;

//     double norm2 = std::sqrt(vec2.first * vec2.first + vec2.second * vec2.second);
//     vec2.first /= norm2;
//     vec2.second /= norm2;
// }

// void Melt::calc_mp_info(const vector<int>& depths, const vector<int>& liq_pts, Grid& grid, const Simdat& sim){
// 	// If not outputting, don't do
// 	if (sim.output.mp_stats == 0) { return;}
    
// 	// Compute means
// 	double x_mean = 0.0;
// 	double y_mean = 0.0;
// 	for (const int& liq_pt:liq_pts){
// 		x_mean += grid.get_x(liq_pt);
// 		y_mean += grid.get_y(liq_pt);
// 	}
// 	x_mean /= liq_pts.size();
//     y_mean /= liq_pts.size();

// 	// Compute the covariance matrix elements
//     double cov_xx = 0.0;
// 	double cov_yy = 0.0;
// 	double cov_xy = 0.0;
// 	for (const int& liq_pt:liq_pts){
// 		double x_diff = grid.get_x(liq_pt); - x_mean;
//         double y_diff = grid.get_y(liq_pt); - y_mean;
//         cov_xx += x_diff * x_diff;
//         cov_yy += y_diff * y_diff;
//         cov_xy += x_diff * y_diff;
// 	}
//     cov_xx /= liq_pts.size();
//     cov_yy /= liq_pts.size();
//     cov_xy /= liq_pts.size();

// 	// Get eigenvalues and eigenvectors
// 	double eig_1, eig_2;
// 	std::pair<double, double> vec_1, vec_2;
// 	findEigenValuesAndVectors(cov_xx, cov_xy, cov_yy, eig_1, eig_2, vec_1, vec_2);

// 	// The eigenvector corresponding to the largest eigenvalue is the major axis
// 	std::pair<double, double> minorAxis, majorAxis;
//     if (eig_1 > eig_2) {
//         majorAxis = vec_1;
//         minorAxis = vec_2;
//     } else {
//         majorAxis = vec_2;
//         minorAxis = vec_1;
//     }

// 	// Initialize min and max projections for major and minor axes
//     double minMajorProj = std::numeric_limits<double>::max();
//     double maxMajorProj = std::numeric_limits<double>::lowest();
//     double minMinorProj = std::numeric_limits<double>::max();
//     double maxMinorProj = std::numeric_limits<double>::lowest();

//     // Project each point onto the major and minor axes
//     for (const int& liq_pt:liq_pts) {
//         double x_diff = grid.get_x(liq_pt) - x_mean;
//         double y_diff = grid.get_y(liq_pt) - y_mean;

//         // Projection onto the major axis
//         double majorProj = x_diff * majorAxis.first + y_diff * majorAxis.second;
//         minMajorProj = std::min(minMajorProj, majorProj);
//         maxMajorProj = std::max(maxMajorProj, majorProj);

//         // Projection onto the minor axis
//         double minorProj = x_diff * minorAxis.first + y_diff * minorAxis.second;
//         minMinorProj = std::min(minMinorProj, minorProj);
//         maxMinorProj = std::max(maxMinorProj, minorProj);
//     }

//     // Calculate width and depth
// 	const double width = maxMinorProj - minMinorProj;
//     const double length = maxMajorProj - minMajorProj;
	
// 	// Get depth
// 	const double depth = sim.domain.xres * (*std::max_element(depths.begin(), depths.end()));

// 	// Now add to all relevant points
// 	for (const int& liq_pt:liq_pts){
// 		// Get i,j and <ij>
// 		const int i = grid.get_i(liq_pt);
// 		const int j = grid.get_j(liq_pt);
// 		int dnum = i * sim.domain.ynum + j;
// 		// For points in depth
// 		for (int d=0;d<=depths[dnum];d++){
// 			const int p_temp = Util::ijk_to_p(i, j, sim.domain.znum - 1 - d, sim);
// 			grid.set_mpWidth(width, p_temp);
// 			grid.set_mpLength(length, p_temp);
// 			grid.set_mpDepth(depth, p_temp);
// 		}
// 	}

// }
