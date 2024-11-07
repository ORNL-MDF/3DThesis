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

#include <cmath>
#include <algorithm>
#include "DataStructs.h"
#include "Calc.h"
#include "Util.h"

#include <iostream>
#include <fstream>

void Calc::Integrate_Parallel(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	const int numThreads = sim.settings.thnum;
	static vector<Nodes> nodes_par;
	// If integrating in parallel
	if (numThreads>1) {
		// if the size is zero, integrate in parallel. Otherise, just take the last one (next time) from the previously calculated vector
		if (!nodes_par.size()) {
			nodes_par.resize(numThreads);
			#pragma omp parallel num_threads(numThreads)
			{
				Nodes th_nodes;
				#pragma omp for schedule(static)
				for (int i = 0; i < numThreads; i++) {
					Calc::Integrate_Serial(th_nodes, sim, t+i*sim.param.dt, isSol);
					nodes_par[numThreads - i - 1] = th_nodes;
				}
			}
		}
		nodes = nodes_par.back();
		nodes_par.pop_back();
	}
	else {
		Calc::Integrate_Serial(nodes, sim, t, isSol);
	}
	return;
}

void Calc::Integrate_Serial(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	if (sim.settings.compress) { 
		Calc::GaussCompressIntegrate(nodes, sim, t, isSol); 
		//If not solidifying, always choose the minimum of the two; otherwise, too expensive
		if (!isSol) {
			Nodes nodes_reg;
			Calc::GaussIntegrate(nodes_reg, sim, t, isSol);
			if (nodes_reg.size <= nodes.size) { nodes = nodes_reg; }
		}
	}
	else { 
		Calc::GaussIntegrate(nodes, sim, t, isSol);
	}

	if (sim.domain.use_BCs) { Calc::AddBCs(nodes, sim.domain); }

	return;
}

void Calc::GaussIntegrate(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	
	// This enables the starting search path segment to be quickly initialized each time. 
	static vector<int> start_seg(sim.paths.size(), 1);

	// Quadrature node locations for order 2, 4, 8, and 16
	static const double locs[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};

	// Quadrature weights for order 2, 4, 8, and 16
	static const double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};
	// For each beam and path
	for (int i = 0; i < sim.paths.size(); i++) {

		// Set beam, path, and seg_temp
		const Beam& beam = sim.beams[i];
		const vector<path_seg>& path = sim.paths[i];
		int seg_temp = start_seg[i];

		// Get beta for beam
		const double beta = pow(3.0 / PI, 1.5) * beam.q / (sim.material.rho * sim.material.cps);

		// Keep incrementing up if t is greater than the end of the path segment but also below the end time of the scan
		while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) { seg_temp++; }
		// Keep incrementing down if t is less than end of previous path segment
		while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) { seg_temp--; }
		
		// If not solidifying, then make this the start next time the function is run
		if (!isSol) {start_seg[i] = seg_temp;}

		// Get minimim time to integrate to
		const double t0 = Util::t0calc(t, beam, sim.material, sim.settings);
		
		// Set maximum starting step size to nonDimensional diffusion time
		double curStep_max_start = beam.nond_dt;
		
		// If solidifying, refine integration time
		if (isSol) { curStep_max_start *= beam.az; }

		// Set max step size to starting step size
		double curStep_max = curStep_max_start;

		// Set step size to max step
		double curStep_use = curStep_max;
		
		// Start quadrature order at 16
		int curOrder = 16;

		// Make 1st segment at the exact time
		// Used for instantaneous heat source additon to laplacian
		int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim); 
		current_beam.phix = (beam.ax * beam.ax + 0.0);
		current_beam.phiy = (beam.ay * beam.ay + 0.0);
		current_beam.phiz = (beam.az * beam.az + 0.0);
		current_beam.qmod *= beta;
		current_beam.dtau = 0.0;
		if (t <= path.back().seg_time && current_beam.qmod > 0.0) { Util::AddToNodes(nodes, current_beam); }

		bool tflag = true;
		double t2 = t;
		double t1 = t;
		double tpp = 0.0;

		if (t > path.back().seg_time) {
			tpp += t - path.back().seg_time;
			t2 = path.back().seg_time;
		}

		while (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}

		while (tflag) {
			bool switchSeg = false;

			// Get maximum integration step for the segment
			const double ref_time = Util::GetRefTime(tpp, seg_temp, path, beam);

			// Sets step to the minimum of the ref time and max allowable time
			curStep_use = (ref_time < curStep_max) ? ref_time : curStep_max;

			// Set ending quadradure time based on step used
			t1 = t2 - curStep_use;

			// Time the next segment ends
			const double next_time = path[seg_temp - 1].seg_time;

			//If we are at the end of a segment, hit the end of it and set the program to jump to the next segment next time
			if (t1 <= next_time) {
				t1 = next_time;
				// If next time is greater than t0, switch segments and keep going
				if (next_time > t0) { switchSeg = true; }
				// Otherwise, it should end
				else { tflag = false; }
			}

			//Add Quadrature Points
			double tau, ct;
			for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
				double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));
				tau = t - tp;
				ct = 12.0 * sim.material.a * tau;

				current_beam = Util::GetBeamLoc(tp, seg_temp, path, sim);
				current_beam.phix = (beam.ax * beam.ax + ct);
				current_beam.phiy = (beam.ay * beam.ay + ct);
				current_beam.phiz = (beam.az * beam.az + ct);
				current_beam.qmod *= beta;
				current_beam.dtau = 0.5 * (t2 - t1) * weights[a];

				if (current_beam.qmod > 0.0 && current_beam.dtau > 0.0) { Util::AddToNodes(nodes, current_beam); }
			}

			//If we are switching segments, increment the start segment down
			if (switchSeg) {seg_temp--;}
			
			// Increment the total time passed
			tpp += (t2 - t1);

			// Set start quadruatre time to current end quadrature time
			t2 = t1;

			//Increase Maximum Step and Decrease Gauss order if it is okay to do so
			if (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}
		}
	}
	return;
}

void Calc::GaussCompressIntegrate(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	// This enables the starting search path segment to be quickly initialized each time. 
	static vector<int> start_seg(sim.paths.size(), 1);

	// Quadrature node locations for order 2, 4, 8, and 16
	static const double locs[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};

	// Quadrature weights for order 2, 4, 8, and 16
	static const double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};

	// For each beam and path
	for (int i = 0; i < sim.paths.size(); i++) {

		// Set beam, path, and seg_temp
		const Beam& beam = sim.beams[i];
		const vector<path_seg>& path = sim.paths[i];
		int seg_temp = start_seg[i];

		// Get beta for beam
		const double beta = pow(3.0 / PI, 1.5) * beam.q / (sim.material.rho * sim.material.cps);

		// Keep incrementing up if t is greater than the end of the path segment but also below the end time of the scan
		while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) { seg_temp++; }
		// Keep incrementing down if t is less than end of previous path segment
		while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) { seg_temp--; }

		// If not solidifying, then make this the start next time the function is run
		if (!isSol) { start_seg[i] = seg_temp; }

		// Get minimim time to integrate to
		const double t0 = Util::t0calc(t, beam, sim.material, sim.settings);

		// Set maximum starting step size to nonDimensional diffusion time
		double curStep_max_start = beam.nond_dt;

		// If solidifying, refine integration time
		if (isSol) { curStep_max_start *= beam.az; }

		// Set max step size to starting step size
		double curStep_max = curStep_max_start;

		// Set step size to max step
		double curStep_use = curStep_max;

		// Start quadrature order at 16
		int curOrder = 16;

		// Make 1st segment at the exact time
		// Used for instantaneous heat source additon to laplacian
		int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim);
		current_beam.phix = (beam.ax * beam.ax + 0.0);
		current_beam.phiy = (beam.ay * beam.ay + 0.0);
		current_beam.phiz = (beam.az * beam.az + 0.0);
		current_beam.qmod *= beta;
		current_beam.dtau = 0.0;
		if (t <= path.back().seg_time && current_beam.qmod>0.0) { Util::AddToNodes(nodes, current_beam); }

		bool tflag = true;
		double t2 = t;
		double t1 = t;
		double tpp = 0.0;

		if (t > path.back().seg_time) {
			tpp += t - path.back().seg_time;
			t2 = path.back().seg_time;
		}

		while (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}

		while (tflag) {

			//Compression variables
			double r2, dist2, xp, yp, xs, ys, dx, dy, ts, dt;
			double sum_t = 0, sum_qmodt = 0, sum_qmodtx = 0, sum_qmodty = 0; //sums of t, qmod*t, qmod*t*x,... to find centers
			int seg_temp_2;
			int num_comb_segs = 0; //number of segments to be combined

			//If the time is less than t0, break the whole thing
			int quit = 0;
			//If outside r, then keep going down until back in r
			//NOTE:Only looks at end of scan path...good for points but not lines
			while (true) {
				xs = path[seg_temp].sx;
				ys = path[seg_temp].sy;
				ts = path[seg_temp].seg_time;
				if (ts <= t0) { quit = 1; break; }
				// If in R, quit loop
				if (Util::InRMax(xs, ys, sim.domain, sim.settings)) { break; }
				// If outside R, add cumulative time, set time to be end of segment
				else { tpp += t2 - ts; t2 = ts; }
				// if (!Util::InRMax(xs, ys, sim)) { seg_temp--; }
				// else { tpp += t2 - ts; spp = tpp / sim.util.nond_dt; t2 = ts; break; }
			}
			if (quit) { break; }

			while (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}

			double ref_time = Util::GetRefTime(tpp, seg_temp, path, beam);
			if (ref_time < curStep_max) { curStep_use = ref_time; }
			else { curStep_use = curStep_max; }

			int_seg current_beam_t2 = Util::GetBeamLoc(t2, seg_temp, path, sim);
			xp = current_beam_t2.xb;
			yp = current_beam_t2.yb;

			// Diffusion Distance Squared (Distance it diffused by *some amount*, squared)
			r2 = log(2.0) / 8.0 * (beam.ax * beam.ax) * (12.0 * (t - t2) * sim.material.a / (beam.ax * beam.ax) + 1.0);

			seg_temp_2 = seg_temp;

			t1 = path[seg_temp - 1].seg_time;

			//If we won't be switching segments, do normal integration
			if (t1 < t2 - curStep_use) {
				num_comb_segs = 0;
				t1 = t2 - curStep_use;
			}
			else {
				bool cflag = true;
				while (cflag) { //If we will be switching segments, do compressed integration
								//If the segment start time is zero, end the loop
					if (path[seg_temp_2 - 1].seg_time <= t0) { tflag = 0; break; }

					xs = path[seg_temp_2 - 1].sx;
					ys = path[seg_temp_2 - 1].sy;
					dist2 = (xs - xp) * (xs - xp) + (ys - yp) * (ys - yp);
					ts = path[seg_temp_2 - 1].seg_time;

					// If the next segment is outside the calculation domain, break the loop
					// NOTE:Only looks at end of scan path...good for points but not lines
					if (!Util::InRMax(xp, yp, sim.domain, sim.settings)) { break; }

					// IF 
					//	the distance between endpoints, or points, is bigger than the diffusion distance
					// OR
					//	the qmod's aren't equal (discontinuity) AND current order or time difference is *as defined* (to minimize error)
					// THEN
					//	stop combining segments and finalize
					// ELSE
					//  If the distance of the next segment is less than the diffusion distance, set end time to next segment, average them, and keep going
					if ((dist2 > r2) || (path[seg_temp_2].sqmod != path[seg_temp_2 - 1].sqmod && (curOrder > 2 || (t2 - ts) > (curStep_max / 64.0)))) {
						//If point, do averaging calculations and set new end time to the next segment
						if (path[seg_temp_2].smode) {
							sum_qmodtx += path[seg_temp_2].sx * path[seg_temp_2].sqmod * (t1 - ts);
							sum_qmodty += path[seg_temp_2].sy * path[seg_temp_2].sqmod * (t1 - ts);
							sum_qmodt += path[seg_temp_2].sqmod * (t1 - ts);
							sum_t += (t1 - ts);
							t1 = ts;
						}
						//If line, find the time when the dist=r. If this time is less than the next segment start time, set end time to next segment; else, set start time to when dist=r
						else {
							dx = path[seg_temp_2].sx - path[seg_temp_2 - 1].sx;
							dy = path[seg_temp_2].sy - path[seg_temp_2 - 1].sy;
							dt = path[seg_temp_2].seg_time - ts;

							double t_int = ts + dt * (sqrt(((xp - xs) * (xp - xs) + (yp - ys) * (yp - ys)) / (dx * dx + dy * dy)) - sqrt(r2 / (dx * dx + dy * dy)));
							if (t_int < t2 - curStep_max) { t_int = t2 - curStep_max; }
							sum_qmodtx += (xs + dx * ((t1 + t_int) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * (t1 - t_int);
							sum_qmodty += (ys + dy * ((t1 + t_int) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * (t1 - t_int);
							sum_qmodt += path[seg_temp_2].sqmod * (t1 - t_int);
							sum_t += (t1 - t_int);
							t1 = t_int;
						}
						cflag = false;
					}
					else {
						//If the next time is less than the minimum allowed time, set appropriate times to the minimum allowed time
						if (ts < t2 - curStep_max) {
							cflag = false;
							if (path[seg_temp_2].smode) { t1 = t2 - curStep_max; }
							else { t1 = t2 - curStep_max; }
						}
						//Otherwise combine the segments and continue
						else {
							t1 = ts;
							num_comb_segs++;
						}

						//If point, 
						if (path[seg_temp_2].smode) {
							dt = (path[seg_temp_2].seg_time - t1);

							sum_qmodtx += path[seg_temp_2].sx * path[seg_temp_2].sqmod * dt;
							sum_qmodty += path[seg_temp_2].sy * path[seg_temp_2].sqmod * dt;
							sum_qmodt += path[seg_temp_2].sqmod * dt;
							sum_t += dt;
						}
						else {
							dx = path[seg_temp_2].sx - path[seg_temp_2 - 1].sx;
							dy = path[seg_temp_2].sy - path[seg_temp_2 - 1].sy;
							dt = (path[seg_temp_2].seg_time - t1);

							sum_qmodtx += (xs + dx * ((path[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * dt;
							sum_qmodty += (ys + dy * ((path[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * dt;
							sum_qmodt += path[seg_temp_2].sqmod * dt;
							sum_t += dt;
						}
					}
					if (t1 == path[seg_temp_2 - 1].seg_time) {
						xs = path[seg_temp_2 - 1].sx;
						ys = path[seg_temp_2 - 1].sy;
						seg_temp_2--;
						if (!Util::InRMax(xs, ys, sim.domain, sim.settings)) { break; }
					}
				}

				if (!num_comb_segs) { //If we didn't combine any segments anyways, then just do normal integration
					t1 = path[seg_temp - 1].seg_time;
					seg_temp_2 = seg_temp - 1;
				}
			}

			//Add Quadrature Points using the same average for x, y, and qmod (CAN BE IMPROVED)
			double tau, ct;
			if (num_comb_segs) {
				int_seg current_beam;
				if (sum_qmodt > 0.0) {
					current_beam.xb = sum_qmodtx / sum_qmodt;
					current_beam.yb = sum_qmodty / sum_qmodt;
					current_beam.qmod = sum_qmodt / sum_t;		
					for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
						double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));
						tau = t - tp;
						ct = 12.0 * sim.material.a * tau;

						current_beam.phix = (beam.ax * beam.ax + ct);
						current_beam.phiy = (beam.ay * beam.ay + ct);
						current_beam.phiz = (beam.az * beam.az + ct);
						current_beam.qmod *= beta;
						current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
						if ((current_beam.qmod > 0.0) && (current_beam.dtau > 0.0)) { Util::AddToNodes(nodes, current_beam); }
					}
				}

			}
			//Add Quadrature Points
			else {
				for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
					double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));	
					tau = t - tp;
					ct = 12.0 * sim.material.a * tau;

					int_seg current_beam = Util::GetBeamLoc(tp, seg_temp, path, sim);
					current_beam.phix = (beam.ax * beam.ax + ct);
					current_beam.phiy = (beam.ay * beam.ay + ct);
					current_beam.phiz = (beam.az * beam.az + ct);
					current_beam.qmod *= beta;
					current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
					if ((current_beam.qmod > 0.0) && (current_beam.dtau > 0.0)) { Util::AddToNodes(nodes, current_beam); }
				}
			}

			tpp += t2 - t1;
			t2 = t1;

			//Increase Maximum Step and Decrease Gauss order if it is "safe" to
			if (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}

			seg_temp = seg_temp_2;
			if (t1 <= t0) { tflag = 0; }
		}
	}
	return;
}

void Calc::AddBCs(Nodes& nodes, const Domain& domain) {

	vector<vector<int>> allCoords;
	vector<vector<int>> newCoords;
	vector<int> coords = { 0,0,0 };
	newCoords.push_back(coords);
	allCoords.push_back(coords);

	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	xmin = domain.BC_xmin;
	xmax = domain.BC_xmax;
	ymin = domain.BC_ymin;
	ymax = domain.BC_ymax;
	zmin = domain.BC_zmin;
	zmax = domain.zmax;

	vector<vector<int>> checkCoords = newCoords;
	newCoords.clear();
	// Strength is how many reflections to use
	int iter = 0; int strength = 0; 
	if (domain.BC_reflections != INT_MAX) { strength = domain.BC_reflections; }
	while (checkCoords.size() > 0 && iter < strength) {
		for (int i = 0; i < checkCoords.size(); i++) {
			if (xmin != DBL_MAX) {
				coords = checkCoords[i];
				coords[0] = (-1 - coords[0]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {;
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (xmax != DBL_MAX) {
				coords = checkCoords[i];
				coords[0] = (1 - coords[0]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (ymin != DBL_MAX) {
				coords = checkCoords[i];
				coords[1] = (-1 - coords[1]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (ymax != DBL_MAX) {
				coords = checkCoords[i];
				coords[1] = (1 - coords[1]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (zmin != DBL_MAX) {
				coords = checkCoords[i];
				coords[2] = (-1 - coords[2]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (zmax != DBL_MAX) {
				coords = checkCoords[i];
				coords[2] = (1 - coords[2]);
				
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end() && (coords[2]!=1)) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
		}
		checkCoords = newCoords;
		newCoords.clear();
		iter += 1;
	}

	int org_size = nodes.size;
	Nodes nodes_org = nodes;
	Util::ClearNodes(nodes);

	// Now just convert between refCoords and realCoords
	for (int i = 0; i < allCoords.size(); i++) {
		int xRef = allCoords[i][0];
		int yRef = allCoords[i][1];
		int zRef = allCoords[i][2];	

		int nX = std::abs(xRef % 2);
		int nY = std::abs(yRef % 2);
		int nZ = std::abs(zRef % 2);	

		Nodes nodes2 = nodes_org;
		for (int i = 0; i < org_size; i++) { 
			nodes2.xb[i] = (xRef + nX) * xmax - (xRef - nX) * xmin + (1 - 2 * nX) * nodes2.xb[i];
			nodes2.yb[i] = (yRef + nY) * ymax - (yRef - nY) * ymin + (1 - 2 * nY) * nodes2.yb[i];
			nodes2.zb[i] = (zRef + nZ) * zmax - (zRef - nZ) * zmin + (1 - 2 * nZ) * nodes2.zb[i];
		}
		Util::CombineNodes(nodes, nodes2);	
	}

	return;
}