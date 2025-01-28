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
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "impl/Structs/DataStructs.hpp"
#include "impl/Calc/Calc.hpp"
#include "impl/Calc/Util.hpp"

#ifdef __GNUC__ // GCC or Clang
#include <cstdint>
inline int my_bit_width(uint64_t x) {
    if (x == 0) return 0;
    return 64 - __builtin_clzll(x); // Count leading zeros
}
#elif defined(_MSC_VER) // Microsoft Visual C++
#include <intrin.h>
#include <cstdint>
inline int my_bit_width(uint64_t x) {
    if (x == 0) return 0;
    unsigned long index;
    _BitScanReverse64(&index, x); // Find index of most significant set bit
    return index + 1;
}
#else
#error "Compiler not supported. Please implement my_bit_width for your compiler."
#endif

namespace Thesis::impl
{

	void Calc::Integrate_Parallel(Nodes& nodes, vector<Nodes>& nodes_par, vector<int>& start_seg, const Simdat& sim, const double t, const bool isSol) {
		const int numThreads = sim.settings.thnum;
		// If it hasn't been initialized yet, initilaize it
		if (!start_seg.size()){
			// This enables the starting search path segment to be quickly initialized each time. 
			start_seg = vector<int>(sim.paths.size(), 1);
		}
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
						Calc::Integrate_Serial(th_nodes, start_seg, sim, t+i*sim.param.dt, isSol);
						nodes_par[numThreads - i - 1] = th_nodes;
					}
				}
			}
			nodes = nodes_par.back();
			nodes_par.pop_back();
		}
		else {
			Calc::Integrate_Serial(nodes, start_seg, sim, t, isSol);
		}
		return;
	}

	void Calc::Integrate_Serial(Nodes& nodes, vector<int>& start_seg, const Simdat& sim, const double t, const bool isSol) {
		// If it hasn't been initialized yet, initilaize it
		if (!start_seg.size()){
			// This enables the starting search path segment to be quickly initialized each time. 
			start_seg = vector<int>(sim.paths.size(), 1);
		}
		
		if (sim.settings.compress) { 
			Calc::GaussCompressIntegrate(nodes, start_seg, sim, t, isSol); 
			//If not solidifying, always choose the minimum of the two; otherwise, too expensive
			if (!isSol) {
				Nodes nodes_reg;
				Calc::GaussIntegrate(nodes_reg, start_seg, sim, t, isSol);
				if (nodes_reg.size <= nodes.size) { nodes = nodes_reg; }
			}
		}
		else { 
			Calc::GaussIntegrate(nodes, start_seg, sim, t, isSol);
		}

		if (sim.domain.use_BCs) { Calc::AddBCs(nodes, sim.domain); }

		return;
	}

	// Add these helper functions in Util namespace
	namespace impl {

		// Fast log2 floor using bit scanning
		inline int log2_floor(double x) {
			return my_bit_width(static_cast<int>(x)) - 1;
		}

		// Fast log2 ceil using bit scanning
		inline int log2_ceil(double x) {
			return my_bit_width(static_cast<int>(x) - 1);
		}

		// Analytic determination for spot mode steps
		inline int spot_steps(double s_start, double s_end) {
			const int delta = static_cast<int>(std::ceil(s_end - s_start));
			const int prevPow = log2_floor(s_start+1.0);
			const int prevMinusOne = (1 << prevPow) - 1;
			return (log2_floor(prevMinusOne + delta)-prevPow) + 1;
		}

		// Upper bound for line mode steps
		inline int line_steps(double s_start, double s_end, double V) {
			//constexpr double lineConst_z = 12.0 * std::sqrt(std::log(std::sqrt(2)));
			constexpr double lineConst_z = 7.064460135092848;
			const double A = lineConst_z / V;
			const double z_start = 12.0 * s_start;
			const double z_end = 12.0 * s_end;

			const double t0 = std::sqrt(z_start + 1.0);
			const double T = std::sqrt(z_end + 1.0);
			const double denominator = std::sqrt(t0 * t0 + A * t0) - t0;
			
			return static_cast<int>(std::ceil((T - t0) / denominator));
		}
	}

	void Calc::GaussIntegrate_Parallel(Nodes& nodes, vector<int>& start_seg, const Simdat& sim, const double t, const bool isSol) {
				
		static constexpr double lineConst_s = 0.5887050112577373;
		static constexpr double piConst = 0.933162059717596;
		static constexpr int ORDER_OFFSETS[4] = {0, 2, 6, 14};
		static constexpr int ORDER_COUNTS[4] = {2, 4, 8, 16};

		// Quadrature node locations for order 2, 4, 8, and 16
		static constexpr double locs[30] = {
			-0.57735027,  0.57735027,
			-0.86113631, -0.33998104,  0.33998104,  0.86113631,
			-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
			-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
		};

		// Quadrature weights for order 2, 4, 8, and 16
		static constexpr double weights[30] = {
			1.0, 1.0,
			0.34785485, 0.65214515, 0.65214515, 0.34785485,
			0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
			0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
		};

		for (int i = 0; i < sim.paths.size(); i++) {
			const Beam& beam = sim.beams[i];
			const vector<path_seg>& path = sim.paths[i];
			int seg_temp = start_seg[i];

			// Set material constant
			const double beta = piConst * beam.q / (sim.material.rho * sim.material.cps);

			while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) {seg_temp++;}
			while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) {seg_temp--;}
			
			if (!isSol && (omp_get_thread_num() == 0)){start_seg[i] = seg_temp;}

			// Reserve size for nodes
			nodes.xb.reserve(1000);
			nodes.yb.reserve(1000);
			nodes.zb.reserve(1000);
			nodes.phix.reserve(1000);
			nodes.phiy.reserve(1000);
			nodes.phiz.reserve(1000);
			nodes.dtau.reserve(1000);
			nodes.expmod.reserve(1000);

			// Add current beam only if solidifying 
			if (isSol){
				int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim);
				current_beam.phix = beam.ax * beam.ax;
				current_beam.phiy = beam.ay * beam.ay;
				current_beam.phiz = beam.az * beam.az;
				current_beam.qmod *= beta;
				current_beam.dtau = 0.0;
				if (t <= path.back().seg_time && current_beam.qmod > 0.0) Util::AddToNodes(nodes, current_beam);
			}

			// Make vectors of points and lines
			const int seg_max = seg_temp;
			vector<int> spot_segs; spot_segs.reserve(seg_max);
			vector<int> line_segs; line_segs.reserve(seg_max);
			for (int seg = 1; seg <= seg_max; seg++) {
				if (path[seg].sqmod == 0.0){continue;}
				if (path[seg].smode == 1){spot_segs.push_back(seg);}
				else{line_segs.push_back(seg);}
			}

			// Set consts and reserve sizes
			const size_t spot_num = spot_segs.size();
			const size_t line_num = line_segs.size();

			// Do just the spots
			#pragma omp parallel num_threads(sim.settings.thnum)
			{
				// Thread local variables
				int th_size = 0;
				vector<double> th_xb; th_xb.reserve(128);
				vector<double> th_yb; th_yb.reserve(128);
				vector<double> th_zb; th_zb.reserve(128);
				vector<double> th_phix; th_phix.reserve(128);
				vector<double> th_phiy; th_phiy.reserve(128);
				vector<double> th_phiz; th_phiz.reserve(128);
				vector<double> th_dtau; th_dtau.reserve(128);
				vector<double> th_expmod; th_expmod.reserve(128);

				#pragma omp for schedule(static)
				for (size_t spot_seg=0; spot_seg<spot_num; spot_seg++){
					// Get actual seg
					const int& seg = spot_segs[spot_seg];
					
					// Get position information
					const double& xb = path[seg].sx;
					const double& yb = path[seg].sy;
					const double& zb = path[seg].sz;
					const double qmod = path[seg].sqmod*beta;
					if (!Util::InRMax(xb,yb,sim.domain,sim.settings)){continue;}

					// Set path information
					const double nond_dt = beam.nond_dt;
					double tp_start = max(t - path[seg].seg_time, 0.0);
					double tp_end = t - path[seg - 1].seg_time;
					double sp_start = tp_start / nond_dt;
					double sp_end = tp_end / nond_dt;

					// Calculate maximum number of steps
					const int steps = impl::spot_steps(sp_start, sp_end);
					// const int steps = (sp_start-sp_end < sp_start) ? 1 : impl::spot_steps(sp_start, sp_end);

					// Start stepping
					double s_current = sp_start;
					for (int step = 0; step < steps; ++step) {
						
						// Order and step calculation
						const int flat = impl::log2_floor(s_current + 1.0);
						const int k = std::max(3-flat,0);
						const int step_size = (1 << flat);
						const int order = ORDER_COUNTS[k];
						const int offset = ORDER_OFFSETS[k];
						const double s_next = std::min(s_current + step_size, sp_end);

						// Precalculate quadrature stuff
						const double s_diff = s_next - s_current;
						const double s_sum = s_next + s_current;
						
						// Quadrature
						for (int q = 0; q < order; ++q) {

							// Quadrature
							const double loc = locs[offset + q];
							const double weight = weights[offset + q];
							const double ct = 6.0 * (s_diff * loc + s_sum) * sim.material.a * nond_dt;
							const double phix = 1.0/(beam.ax * beam.ax + ct);	// Precompute division
							const double phiy = 1.0/(beam.ay * beam.ay + ct);	// Precompute division
							const double phiz = 1.0/(beam.az * beam.az + ct);	// Precompute division
							const double dtau = 0.5 * s_diff * nond_dt * weight;

							if (dtau > 0.0) {
								// Add information to nodes
								th_size++;
								th_xb.push_back(xb);
								th_yb.push_back(yb);
								th_zb.push_back(zb);
								th_phix.push_back(phix);	
								th_phiy.push_back(phiy);	
								th_phiz.push_back(phiz);	
								th_dtau.push_back(dtau);
								th_expmod.push_back(0.5*log(qmod*qmod*phix*phiy*phiz));
							}
						}

						// Set s for next loop
						s_current = s_next;
					}
				}

				#pragma omp for schedule(static)
				// Do just the lines
				for (size_t line_seg=0; line_seg<line_num; line_seg++){
					const int& seg = line_segs[line_seg];

					// Set object and variables which don't change based on node				
					const double px =  path[seg - 1].sx;
					const double py =  path[seg - 1].sy;
					const double pz =  path[seg - 1].sz;
					const double pt =  path[seg - 1].seg_time;
					const double dx = path[seg].sx - path[seg - 1].sx;
					const double dy = path[seg].sy - path[seg - 1].sy;
					const double dz = path[seg].sz - path[seg - 1].sz;
					const double dt_cur = path[seg].seg_time - path[seg - 1].seg_time;
					const double qmod = path[seg].sqmod*beta;

					// Set path information
					const double nond_dt = beam.nond_dt;
					double tp_start = max(t - path[seg].seg_time, 0.0);
					double tp_end = t - path[seg - 1].seg_time;
					double sp_start = tp_start / nond_dt;
					double sp_end = tp_end / nond_dt;
					const double V = path[seg].sparam * (beam.ax / sim.material.a);

					// Calculate maximum number of steps
					const int steps = impl::line_steps(sp_start, sp_end, V);

					// Start stepping
					double s_current = sp_start;
					for (int step = 0; step < steps; ++step) {
						if (s_current == sp_end) {continue;}

						// Order and step calculation
						const int flat = impl::log2_floor(s_current + 1.0);
						const int k = std::max(3-flat,0);
						const double step_size = lineConst_s / V * sqrt(12.0 * s_current + 1.0);
						const int order = ORDER_COUNTS[k];
						const int offset = ORDER_OFFSETS[k];
						const double s_next = std::min(s_current + step_size, sp_end);

						// Precalculate quadrature stuff
						const double s_diff = s_next - s_current;
						const double s_sum = s_next + s_current;

						// Quadrature
						for (int q = 0; q < order; ++q) {					

							// Quadrature
							const double loc = locs[offset + q];
							const double weight = weights[offset + q];
							const double tau = 0.5 * (s_diff * loc + s_sum) * nond_dt;
							const double ct = 12.0 * sim.material.a * tau;
							const double phix = 1.0/(beam.ax * beam.ax + ct);	// Precompute division
							const double phiy = 1.0/(beam.ay * beam.ay + ct);	// Precompute division
							const double phiz = 1.0/(beam.az * beam.az + ct);	// Precompute division
							const double dtau = 0.5 * s_diff * nond_dt * weight;

							// Location
							const double tcur = (t - tau) - path[seg - 1].seg_time;
							const double xb = path[seg - 1].sx + (tcur / dt_cur)*dx;
							const double yb = path[seg - 1].sy + (tcur / dt_cur)*dy;
							const double zb = path[seg - 1].sz + (tcur / dt_cur)*dz;

							if (Util::InRMax(xb,yb,sim.domain,sim.settings) && dtau > 0.0) {	
								// Add information to nodes
								th_size++;
								th_xb.push_back(xb);
								th_yb.push_back(yb);
								th_zb.push_back(zb);
								th_phix.push_back(phix);	
								th_phiy.push_back(phiy);	
								th_phiz.push_back(phiz);	
								th_dtau.push_back(dtau);
								th_expmod.push_back(0.5*log(qmod*qmod*phix*phiy*phiz));
							}
						}

						// Set s for next loop
						s_current = s_next;
					}
				}

				// Now add to nodes
				#pragma omp critical
				{
					nodes.size += th_size;
					nodes.xb.insert(nodes.xb.end(), th_xb.begin(), th_xb.end());
					nodes.yb.insert(nodes.yb.end(), th_yb.begin(), th_yb.end());
					nodes.zb.insert(nodes.zb.end(), th_zb.begin(), th_zb.end());
					nodes.phix.insert(nodes.phix.end(), th_phix.begin(), th_phix.end());
					nodes.phiy.insert(nodes.phiy.end(), th_phiy.begin(), th_phiy.end());
					nodes.phiz.insert(nodes.phiz.end(), th_phiz.begin(), th_phiz.end());
					nodes.dtau.insert(nodes.dtau.end(), th_dtau.begin(), th_dtau.end());
					nodes.expmod.insert(nodes.expmod.end(), th_expmod.begin(), th_expmod.end());
				}
			}
		}
	}

	void Calc::GaussIntegrate(Nodes& nodes, vector<int>& start_seg, const Simdat& sim, const double t, const bool isSol) {
				
		static constexpr double lineConst_s = 0.5887050112577373;
		static constexpr double piConst = 0.933162059717596;
		static constexpr int ORDER_OFFSETS[4] = {0, 2, 6, 14};
		static constexpr int ORDER_COUNTS[4] = {2, 4, 8, 16};

		// Quadrature node locations for order 2, 4, 8, and 16
		static constexpr double locs[30] = {
			-0.57735027,  0.57735027,
			-0.86113631, -0.33998104,  0.33998104,  0.86113631,
			-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
			-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
		};

		// Quadrature weights for order 2, 4, 8, and 16
		static constexpr double weights[30] = {
			1.0, 1.0,
			0.34785485, 0.65214515, 0.65214515, 0.34785485,
			0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
			0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
		};

		for (int i = 0; i < sim.paths.size(); i++) {
			const Beam& beam = sim.beams[i];
			const vector<path_seg>& path = sim.paths[i];
			int seg_temp = start_seg[i];

			// Set material constant
			const double beta = piConst * beam.q / (sim.material.rho * sim.material.cps);

			while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) {seg_temp++;}
			while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) {seg_temp--;}
			
			if (!isSol && (omp_get_thread_num() == 0)){start_seg[i] = seg_temp;}

			// Reserve size for nodes
			nodes.xb.reserve(1000);
			nodes.yb.reserve(1000);
			nodes.zb.reserve(1000);
			nodes.phix.reserve(1000);
			nodes.phiy.reserve(1000);
			nodes.phiz.reserve(1000);
			nodes.dtau.reserve(1000);
			nodes.expmod.reserve(1000);

			// Add current beam only if solidifying 
			if (isSol){
				int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim);
				current_beam.phix = beam.ax * beam.ax;
				current_beam.phiy = beam.ay * beam.ay;
				current_beam.phiz = beam.az * beam.az;
				current_beam.qmod *= beta;
				current_beam.dtau = 0.0;
				if (t <= path.back().seg_time && current_beam.qmod > 0.0) Util::AddToNodes(nodes, current_beam);
			}

			// Make vectors of points and lines
			const int seg_max = seg_temp;
			vector<int> spot_segs; spot_segs.reserve(seg_max);
			vector<int> line_segs; line_segs.reserve(seg_max);
			for (int seg = 1; seg <= seg_max; seg++) {
				if (path[seg].sqmod == 0.0){continue;}
				if (path[seg].smode == 1){spot_segs.push_back(seg);}
				else{line_segs.push_back(seg);}
			}

			// Set consts and reserve sizes
			const size_t spot_num = spot_segs.size();
			const size_t line_num = line_segs.size();

			// Do just the spots
			for (size_t spot_seg=0; spot_seg<spot_num; spot_seg++){
				
				// Get actual seg
				const int& seg = spot_segs[spot_seg];
				
				// Get position information
				const double& xb = path[seg].sx;
				const double& yb = path[seg].sy;
				const double& zb = path[seg].sz;
				const double qmod = path[seg].sqmod*beta;
				if (!Util::InRMax(xb,yb,sim.domain,sim.settings)){continue;}

				// Set path information
				const double nond_dt = beam.nond_dt;
				double tp_start = max(t - path[seg].seg_time, 0.0);
				double tp_end = t - path[seg - 1].seg_time;
				double sp_start = tp_start / nond_dt;
				double sp_end = tp_end / nond_dt;

				// Calculate maximum number of steps
				const int steps = impl::spot_steps(sp_start, sp_end);
				// const int steps = (sp_start-sp_end < sp_start) ? 1 : impl::spot_steps(sp_start, sp_end);

				// Start stepping
				double s_current = sp_start;
				for (int step = 0; step < steps; ++step) {
					
					// Order and step calculation
					const int flat = impl::log2_floor(s_current + 1.0);
					const int k = std::max(3-flat,0);
					const int step_size = (1 << flat);
					const int order = ORDER_COUNTS[k];
					const int offset = ORDER_OFFSETS[k];
					const double s_next = std::min(s_current + step_size, sp_end);

					// Precalculate quadrature stuff
					const double s_diff = s_next - s_current;
					const double s_sum = s_next + s_current;
					
					// Quadrature
					for (int q = 0; q < order; ++q) {

						// Quadrature
						const double loc = locs[offset + q];
						const double weight = weights[offset + q];
						const double ct = 6.0 * (s_diff * loc + s_sum) * sim.material.a * nond_dt;
						const double phix = 1.0/(beam.ax * beam.ax + ct);	// Precompute division
						const double phiy = 1.0/(beam.ay * beam.ay + ct);	// Precompute division
						const double phiz = 1.0/(beam.az * beam.az + ct);	// Precompute division
						const double dtau = 0.5 * s_diff * nond_dt * weight;

						if (dtau > 0.0) {
							// Add information to nodes
							nodes.size++;
							nodes.xb.push_back(xb);
							nodes.yb.push_back(yb);
							nodes.zb.push_back(zb);
							nodes.phix.push_back(phix);	
							nodes.phiy.push_back(phiy);	
							nodes.phiz.push_back(phiz);	
							nodes.dtau.push_back(dtau);
							nodes.expmod.push_back(0.5*log(qmod*qmod*phix*phiy*phiz));
						}
					}

					// Set s for next loop
					s_current = s_next;
				}
			}

			// Do just the lines
			for (size_t line_seg=0; line_seg<line_num; line_seg++){
				const int& seg = line_segs[line_seg];

				// Set object and variables which don't change based on node				
				const double px =  path[seg - 1].sx;
				const double py =  path[seg - 1].sy;
				const double pz =  path[seg - 1].sz;
				const double pt =  path[seg - 1].seg_time;
				const double dx = path[seg].sx - path[seg - 1].sx;
				const double dy = path[seg].sy - path[seg - 1].sy;
				const double dz = path[seg].sz - path[seg - 1].sz;
				const double dt_cur = path[seg].seg_time - path[seg - 1].seg_time;
				const double qmod = path[seg].sqmod*beta;

				// Set path information
				const double nond_dt = beam.nond_dt;
				double tp_start = max(t - path[seg].seg_time, 0.0);
				double tp_end = t - path[seg - 1].seg_time;
				double sp_start = tp_start / nond_dt;
				double sp_end = tp_end / nond_dt;
				const double V = path[seg].sparam * (beam.ax / sim.material.a);

				// Calculate maximum number of steps
				const int steps = impl::line_steps(sp_start, sp_end, V);

				// Start stepping
				double s_current = sp_start;
				for (int step = 0; step < steps; ++step) {
					if (s_current == sp_end) {continue;}

					// Order and step calculation
					const int flat = impl::log2_floor(s_current + 1.0);
					const int k = std::max(3-flat,0);
					const double step_size = lineConst_s / V * sqrt(12.0 * s_current + 1.0);
					const int order = ORDER_COUNTS[k];
					const int offset = ORDER_OFFSETS[k];
					const double s_next = std::min(s_current + step_size, sp_end);

					// Precalculate quadrature stuff
					const double s_diff = s_next - s_current;
					const double s_sum = s_next + s_current;

					// Quadrature
					for (int q = 0; q < order; ++q) {					

						// Quadrature
						const double loc = locs[offset + q];
						const double weight = weights[offset + q];
						const double tau = 0.5 * (s_diff * loc + s_sum) * nond_dt;
						const double ct = 12.0 * sim.material.a * tau;
						const double phix = 1.0/(beam.ax * beam.ax + ct);	// Precompute division
						const double phiy = 1.0/(beam.ay * beam.ay + ct);	// Precompute division
						const double phiz = 1.0/(beam.az * beam.az + ct);	// Precompute division
						const double dtau = 0.5 * s_diff * nond_dt * weight;

						// Location
						const double tcur = (t - tau) - path[seg - 1].seg_time;
						const double xb = path[seg - 1].sx + (tcur / dt_cur)*dx;
						const double yb = path[seg - 1].sy + (tcur / dt_cur)*dy;
						const double zb = path[seg - 1].sz + (tcur / dt_cur)*dz;

						if (Util::InRMax(xb,yb,sim.domain,sim.settings) && dtau > 0.0) {
							
							// Add information to nodes
							nodes.size++;
							nodes.xb.push_back(xb);
							nodes.yb.push_back(yb);
							nodes.zb.push_back(zb);
							nodes.phix.push_back(phix);	
							nodes.phiy.push_back(phiy);	
							nodes.phiz.push_back(phiz);	
							nodes.dtau.push_back(dtau);
							nodes.expmod.push_back(0.5*log(qmod*qmod*phix*phiy*phiz));
						}
					}

					// Set s for next loop
					s_current = s_next;
				}
			}
		}
	}

	void Calc::GaussCompressIntegrate(Nodes& nodes, vector<int>& start_seg, const Simdat& sim, const double t, const bool isSol) {
		
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
}