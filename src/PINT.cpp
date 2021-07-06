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

#include <vector>
#include <omp.h>
#include <cmath>

#include "PINT.h"
#include "Util.h"
#include "DataStructs.h"
#include "Point.h"

void PINT::GodCheck(std::deque<Point>& ptv, vector<int>& god_pts, vector<omp_lock_t>& lock, Simdat& sim, int point_num) {
	if (god_pts[point_num] < 0) {
		Point point;
		Util::MakePoint(point, sim, point_num);
		point.Initialize(sim);
		omp_set_lock(&(lock[sim.param.pnum]));
		if (god_pts[point_num] < 0) {
			ptv.push_back(point);
			god_pts[point_num] = ptv.size()-1;
		}
		omp_unset_lock(&(lock[sim.param.pnum]));
	}
	return;
}

void PINT::GodToPtv(Point * const ptv_master, vector<int>& god_pts, std::deque<Point>& ptv, Simdat& sim, vector<omp_lock_t>& lock) {
	for (int p = 0; p < ptv.size(); p++) {
		if (ptv[p].get_t_last_liq()>0.0) {
			int p_real = Util::ijk_to_p(ptv[p].get_i(), ptv[p].get_j(), ptv[p].get_k(), sim);
			// only one thread in here at a time (per ptv_master point); others block at set_lock()
			omp_set_lock(&(lock[p_real]));
			if (ptv_master[p_real].get_t_last_liq() <= ptv[p].get_t_last_liq()) { ptv_master[p_real] = ptv[p]; } 
			omp_unset_lock(&(lock[p_real]));
		}
	}
}

void PINT::beam_trace(vector<int>& test_pts, vector<int>& god_pts, std::deque<Point>& ptv, vector<omp_lock_t>& lock, vector<path_seg>& segv, Simdat& sim, vector<int>& seg_num, int itert_start, int itert_end) {
	int seg_temp_prev, seg_temp_now;

	if (itert_end >= seg_num.size()) { return; }
	else { seg_temp_prev = seg_num[itert_end]; }

	if (itert_start >= seg_num.size()) { seg_temp_now = seg_num[seg_num.size() - 1]; }
	else{ seg_temp_now = seg_num[itert_start]; }

	int pnum = 0;
	int x_grid_num = 0, y_grid_num = 0, z_grid_num = 0;
	int x_flat = 0, y_flat = 0, z_flat = 0;
	for (int i = seg_temp_prev - 1; i <= seg_temp_now; i++) {
		if (i < 0) { i = 0; }
		if (sim.param.xnum - 1) { x_grid_num = (int)((segv[i].sx - sim.param.xmin) / sim.param.xres); x_flat = 1; }
		if (sim.param.ynum - 1) { y_grid_num = (int)((segv[i].sy - sim.param.ymin) / sim.param.yres); y_flat = 1; }
		if (sim.param.znum - 1) { z_grid_num = (int)((segv[i].sz - sim.param.zmin) / sim.param.zres); z_flat = 1; }

		if (x_grid_num < 0) { x_grid_num = 0; x_flat = 0; }
		else if (x_grid_num >= (sim.param.xnum - 1)) { x_grid_num = sim.param.xnum - 1; x_flat = 0; }
		if (y_grid_num < 0) { y_grid_num = 0; y_flat = 0; }
		else if (y_grid_num >= (sim.param.ynum - 1)) { y_grid_num = sim.param.ynum - 1; y_flat = 0; }
		if (z_grid_num < 0) { z_grid_num = 0; z_flat = 0; }
		else if (z_grid_num >= (sim.param.znum - 1)) { z_grid_num = sim.param.znum - 1; z_flat = 0; }

		vector<int> point_nums;
		for (int a = 0; a <= z_flat; a++) {
			for (int b = 0; b <= y_flat; b++) {
				for (int c = 0; c <= x_flat; c++) {
					pnum = (z_grid_num + a) + sim.param.znum*(y_grid_num + b) + sim.param.znum*sim.param.ynum*(x_grid_num + c);
					PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
					if (!ptv[god_pts[pnum]].get_T_calc_flag()) {
						test_pts.push_back(pnum);
						ptv[god_pts[pnum]].set_T_calc_flag();
					}
				}
			}
		}
	}

	//Same Thing but With Current Time too if 1st iteration
	if (!itert_end) {
		int_seg current_beam = Util::GetBeamLoc(itert_start*sim.param.dt, segv, sim, sim.util.start_seg);
		if (sim.param.xnum - 1) { x_grid_num = (int)((current_beam.xb - sim.param.xmin) / sim.param.xres); x_flat = 1; }
		if (sim.param.ynum - 1) { y_grid_num = (int)((current_beam.yb - sim.param.ymin) / sim.param.yres); y_flat = 1; }
		if (sim.param.znum - 1) { z_grid_num = (int)((current_beam.zb - sim.param.zmin) / sim.param.zres); z_flat = 1; }

		if (x_grid_num < 0) { x_grid_num = 0; x_flat = 0; }
		else if (x_grid_num >= (sim.param.xnum - 1)) { x_grid_num = sim.param.xnum - 1; x_flat = 0; }
		if (y_grid_num < 0) { y_grid_num = 0; y_flat = 0; }
		else if (y_grid_num >= (sim.param.ynum - 1)) { y_grid_num = sim.param.ynum - 1; y_flat = 0; }
		if (z_grid_num < 0) { z_grid_num = 0; z_flat = 0; }
		else if (z_grid_num >= (sim.param.znum - 1)) { z_grid_num = sim.param.znum - 1; z_flat = 0; }

		for (int a = 0; a <= z_flat; a++) {
			for (int b = 0; b <= y_flat; b++) {
				for (int c = 0; c <= x_flat; c++) {
					pnum = (z_grid_num + a) + sim.param.znum*(y_grid_num + b) + sim.param.znum*sim.param.ynum*(x_grid_num + c);
					PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
					if (!ptv[god_pts[pnum]].get_T_calc_flag()) {
						test_pts.push_back(pnum);
						ptv[god_pts[pnum]].set_T_calc_flag();
					}
				}
			}
		}
	}

	return;
}

void PINT::neighbor_check(vector<int>& test_pts, vector<int>& liq_pts, vector<int>& reset_pts, vector<int>& god_pts, std::deque<Point>& ptv, vector<omp_lock_t>& lock, double& t, vector<int_seg>& isegv, Simdat& sim, int num_free, int surface_only) {
	vector<int> test_tmp;

	//identify neighbors of liquid points ONLY ON SURFACE
	#pragma omp parallel num_threads(1+num_free) if(num_free)
	{
		vector<int> th_liq_pts;
		vector<int> th_test_tmp;
		vector<int> th_reset_pts;
		int Tflag = 0;
		#pragma omp for schedule(dynamic,1+test_pts.size()/(1+num_free)/64) 
		//Find neighbors in test_pts for checking
		for (int it = 0; it < test_pts.size(); it++) {
			//Get i, j, k location of current point, then construct array of neighbors
			int i = ptv[god_pts[test_pts[it]]].get_i();
			int j = ptv[god_pts[test_pts[it]]].get_j();
			int k = ptv[god_pts[test_pts[it]]].get_k();

			int n = sim.setting.neighborhood;
			std::vector <int> ijkminmax;
			ijkminmax.push_back(n - i);
			ijkminmax.push_back(i + 1 + n - sim.param.xnum);
			ijkminmax.push_back(n - j);
			ijkminmax.push_back(j + 1 + n - sim.param.ynum);
			ijkminmax.push_back(n - k);
			ijkminmax.push_back(k + 1 + n - sim.param.znum);

			for (int temp = 0; temp < ijkminmax.size(); temp++) {
				if (ijkminmax[temp] < 0) { ijkminmax[temp] = 0; }
			}

			vector<int> nbs;
			int pnum;
			for (int di = -n + ijkminmax[0]; di <= (n - ijkminmax[1]); di++) {
				for (int dj = -n + ijkminmax[2]; dj <= (n - ijkminmax[3]); dj++) {
					for (int dk = -n + ijkminmax[4]; dk <= (n - ijkminmax[5]); dk++) {
						if (surface_only) { pnum = Util::ijk_to_p(i + di, j + dj, k, sim); }
						else { pnum = Util::ijk_to_p(i + di, j + dj, k + dk, sim); }
						PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
						omp_set_lock(&(lock[pnum]));
						if (!ptv[god_pts[pnum]].get_T_calc_flag()) { ptv[god_pts[pnum]].set_T_calc_flag(); Tflag = 1; }
						omp_unset_lock(&(lock[pnum]));
						if (Tflag) {
							if (ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
								th_liq_pts.push_back(pnum);
								th_test_tmp.push_back(pnum);
							}
							th_reset_pts.push_back(pnum);
						}
						Tflag = 0;
						if (surface_only) { break; }
					}
				}
			}
		}
		omp_set_lock(&(lock[sim.param.pnum]));
		{
			liq_pts.insert(liq_pts.end(), th_liq_pts.begin(), th_liq_pts.end());
			test_tmp.insert(test_tmp.end(), th_test_tmp.begin(), th_test_tmp.end());
			reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
		}
		omp_unset_lock(&(lock[sim.param.pnum]));
	}
	test_pts.clear();
	test_pts = test_tmp;
}

void PINT::calc_depth(vector<int>& depths, vector<int>& liq_pts, vector<int>& reset_pts, vector<int>& god_pts, std::deque<Point>& ptv, vector<omp_lock_t>& lock, double& t, vector<int_seg>& isegv, vector<int_seg>& isegv_last, vector<path_seg>& segv, Simdat& sim, int num_free) {
	if (sim.param.znum == 1) { return; }
	#pragma omp parallel num_threads(1+num_free) if(num_free)
	{
		vector<int> th_reset_pts;
		#pragma omp for schedule(dynamic)
		for (int it = 0; it < liq_pts.size(); it++) {
			int i = ptv[god_pts[liq_pts[it]]].get_i();
			int j = ptv[god_pts[liq_pts[it]]].get_j();
			int dnum = i * sim.param.ynum + j;
			int depth = depths[dnum];
			int flag = 1;
			while (true) {
				depth++;
				if (depth == sim.param.znum) { break; }
				int pnum = Util::ijk_to_p(i, j, sim.param.znum - 1 - depth, sim);
				PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
				th_reset_pts.push_back(pnum);
				if (ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) { break; }
				else { flag = 0; }
			}
			while (true) {
				depth--;
				if (depth == 0 || flag == 0) { break; }
				int pnum = Util::ijk_to_p(i, j, sim.param.znum - 1 - depth, sim);
				PINT::GodCheck(ptv, god_pts, lock, sim, pnum);
				th_reset_pts.push_back(pnum);
				if (ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) < sim.mat.T_liq) {
					ptv[god_pts[pnum]].set_output_flag(1);
					if (depth != depths[dnum]) { ptv[god_pts[pnum]].Temp_Calc_Pre_Path(t - sim.param.dt, isegv_last, sim, -1, 0); }
					ptv[god_pts[pnum]].Solidify(t, segv, sim);
				}
				else { break; }
			}
			if (depth == sim.param.znum) { depth--; }
			depths[dnum] = depth;
		}
		omp_set_lock(&(lock[sim.param.pnum]));
		{
			reset_pts.insert(reset_pts.end(), th_reset_pts.begin(), th_reset_pts.end());
		}
		omp_unset_lock(&(lock[sim.param.pnum]));
	}
	return;
}

void PINT::calcSpeedPow(double& power, Simdat& sim) {
	int tot_threads = sim.setting.thnum;
	int	endItert = int(sim.util.approxEndTime/sim.param.dt);

	while (true) {
		int itert_start_last = int(endItert*pow((tot_threads - 1) / double(tot_threads), power));
		int itert_end_last = int(endItert*pow((tot_threads) / double(tot_threads), power));
		if (power >= 1.0) { break; }
		if (itert_start_last == itert_end_last) { power += 0.001; }
		else { break; }
	}
	return;
}

void PINT::calc_iterts(int& itert, int& itert_end, Simdat&  sim, int thread, double power) {
	int tot_threads = sim.setting.thnum;
	int	endItert = int(sim.util.approxEndTime/sim.param.dt);
	
	itert = int(endItert*pow(thread/double(tot_threads), power));
	itert_end = int(endItert*pow((thread + 1) / double(tot_threads), power));

	return;
}