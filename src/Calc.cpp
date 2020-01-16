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

#include <cmath>
#include <cfloat>
#include "DataStructs.h"
#include "Calc.h"
#include "Util.h"

void Calc::Integrate(std::vector<int_seg>& isegv, std::vector<std::vector<int_seg>>& isegv_par, std::vector<path_seg>& segv, Simdat& sim, std::vector<int>& seg_num, int itert, double t, int sol, int par_num) {
	Util::GetStartSeg(sim, seg_num, itert);
	if (par_num>1) {
		if (!isegv_par.size()) {
			isegv_par.resize(par_num);
			#pragma omp parallel num_threads(par_num)
			{
				std::vector<int_seg> th_isegv;
				Simdat th_sim = sim;
				#pragma omp for schedule(static)
				for (int i = 0; i < par_num; i++) {
					Util::GetStartSeg(th_sim, seg_num, itert+i);
					Calc::Integrate_thread(th_isegv, segv, th_sim, t+i*sim.param.dt, sol);
					isegv_par[par_num - i - 1] = th_isegv;
				}
			}
		}
		isegv = isegv_par.back();
		isegv_par.pop_back();
	}
	else {Calc::Integrate_thread(isegv, segv, sim, t, sol);}
	return;
}

void Calc::Integrate_thread(std::vector<int_seg>& isegv, std::vector<path_seg>& segv, Simdat& sim, double t, int sol) {
	if (sim.setting.compress) { 
		Calc::GaussCompressIntegrate(isegv, segv, sim, t, sol); 
		//If not solidifying, always choose the minimum of the two; otherwise, too expensive
		/*if (!sol) {
			std::vector<int_seg> isegv_reg;
			Calc::GaussIntegrate(isegv_reg, segv, sim, t, sol);
			if (isegv_reg.size() <= isegv.size()) { isegv = isegv_reg; }
		}*/
	}
	else { Calc::GaussIntegrate(isegv, segv, sim, t, sol); }

	if (sim.setting.parBeams) { Calc::UseParBeams(isegv, sim); }
	if (sim.setting.use_BCs){ Calc::AddBCs(isegv, sim); }

	return;
}

void Calc::GaussIntegrate(std::vector<int_seg>& isegv, std::vector<path_seg>& segv, Simdat& sim, double t, int sol) {
	int seg_temp = sim.util.start_seg; //Stores the starting seg to change it back
	
	while ((t > segv[seg_temp].seg_time) && (t < sim.util.scanEndTime)) { seg_temp++; }
	while (t < segv[seg_temp - 1].seg_time) { seg_temp--; }

	double t0 = Util::t0calc(t, sim); //the time you want to integrate back too

	double nodes[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};
	double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};

	double curStep_max_start = sim.util.nond_dt;
	if (sol) { curStep_max_start *= sim.beam.az; }

	double curStep_max = curStep_max_start;
	double curStep_use = curStep_max;
	int curOrder = 16;

	int_seg current_beam = Util::GetBeamLoc(t, segv, sim, seg_temp); //Make 1st segment at the exact time...for instantaneous heat source additon to laplacian
	current_beam.taui = t;
	current_beam.dtau = 0.0;
	if (t > sim.util.scanEndTime) { current_beam.qmod = 0.0; }
	isegv.push_back(current_beam);

	int tflag = 1;
	double t2 = t;
	double tpp = 0.0;
	double spp = tpp / sim.util.nond_dt;

	if (t > sim.util.scanEndTime) {
		tpp += t- sim.util.scanEndTime;
		spp = tpp / sim.util.nond_dt;
		t2 = sim.util.scanEndTime;
	}

	while (tpp >= 2 * curStep_max - curStep_max_start) {
		curStep_max *= 2.0;
		if (curOrder != 2) {
			curOrder = (curOrder / 2);
		}
	}

	while (tflag) {
		int sflag = 0;
		spp = tpp / sim.util.nond_dt;
		
		double ref_time = Util::GetRefTime(spp, segv, sim, seg_temp);

		if (ref_time < curStep_max) {curStep_use = ref_time;}
		else {curStep_use = curStep_max;}

		double t1 = t2 - curStep_use;
		double next_time = segv[seg_temp - 1].seg_time;

		//If we are at the end of a segment, hit the end of it and set the program to jump to the next segment next time
		if (t1 < next_time) {	
			t1 = next_time;
			if (next_time > t0) {sflag = 1;}
			else {tflag = 0;}
		}

		//Add Quadrature Points
		for (int a = (2 * curOrder - 3); a >(curOrder - 3); a--) {
			double tp = 0.5 * ((t2 - t1)*nodes[a] + (t2 + t1));
			int_seg current_beam = Util::GetBeamLoc(tp, segv, sim, seg_temp);
			current_beam.taui = tp;
			current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
			if (current_beam.qmod > 0.0 && current_beam.dtau) { isegv.push_back(current_beam); }
		}

		//If we are switching segments, set start time to start of next segment and increment the start segment down
		if (sflag) {	
			t2 = next_time;
			tpp += (t2 - next_time);
			seg_temp--;	
		}
		//If we are not switching segments, do the normal thing
		else {	
			t2 -= curStep_use;
			tpp += curStep_use;
		}

		//Increase Maximum Step and Decrease Gauss order if it is "safe" to
		if (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}
	}
	return;
}

void Calc::GaussCompressIntegrate(std::vector<int_seg>& isegv, std::vector<path_seg>& segv, Simdat& sim, double t, int sol) {
	int seg_temp = sim.util.start_seg; //Stores the starting seg to change it back

	while ((t > segv[seg_temp].seg_time) && (t < sim.util.scanEndTime)) { seg_temp++; }
	while (t < segv[seg_temp - 1].seg_time) { seg_temp--; }

	double t0 = Util::t0calc(t, sim); //the time you want to integrate back too

	double nodes[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};
	double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};

	double curStep_max_start = sim.util.nond_dt;
	if (sol) { curStep_max_start *= sim.beam.az; }

	double curStep_max = curStep_max_start;
	double curStep_use = curStep_max;
	int curOrder = 16;

	int_seg current_beam = Util::GetBeamLoc(t, segv, sim, seg_temp); //Make 1st segment at the exact time...for instantaneous heat source additon to laplacian
	current_beam.taui = t;
	current_beam.dtau = 0.0;
	if (t > sim.util.scanEndTime) { current_beam.qmod = 0.0; }
	isegv.push_back(current_beam);

	int tflag = 1;
	double t2 = t;
	double tpp = 0.0;
	double spp = tpp / sim.util.nond_dt;

	if (t > sim.util.scanEndTime) {
		tpp += t - sim.util.scanEndTime;
		spp = tpp / sim.util.nond_dt;
		t2 = sim.util.scanEndTime;
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
		int cflag = 1;

		//If the time is less than t0, break the whole thing
		int quit = 0;
		//If outside r, then keep going down until back in r
		while (true) {
			xs = segv[seg_temp].sx;
			ys = segv[seg_temp].sy;
			ts = segv[seg_temp].seg_time;
			if (ts <= t0) { quit = 1; break; }
			if (!Util::InRMax(xs, ys, sim)) { seg_temp--; }
			else { tpp += t2 - ts; spp = tpp / sim.util.nond_dt; t2 = ts; break; }
		}
		if (quit) { break; }

		while (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}

		double ref_time = Util::GetRefTime(spp, segv, sim, seg_temp);
		if (ref_time < curStep_max) { curStep_use = ref_time; }
		else { curStep_use = curStep_max; }	

		int_seg current_beam_t2 = Util::GetBeamLoc(t2, segv, sim, seg_temp);
		xp = current_beam_t2.xb;
		yp = current_beam_t2.yb;

		// Diffusion Distance Squared (Distance it diffused by *some amount*, squared)
		r2 = log(2.0) / 8.0 * (sim.beam.ax*sim.beam.ax) *(12.0 * (t - t2) * sim.mat.a / (sim.beam.ax*sim.beam.ax) + 1.0);

		seg_temp_2 = seg_temp;

		double t1 = segv[seg_temp - 1].seg_time;

		//If we won't be switching segments, do normal integration
		if (t1 < t2 - curStep_use) { num_comb_segs = 0; t1 = t2 - curStep_use; }
		else {
			int cflag = 1;
			while (cflag) { //If we will be switching segments, do compressed integration
							//If the segment start time is zero, end the loop
				if (segv[seg_temp_2 - 1].seg_time <= t0) { tflag = 0; break; }

				xs = segv[seg_temp_2 - 1].sx;
				ys = segv[seg_temp_2 - 1].sy;
				dist2 = (xs - xp)*(xs - xp) + (ys - yp)* (ys - yp);
				ts = segv[seg_temp_2 - 1].seg_time;

				// If the next segment is outside the calculation domain, break the loop
				if (!Util::InRMax(xp, yp, sim)) { break; }

				// IF 
				//	the distance between endpoints, or points, is bigger than the diffusion distance
				// OR
				//	the qmod's aren't equal (discontinuity) AND current order or time difference is *as defined* (to minimize error)
				// THEN
				//	stop combining segments and finalize
				// ELSE
				//  If the distance of the next segment is less than the diffusion distance, set end time to next segment, average them, and keep going
				if (dist2 > r2 || (segv[seg_temp_2].sqmod != segv[seg_temp_2 - 1].sqmod && (curOrder > 2 || (t2 - ts) > (curStep_max / 64.0)))) {
					//If point, do averaging calculations and set new end time to the next segment
					if (segv[seg_temp_2].smode) {
						sum_qmodtx += segv[seg_temp_2].sx*segv[seg_temp_2].sqmod*(t1 - ts);
						sum_qmodty += segv[seg_temp_2].sy*segv[seg_temp_2].sqmod*(t1 - ts);
						sum_qmodt += segv[seg_temp_2].sqmod*(t1 - ts);
						sum_t += (t1 - ts);
						t1 = ts;
					}
					//If line, find the time when the dist=r. If this time is less than the next segment start time, set end time to next segment; else, set start time to when dist=r
					else {
						dx = segv[seg_temp_2].sx - segv[seg_temp_2 - 1].sx;
						dy = segv[seg_temp_2].sy - segv[seg_temp_2 - 1].sy;
						dt = segv[seg_temp_2].seg_time - ts;

						double t_int = ts + dt * (sqrt(((xp - xs)*(xp - xs) + (yp - ys)*(yp - ys)) / (dx*dx + dy * dy)) - sqrt(r2 / (dx*dx + dy * dy)));
						if (t_int < t2 - curStep_max) { t_int = t2 - curStep_max; }
						sum_qmodtx += (xs + dx * ((t1 + t_int) / 2.0 - ts) / dt)*segv[seg_temp_2].sqmod*(t1 - t_int);
						sum_qmodty += (ys + dy * ((t1 + t_int) / 2.0 - ts) / dt)*segv[seg_temp_2].sqmod*(t1 - t_int);
						sum_qmodt += segv[seg_temp_2].sqmod*(t1 - t_int);
						sum_t += (t1 - t_int);
						t1 = t_int;
					}
					cflag = 0;
				}
				else {
					//If the next time is less than the minimum allowed time, set appropriate times to the minimum allowed time
					if (ts < t2 - curStep_max) {
						cflag = 0;
						if (segv[seg_temp_2].smode) { t1 = t2 - curStep_max; }
						else { t1 = t2 - curStep_max; }
					}
					//Otherwise combine the segments and continue
					else {
						t1 = ts;
						num_comb_segs++;
					}

					//If point, 
					if (segv[seg_temp_2].smode) {
						dt = (segv[seg_temp_2].seg_time - t1);
						
						sum_qmodtx += segv[seg_temp_2].sx*segv[seg_temp_2].sqmod*dt;
						sum_qmodty += segv[seg_temp_2].sy*segv[seg_temp_2].sqmod*dt;
						sum_qmodt += segv[seg_temp_2].sqmod*dt;
						sum_t += dt;
					}
					else {
						dx = segv[seg_temp_2].sx - segv[seg_temp_2 - 1].sx;
						dy = segv[seg_temp_2].sy - segv[seg_temp_2 - 1].sy;
						dt = (segv[seg_temp_2].seg_time - t1);

						sum_qmodtx += (xs + dx * ((segv[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt)*segv[seg_temp_2].sqmod*dt;
						sum_qmodty += (ys + dy * ((segv[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt)*segv[seg_temp_2].sqmod*dt;
						sum_qmodt += segv[seg_temp_2].sqmod*dt;
						sum_t += dt;
					}
				}
				if (t1 == segv[seg_temp_2 - 1].seg_time) {
					xs = segv[seg_temp_2 - 1].sx;
					ys = segv[seg_temp_2 - 1].sy;
					seg_temp_2--;
					if (!Util::InRMax(xs, ys, sim)) { break; }
				}
			}

			if (!num_comb_segs) { //If we didn't combine any segments anyways, then just do normal integration
				t1 = segv[seg_temp - 1].seg_time;
				seg_temp_2 = seg_temp - 1;
			}
		}

		//Add Quadrature Points using the same average for x, y, and qmod (CAN BE IMPROVED)
		if (num_comb_segs) {
			int_seg current_beam;
			if (sum_qmodt > 0.0) {
				current_beam.xb = sum_qmodtx / sum_qmodt;
				current_beam.yb = sum_qmodty / sum_qmodt;
				current_beam.qmod = sum_qmodt / sum_t;
				for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
					double tp = 0.5 * ((t2 - t1)*nodes[a] + (t2 + t1));
					current_beam.taui = tp;
					current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
					if (current_beam.qmod > 0.0 && current_beam.dtau) { isegv.push_back(current_beam); }
				}
			}

		}
		//Add Quadrature Points
		else {
			for (int a = (2 * curOrder - 3); a >(curOrder - 3); a--) {
				double tp = 0.5 * ((t2 - t1)*nodes[a] + (t2 + t1));
				int_seg current_beam = Util::GetBeamLoc(tp, segv, sim, seg_temp);
				current_beam.taui = tp;
				current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
				if (current_beam.qmod > 0.0 && current_beam.dtau) { isegv.push_back(current_beam); }
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
	return;
}

void Calc::AddBCs(std::vector<int_seg>& isegv, Simdat& sim) {
	int org_size = isegv.size();
	std::vector<int_seg> isegv_org = isegv;
	std::vector<int_seg> isegv2=isegv_org;
	
	double xmin, xmax, xStr;
	double ymin, ymax, yStr;
	
	xmin = sim.param.BC_xmin;
	xmax = sim.param.BC_xmax;
	ymin = sim.param.BC_ymin;
	ymax = sim.param.BC_ymax;
	
	if (xmin != DBL_MAX) {
		isegv2 = isegv_org;
		for (int i = 0; i < org_size; i++) {isegv2[i].xb = 2 * xmin - isegv2[i].xb;}
		isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		if (ymin != DBL_MAX) {
			for (int i = 0; i < org_size; i++) {isegv2[i].yb = 2 * ymin - isegv2[i].yb;}
			isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		}
	}
	if (xmax != DBL_MAX) {
		isegv2 = isegv_org;
		for (int i = 0; i < org_size; i++) {isegv2[i].xb = 2 * xmax - isegv2[i].xb;}
		isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		if (ymax != DBL_MAX) {
			for (int i = 0; i < org_size; i++) {isegv2[i].yb = 2 * ymax - isegv2[i].yb;}
			isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		}
	}
	if (ymin != DBL_MAX) {
		isegv2 = isegv_org;
		for (int i = 0; i < org_size; i++) {isegv2[i].yb = 2 * ymin - isegv2[i].yb;}
		isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		if (xmax != DBL_MAX) {
			for (int i = 0; i < org_size; i++) {isegv2[i].xb = 2 * xmax - isegv2[i].xb;}
			isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		}
	}
	if (ymax != DBL_MAX) {
		isegv2 = isegv_org;
		for (int i = 0; i < org_size; i++) {isegv2[i].yb = 2 * ymax - isegv2[i].yb;}
		isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		if (xmin != DBL_MAX) {
			for (int i = 0; i < org_size; i++) {isegv2[i].xb = 2 * xmin - isegv2[i].xb;}
			isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
		}
	}

	return;
}

void Calc::UseParBeams(std::vector<int_seg>& isegv, Simdat& sim) {
	int org_size = isegv.size();
	std::vector<int_seg> isegv_org = isegv;
	std::vector<int_seg> isegv2 = isegv_org;
	isegv.clear();

	int num_beams = sim.beams.size();

	for (int b = 0; b < num_beams; b++) {
		isegv2 = isegv_org;
		for (int i = 0; i < org_size; i++) {
			isegv2[i].xb += sim.beams[b].Xr;
			isegv2[i].yb += sim.beams[b].Yr;
			isegv2[i].qmod *= sim.beams[b].Pmod;
		}
		isegv.insert(isegv.end(), isegv2.begin(), isegv2.end());
	}

	return;
}