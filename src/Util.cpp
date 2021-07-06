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

#include <iomanip>
#include <sstream>
#include <cmath>

#include <omp.h>
#include <cmath>

#include "Util.h"
#include "Calc.h"
#include "DataStructs.h"
#include "Point.h"

string Util::ZeroPadNumber(int num)
{
	std::ostringstream ss;
	ss << std::setw(7) << std::setfill('0') << num;
	return ss.str();
}

int		Util::ijk_to_p(int i, int j, int k, Simdat& sim) {
	return i * (sim.param.znum * sim.param.ynum) + j * sim.param.znum + k;
}

void	Util::MakePoint(Point& point, Simdat& sim, int p) {
	//Get i,j,k from p
	double xp, yp, zp;
	int i, j, k;

	k = p % sim.param.znum;
	p = (p - k) / sim.param.znum;
	j = p % sim.param.ynum;
	p = (p - j) / sim.param.ynum;
	i = p;

	if (sim.param.xnum == 1) { xp = sim.param.xmax; }
	else { xp = sim.param.xmin + ((float)i * ((sim.param.xmax - sim.param.xmin) / ((float)sim.param.xnum - 1))); }

	if (sim.param.ynum == 1) { yp = sim.param.ymax; }
	else { yp = sim.param.ymin + ((float)j * ((sim.param.ymax - sim.param.ymin) / ((float)sim.param.ynum - 1))); }

	if (sim.param.znum == 1) { zp = sim.param.zmax; }
	else { zp = sim.param.zmin + ((float)k * ((sim.param.zmax - sim.param.zmin) / ((float)sim.param.znum - 1))); }

	point.set_i(i);
	point.set_j(j);
	point.set_k(k);
	point.set_xloc(xp);
	point.set_yloc(yp);
	point.set_zloc(zp);
}

void	Util::SetLocks(vector<omp_lock_t>& lock, Simdat& sim) {
	for (int p = 0; p < lock.size(); p++) { omp_init_lock(&(lock[p])); }
	return;
}

void	Util::CalcRMax (Simdat& sim){
	sim.setting.t_hist = 1.0 / sim.setting.t_hist;
	if (sim.setting.r_max<0.0) {
		//If the temperature never gets to 1/t_hist the peak temperature
		if (sim.setting.t_hist < exp(3.0 / 2.0)) { sim.setting.r_max = sim.beam.ax*sqrt(log(sim.setting.t_hist) / 3.0); } 
		else { sim.setting.r_max = sim.beam.ax*pow(sim.setting.t_hist, (1.0 / 3.0)) / sqrt(2.0*exp(1.0)); }
		//If the power never gets to x (K/s)
		double beta = pow(3 / 3.14159, 1.5) * sim.beam.q / (sim.mat.rho * sim.mat.cps); 
		double temp_diff = sim.mat.T_liq - sim.mat.Tinit;
		double x = temp_diff * sim.setting.p_hist;
		double r_max_2;
		if (beta / (x*sim.beam.ax*sim.beam.ax*sim.beam.ax) < exp(3.0 / 2.0)) { r_max_2 = sim.beam.ax*sqrt(log(beta / (x*sim.beam.ax*sim.beam.ax*sim.beam.ax)) / 3.0); }
		else { r_max_2 = pow(beta / x, (1.0 / 3.0)) / sqrt(2.0*exp(1.0)); }

		//Choose the greater of the two
		sim.setting.r_max = fmax(sim.setting.r_max, r_max_2);
		//sim.r_max = sim.ax*pow(sim.t_hist, (1.0 / 3.0)) / sqrt(2.0*exp(1.0)); 
	}
	return;
}

bool	Util::InRMax(double x, double y, Simdat& sim) {
	if ((x > (sim.param.xmax + sim.setting.r_max)) || (x < (sim.param.xmin - sim.setting.r_max))) {return false;}
	else if ((y >(sim.param.ymax + sim.setting.r_max)) || (y < (sim.param.ymin - sim.setting.r_max))) {return false; }
	else {return true;}
}

void	Util::InitStartSeg(vector<int>& seg_num, vector<path_seg> segv, Simdat& sim) {
	double f_time = sim.util.scanEndTime;
	int seg_n = 0;
	double a = 0.0;
	while (a <= f_time) {
		while (segv[seg_n].seg_time<a) { // Toleratnce of 1ns
			seg_n++;
		}
		seg_num.push_back(seg_n);
		a += sim.param.dt;
	}
	seg_num.push_back(segv.size() - 1);
	return;
}

void	Util::GetStartSeg(Simdat& sim,vector<int>& seg_num, int itert) {
	if (itert < seg_num.size()) { sim.util.start_seg = seg_num[itert]; }
	else { sim.util.start_seg = seg_num[seg_num.size() - 1]; }
	if (!sim.util.start_seg) { sim.util.start_seg = 1; }
	return;
}

double	Util::t0calc(double t, Simdat& sim) {
	//Time for beam peak to decay to a fraction (t_hist) of it's initial power
	double t_hist_t = sim.util.nond_dt / 12.0*(pow(sim.setting.t_hist, (2.0 / 3.0)) - 1); 

	//Time for beam to never exert more than a fraction (p_hist) of the difference between the preheat and solidus temperature
	double beta = pow(3 / 3.14159, 1.5) * sim.beam.q / (sim.mat.rho * sim.mat.cps);
	double temp_diff = sim.mat.T_liq - sim.mat.Tinit;
	double x = temp_diff * sim.setting.p_hist;
	double y = 432.0*t*(x*x)*(sim.mat.a*sim.mat.a*sim.mat.a) / (beta*beta);
	double p_hist_t = t / ((1.0 + sqrt(y))*(1.0 + sqrt(y)));  

	double t0 = t - fmax(t_hist_t, p_hist_t);
	if (t0 < 0.0) { t0 = 0.0; }
	t0 = 0.0;

	return t0;
}

double	Util::GetRefTime(double& spp, vector<path_seg>& segv, Simdat& sim, int& seg) {
	
	if (spp < 0) { spp = 0; }

	double ref_t, dt_cur;

	dt_cur = segv[seg].seg_time - segv[seg - 1].seg_time;

	if (segv[seg].smode) {	//Sets maximum time for spot mode (equal to spot time)
		ref_t = dt_cur;
		if (ref_t == 0) {ref_t = 1.0e-9;}
	}
	else {	//Sets maximum time for line mode (derived from diffusion distance)
		double t0 = 0.59 *sim.beam.ax / segv[seg].sparam; // sqrt(log(sqrt(2)))~0.59
		ref_t = t0 * sqrt(12 * spp + 1);
	}
	return ref_t;
}

double	Util::GetRefTimeShape(double& spp, infBeam& beam, Simdat& sim, int& seg) {

	if (spp < 0) { spp = 0; }

	double ref_t, dt_cur;

	dt_cur = beam.ssegv[seg].seg_time - beam.ssegv[seg - 1].seg_time;

	if (beam.ssegv[seg].smode) {	//Sets maximum time for spot mode (equal to spot time)
		ref_t = dt_cur;
		if (ref_t == 0) { ref_t = 1.0e-9; }
	}
	else {	//Sets maximum time for line mode (derived from diffusion distance)
		double t0 = 0.59 * beam.min_axy / beam.ssegv[seg].sparam; // sqrt(log(sqrt(2)))~0.59
		ref_t = t0 * sqrt(12 * spp + 1);
	}
	return ref_t;
}

int_seg	Util::GetBeamLoc(double time, vector<path_seg>& segv, Simdat& sim, int& ref_seg_start) {

	int ref_seg = ref_seg_start;
	int flag = 1;
	double dx, dy, dz, tcur, dt_cur;
	int_seg current_seg;
	//int ref_seg = sim.start_seg;
	int seg = 0;

	//Skips binary search, starts at known segment
	while (flag && ref_seg) {
		if (time < segv[ref_seg - 1].seg_time) {
			ref_seg--;
		}
		else {
			seg = ref_seg;
			flag = 0;
		}
	}

	if (segv[seg].smode) {	//Location calculation for spot mode
		current_seg.xb = segv[seg].sx;
		current_seg.yb = segv[seg].sy;
		current_seg.zb = segv[seg].sz;
	}
	else {							//Location calculation for line mode
		dx = segv[seg].sx - segv[seg - 1].sx;
		dy = segv[seg].sy - segv[seg - 1].sy;
		dz = segv[seg].sz - segv[seg - 1].sz;
		tcur = time - segv[seg - 1].seg_time;
		dt_cur = segv[seg].seg_time - segv[seg - 1].seg_time;
		current_seg.xb = segv[seg - 1].sx + (tcur / dt_cur)*dx;
		current_seg.yb = segv[seg - 1].sy + (tcur / dt_cur)*dy;
		current_seg.zb = segv[seg - 1].sz + (tcur / dt_cur)*dz;
	}

	if (Util::InRMax(current_seg.xb,current_seg.yb,sim)){ current_seg.qmod = segv[seg].sqmod; }
	else { current_seg.qmod = 0.0; }

	return current_seg;
}

int_shape_seg	Util::GetBeamLocShape(double time, vector<path_shape_seg>& segv, Simdat& sim, int& ref_seg_start) {
	int ref_seg = ref_seg_start;
	int flag = 1;
	double dx, dy, dz, tcur, dt_cur;
	double dax, day, daz;
	int_shape_seg current_seg;
	//int ref_seg = sim.start_seg;
	int seg = 0;

	//Skips binary search, starts at known segment
	while (flag && ref_seg) {
		if (time < segv[ref_seg - 1].seg_time) {
			ref_seg--;
		}
		else {
			seg = ref_seg;
			flag = 0;
		}
	}

	if (segv[seg].smode) {	//Location calculation for spot mode
		current_seg.xb = segv[seg].sx;
		current_seg.yb = segv[seg].sy;
		current_seg.zb = segv[seg].sz;
		current_seg.ax = segv[seg].ax;
		current_seg.ay = segv[seg].ay;
		current_seg.az = segv[seg].az;
	}
	else {							//Location calculation for line mode
		dx = segv[seg].sx - segv[seg - 1].sx;
		dy = segv[seg].sy - segv[seg - 1].sy;
		dz = segv[seg].sz - segv[seg - 1].sz;
		tcur = time - segv[seg - 1].seg_time;
		dt_cur = segv[seg].seg_time - segv[seg - 1].seg_time;
		current_seg.xb = segv[seg - 1].sx + (tcur / dt_cur) * dx;
		current_seg.yb = segv[seg - 1].sy + (tcur / dt_cur) * dy;
		current_seg.zb = segv[seg - 1].sz + (tcur / dt_cur) * dz;
		dax = segv[seg].ax - segv[seg - 1].ax;
		day = segv[seg].ay - segv[seg - 1].ay;
		daz = segv[seg].az - segv[seg - 1].az;
		current_seg.ax = segv[seg - 1].ax + (tcur / dt_cur) * dax;
		current_seg.ay = segv[seg - 1].ay + (tcur / dt_cur) * day;
		current_seg.az = segv[seg - 1].az + (tcur / dt_cur) * daz;
	}

	if (Util::InRMax(current_seg.xb, current_seg.yb, sim)) { current_seg.qmod = segv[seg].sqmod; }
	else { current_seg.qmod = 0.0; }

	return current_seg;
}

void	Util::EstimateEndTime(Simdat& sim, vector<path_seg>& segv) {
	if (!sim.param.use_PINT) { sim.util.approxEndTime = sim.util.scanEndTime; return; }
	vector<Point> points;
	//Add centroid point
	double sum_tp = 0, sum_xtp = 0, sum_ytp = 0;
	double x_av, y_av;
	for (int seg = 1; seg < segv.size(); seg++) {
		if (segv[seg].sqmod > 0.0) {
			double dt = segv[seg].seg_time - segv[seg - 1].seg_time;
			sum_tp += segv[seg].sqmod*dt;
			if (segv[seg].smode) {
				sum_xtp += segv[seg].sx*segv[seg].sqmod*dt;
				sum_ytp += segv[seg].sy*segv[seg].sqmod*dt;
			}
			else {
				sum_xtp += (segv[seg].sx + segv[seg - 1].sx) / 2.0*segv[seg].sqmod*dt;
				sum_ytp += (segv[seg].sy + segv[seg - 1].sy) / 2.0*segv[seg].sqmod*dt;
			}
		}
	}
	x_av = sum_xtp / sum_tp;
	y_av = sum_ytp / sum_tp;

	Point pt_temp;
	pt_temp.set_xloc(x_av);
	pt_temp.set_yloc(y_av);
	pt_temp.set_zloc(0.0);
	pt_temp.Initialize(sim);
	points.push_back(pt_temp);

	//Add centroid of last THNUM path segments with power
	for (int seg = segv.size() - 1; seg > 0; seg--) {
		if (points.size() == sim.setting.thnum) { break; }
		if (segv[seg].sqmod > 0.0) {
			if (segv[seg].smode) {
				pt_temp.set_xloc(segv[seg].sx);
				pt_temp.set_yloc(segv[seg].sy);
			}
			else {
				pt_temp.set_xloc((segv[seg].sx + segv[seg - 1].sx) / 2.0);
				pt_temp.set_yloc((segv[seg].sy + segv[seg - 1].sy) / 2.0);
			}
			points.push_back(pt_temp);
		}
	}

	//Simulate them
	sim.util.approxEndTime = sim.util.scanEndTime;
	int p = 1;
	int pastEnd = 0;
	sim.util.start_seg = segv.size() - 1;
	while (true) {
		double t = sim.util.approxEndTime + sim.param.dt*pow(2, p);
		int liq_num = 0;

		//Pre-calculate integration loop information
		vector<int_seg> isegv;
		Calc::GaussIntegrate(isegv, segv, sim, t, 0);
		#pragma omp parallel for num_threads(sim.setting.thnum) schedule(static)
		for (int pnum = 0; pnum < points.size(); pnum++) {
			if (points[pnum].Temp_Calc_Pre_Path(t, isegv, sim, 1, 0) >= sim.mat.T_liq) {
				#pragma omp atomic
				liq_num++;
			}
		}
		if (!liq_num) { 
			if (pastEnd) { p--; }
			else { sim.util.approxEndTime += sim.param.dt*pow(2, p - 1); }
			pastEnd = 1;
		}
		else {
			sim.util.approxEndTime = t;
			if (pastEnd) { p--; }
			else {p++;}
		}
		if (!p) {
			sim.util.approxEndTime += sim.param.dt; // May not be right
			break;
		}	
	}
	return;
}

bool	Util::sim_finish(double t, Simdat& sim, int liq_num) {
	if (t > sim.util.scanEndTime) {
		if (liq_num) {return false;}
		else { return true; }
	}
	return false;
}