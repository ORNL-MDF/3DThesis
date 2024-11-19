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

string Util::ZeroPadNumber(const int num)
{
	return ZeroPadNumber(num, 7);
}

string Util::ZeroPadNumber(const int num, const int width)
{
	std::ostringstream ss;
	ss << std::setw(width) << std::setfill('0') << num;
	return ss.str();
}

void Util::SetLocks(vector<omp_lock_t>& lock, const Simdat& sim) {
	for (int p = 0; p < lock.size(); p++) { omp_init_lock(&(lock[p])); }
	return;
}

void Util::AddToNodes(Nodes& nodes, const int_seg seg) {
	nodes.size++;
	nodes.xb.push_back(seg.xb);
	nodes.yb.push_back(seg.yb);
	nodes.zb.push_back(seg.zb);
	// TODO::FEATURE
	nodes.phix.push_back(1.0/seg.phix);
	nodes.phiy.push_back(1.0/seg.phiy);
	nodes.phiz.push_back(1.0/seg.phiz);
	nodes.dtau.push_back(seg.dtau);
	nodes.expmod.push_back(log(seg.qmod)-log(seg.phix * seg.phiy * seg.phiz) / 2.0);
}

void Util::CombineNodes(Nodes& nodes, const Nodes& nodes2) {
	nodes.size += nodes2.size;
	nodes.xb.insert(nodes.xb.end(), nodes2.xb.begin(), nodes2.xb.end());
	nodes.yb.insert(nodes.yb.end(), nodes2.yb.begin(), nodes2.yb.end());
	nodes.zb.insert(nodes.zb.end(), nodes2.zb.begin(), nodes2.zb.end());
	nodes.phix.insert(nodes.phix.end(), nodes2.phix.begin(), nodes2.phix.end());
	nodes.phiy.insert(nodes.phiy.end(), nodes2.phiy.begin(), nodes2.phiy.end());
	nodes.phiz.insert(nodes.phiz.end(), nodes2.phiz.begin(), nodes2.phiz.end());
	nodes.dtau.insert(nodes.dtau.end(), nodes2.dtau.begin(), nodes2.dtau.end());
	nodes.expmod.insert(nodes.expmod.end(), nodes2.expmod.begin(), nodes2.expmod.end()); 
}

void Util::ClearNodes(Nodes& nodes) {
	nodes.size = 0;
	nodes.xb.clear();
	nodes.yb.clear();
	nodes.zb.clear();
	nodes.phix.clear();
	nodes.phiy.clear();
	nodes.phiz.clear();
	nodes.dtau.clear();
	nodes.expmod.clear();
}

void Util::Calc_AllScansEndTime(Simdat& sim) {
	for (const vector<path_seg>& path : sim.paths) {
		sim.util.allScansEndTime = max(sim.util.allScansEndTime, path.back().seg_time);
	}
}

void Util::Calc_ScanBounds(Domain& domain, const vector<vector<path_seg>>& paths) {
	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;
	double zmin = DBL_MAX;
	double zmax = -DBL_MAX;

	for (const vector<path_seg>& path : paths) {
		for (const path_seg& seg : path) {
			if (seg.sqmod > 0) {
				if (seg.sx < xmin) { xmin = seg.sx; }
				if (seg.sx > xmax) { xmax = seg.sx; }
				if (seg.sy < ymin) { ymin = seg.sy; }
				if (seg.sy > ymax) { ymax = seg.sy; }
				if (seg.sz < zmin) { zmin = seg.sz; }
				if (seg.sz > zmax) { zmax = seg.sz; }
			}
		}
	}

	if (domain.xmin == -DBL_MAX) { domain.xmin = xmin - 500e-6; }
	if (domain.xmax == DBL_MAX) { domain.xmax = xmax + 500e-6; }
	if (domain.ymin == -DBL_MAX) { domain.ymin = ymin - 500e-6; }
	if (domain.ymax == DBL_MAX) { domain.ymax = ymax + 500e-6; }
	if (domain.zmin == -DBL_MAX) { domain.zmin = zmin - 250e-6; }
	if (domain.zmax == DBL_MAX) { domain.zmax = zmax; }

}

void Util::Calc_NonD_dt(vector<Beam>& beams, const Material& material) {
	for (Beam& beam : beams) {
		beam.nond_dt = beam.ax * beam.ax / material.a;
	}
	return;
}

void Util::Calc_RMax (Simdat& sim){
	sim.settings.t_hist = 1.0 / sim.settings.t_hist;
	if (sim.settings.r_max<0.0) {
		for (const Beam& beam : sim.beams) {
			//If the temperature never gets to 1/t_hist the peak temperature
			if (sim.settings.t_hist < exp(3.0 / 2.0)) { sim.settings.r_max = beam.ax * sqrt(log(sim.settings.t_hist) / 3.0); }
			else { sim.settings.r_max = beam.ax * pow(sim.settings.t_hist, (1.0 / 3.0)) / sqrt(2.0 * exp(1.0)); }
			//If the power never gets to x (K/s)
			double beta = pow(3.0 / 3.14159, 1.5) * beam.q / (sim.material.rho * sim.material.cps);
			double temp_diff = sim.material.T_liq - sim.material.T_init;
			double x = temp_diff * sim.settings.p_hist;
			double r_max_2;
			if (beta / (x * beam.ax * beam.ax * beam.ax) < exp(3.0 / 2.0)) { r_max_2 = beam.ax * sqrt(log(beta / (x * beam.ax * beam.ax * beam.ax)) / 3.0); }
			else { r_max_2 = pow(beta / x, (1.0 / 3.0)) / sqrt(2.0 * exp(1.0)); }

			//Choose the greater of the two
			sim.settings.r_max = max(sim.settings.r_max, r_max_2);
			//sim.r_max = sim.ax*pow(sim.t_hist, (1.0 / 3.0)) / sqrt(2.0*exp(1.0)); 
		}
	}
	return;
}

bool Util::InRMax(const double x, const double y, const Domain& domain, const Settings& settings) {
	if ((x > (domain.xmax + settings.r_max)) || (x < (domain.xmin - settings.r_max))) {return false;}
	else if ((y >(domain.ymax + settings.r_max)) || (y < (domain.ymin - settings.r_max))) {return false; }
	else {return true;}
}

double Util::t0calc(const double t, const Beam& beam, const Material& material, const Settings& settings) {
	//Time for beam peak to decay to a fraction (t_hist) of it's initial power
	const double t_hist_t = beam.nond_dt / 12.0*(pow(settings.t_hist, (2.0 / 3.0)) - 1); 

	//Time for beam to never exert more than a fraction (p_hist) of the difference between the preheat and solidus temperature
	const double beta = pow(3 / 3.14159, 1.5) * beam.q / (material.rho * material.cps);
	const double temp_diff = material.T_liq - material.T_init;
	const double x = temp_diff * settings.p_hist;
	const double y = 432.0*t*(x*x)*(material.a*material.a*material.a) / (beta*beta);
	const double p_hist_t = t / ((1.0 + sqrt(y))*(1.0 + sqrt(y)));  

	double t0 = t - max(t_hist_t, p_hist_t);
	if (t0 < 0.0) { t0 = 0.0; }
	t0 = 0.0;

	return t0;
}

double Util::GetRefTime(const double tpp, const int seg, const vector<path_seg>& path, const Beam& beam) {
	
	double ref_t;
	const double spp = max(tpp / beam.nond_dt, 0.0);
	
	//Sets maximum time for line mode (derived from diffusion distance)
	if (path[seg].smode == 0) {
		double t0 = 0.58870501125 * beam.ax / path[seg].sparam; // sqrt(log(sqrt(2)))~0.58870501125
		ref_t = t0 * sqrt(12.0 * spp + 1.0);
	}
	//Sets maximum time for spot mode (equal to spot time)
	else if (path[seg].smode == 1) {
		ref_t = path[seg].seg_time - path[seg - 1].seg_time;
		if (ref_t == 0) {ref_t = 1.0e-9;}
	}

	return ref_t;
}

int_seg	Util::GetBeamLoc(const double time, const int seg, const vector<path_seg>& path, const Simdat& sim) {

	int_seg current_seg;
	//Location calculation for spot mode
	if (path[seg].smode) {	
		current_seg.xb = path[seg].sx;
		current_seg.yb = path[seg].sy;
		current_seg.zb = path[seg].sz;
	}
	//Location calculation for line mode
	else {							
		const double dx = path[seg].sx - path[seg - 1].sx;
		const double dy = path[seg].sy - path[seg - 1].sy;
		const double dz = path[seg].sz - path[seg - 1].sz;
		const double tcur = time - path[seg - 1].seg_time;
		const double dt_cur = path[seg].seg_time - path[seg - 1].seg_time;
		current_seg.xb = path[seg - 1].sx + (tcur / dt_cur)*dx;
		current_seg.yb = path[seg - 1].sy + (tcur / dt_cur)*dy;
		current_seg.zb = path[seg - 1].sz + (tcur / dt_cur)*dz;
	}

	// If we are sufficiently outside the domain, set power to zero (so it won't be added to integration)
	if (Util::InRMax(current_seg.xb,current_seg.yb,sim.domain,sim.settings)){ current_seg.qmod = path[seg].sqmod; }
	else { current_seg.qmod = 0.0; }

	return current_seg;
}

bool Util::sim_finish(const double t, const Simdat& sim, const int liq_num) {
	bool isDone = false;
	if (t > sim.util.allScansEndTime && liq_num == 0) {isDone = true;}
	return isDone;
}

vector<vector<double>> Util::rotateField(const vector<vector<double>>& df, double angle, const int x, const int y){
    vector<vector<double>> df_rot(df);
    for (size_t i = 0; i < df_rot.size(); i++){
        df_rot[i][x] = df[i][x] * cos(-angle) - df[i][y] * sin(-angle);
        df_rot[i][y] = df[i][x] * sin(-angle) + df[i][y] * cos(-angle);
    }
    return df_rot;
}

double Util::getMax(const vector<vector<double>>& df, int index){
    return (*max_element(df.begin(), df.end(), [index](const vector<double>& a, const vector<double>& b) {
        return a[index] < b[index];
    }))[index];

}

double Util::getMin(const vector<vector<double>>& df, int index) {
    return (*min_element(df.begin(), df.end(), [index](const vector<double>& a, const vector<double>& b) {
        return a[index] < b[index];
    }))[index];
}

std::array<double, 4> Util::getLengthWidthOrigin(const vector<vector<double>>& df, double resolution, const int x, const int y){
    double length = Util::getMax(df, x) - Util::getMin(df, x) + resolution;
    double width = Util::getMax(df, y) - Util::getMin(df, y) + resolution;
    std::array<double, 2> origin = {getMin(df, x), getMin(df, y)};
    std::array<double, 4> stats = {length, width, origin[0], origin[1]};
    return stats;
}


double Util::getPerBoxMelted(const vector<vector<double>>& df, double length, double width, double resolution){
    double melted_area = df.size() * pow(resolution, 2);
    double box_area = length * width;
    double per_box_melted = melted_area / box_area * 100;
    if (per_box_melted <= 100.0) {return per_box_melted;}
    else {return std::numeric_limits<double>::quiet_NaN();}
}
