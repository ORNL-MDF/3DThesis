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

#pragma once
#include <string>
#include <vector>

struct int_seg{
	double xb, yb, zb, taui, dtau, qmod;
};

struct path_seg{
	int smode;
	double sx, sy, sz;	//Segment start coordinates
	double sqmod;		//Segmod power modulation
	double sparam;
	double seg_time;
};

struct parBeam
{
	double Xr, Yr, Pmod;
};

struct FileNames {
	std::string name, mat, beam, path;
	std::string domain, settings;
	std::string	parBeams;
};

struct Material {
	double kon, rho, cps, T_liq, T_sol, Tinit;
	double a;

	//Parameter for columnar to equiaxed transition
	int	calc_CET;
	double cet_a, cet_n, cet_N0;
};

struct Beam {
	double ax, ay, az;
	double eff, q;
};

struct SimParams {
	int xnum, ynum, znum, pnum;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double BC_xmin, BC_xmax, BC_ymin, BC_ymax;
	double xres, yres, zres;
	int mode, sol_finish, out_freq, use_PINT;
	double dt;
};

struct Settings {
	int thnum;
	double dttest, dtau_min, t_hist, p_hist, r_max;
	int max_iter;
	int neighborhood;

	int out_mode;

	int compress;
	int use_BCs;
	int parBeams;
};

struct Utility {
	double nond_dt; 
	int np;
	double scanEndTime;
	double approxEndTime;
	int start_seg;
	int prog_print_last=0;
};

struct Simdat{
	FileNames	file;
	Material	mat;
	Beam		beam;
	SimParams	param;
	Settings	setting;
	Utility		util;

	std::vector<parBeam> beams;
};