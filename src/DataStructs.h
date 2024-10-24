//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) � Transformer Resilience and Advanced Components (TRAC) Program.

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

#pragma once
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <climits>
#include <cfloat>

using std::vector;
using std::string;
using std::deque;
using std::list;

#define PI 3.14159265358979323846

struct int_seg {
	double xb, yb, zb, phix, phiy, phiz, dtau, qmod;
};

// What is used to integrate
struct Nodes {
	size_t size = 0;
	// {xb,yb,zb} = coords
	// {phix,phiy,phiz} = diffusion
	// {dtau} = node weight
	// {expmod} = frontloads computation
	vector<double> xb, yb, zb, phix, phiy, phiz, dtau, expmod;
};

// What is read in from the paths
struct path_seg{
	int smode;			//Segment mode (line melt | spot melt)
	double sx, sy, sz;	//Segment end coordinates
	double sqmod;		//Segment power modulation
	double sparam;		//Segment time parameter (speed | spot time )
	double seg_time;	//Segment end time
};

// 3D coordinate
struct coord
{
	double x, y, z;
};

// Filename read into simulation
struct FileNames {
	string	name, dataDir, mode, material, beam, path;
	string	domain, output, settings;
};

// Material constants
struct Material {
	double kon; // Thermal Conductivty
	double rho; // Density
	double cps; // Specifc Heat
	double T_liq; // Liquidus Temperature
	double T_init; // Inital Temperature (Preheat/Ambient)
	double a; // Thermal Diffusivity
	double cet_a, cet_n, cet_N0; // Parameters for CET
};

// Beam specific parameters
struct Beam {
	double ax, ay, az; // Beam Shape
	double eff; // Absoprtion Efficiency
	double q; // Beam Power
	double nond_dt; // Nondimensional Time
};

// Collection of simulation parameters
struct SimParams {
	// Mode
	string mode; // Snapshots, TemperatureHistory, Solidification
	
	// Snapshots
	vector<double> SnapshotTimes;

	// Solidification
	string tracking = "None";		// meltpool tracking mode
	double dt = 1e-5;			// timestep
	int out_freq = 1;		// output frequency
	bool secondary = false;		// calculate secondary solidifiaction
};

// Domain paramters
struct Domain {
	// Domain numbers
	int xnum, ynum, znum, pnum;

	// Domain bounds 
	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;
	double zmin = DBL_MAX;
	double zmax = -DBL_MAX;

	// Domain resolution
	double xres, yres, zres;

	// Domain reflections
	bool use_BCs;
	int	BC_reflections;
	double BC_xmin, BC_xmax, BC_ymin, BC_ymax, BC_zmin;

	// Domain point file
	string pointsFile;
	vector<coord> points;
	bool customPoints;
};

// Fundamental settings
struct Settings {
	
	// Itegration settings
	double dtau_min, t_hist, p_hist, r_max;
	
	// Liquidus search via Newton method settings
	int max_iter;
	double dttest;
	
	// Scan Path Compression
	int compress;

	// Compute
	int thnum;
	bool use_PINT;

	// MPI
	bool mpi_overlap;
};

// Some utility variables
struct Utility { 
	double allScansEndTime = 0;		// Time when all scans are done
	double approxEndTime = 0;		// Approximate end time to simulation
	bool sol_finish = false;				// Has solidification finished
	bool do_sol = false;					// Do solidification calculation
};

// What variables to output
struct Output {
	bool x, y, z;
	bool T, T_hist;
	bool tSol, G, Gx, Gy, Gz, V, dTdt, eqFrac, depth, numMelt;
	bool RDF, mp_stats;
	bool H, Hx, Hy, Hz;
};

// Entire simulation data structure
struct Simdat{
	FileNames	files;
	Output		output;

	Material	material;

	vector<Beam> beams;
	vector<vector<path_seg>> paths;

	Domain		domain;
	SimParams	param;
	
	Settings	settings;

	Utility		util;
};