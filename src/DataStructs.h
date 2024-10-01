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

#pragma once
#include <string>
#include <vector>
#include <climits>
#include <cfloat>

using std::vector;
using std::string;

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
	bool RDF;
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