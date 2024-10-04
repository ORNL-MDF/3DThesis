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
#include <array>
#include <vector>
#include <array>
#include <algorithm>
#include <omp.h>
#include "DataStructs.h"

using std::vector;
using std::array;
using std::string;
using std::max;

namespace Util {
	// Turns an integer into a zero-padded string
	string ZeroPadNumber(const int);
	string ZeroPadNumber(const int, const int);
	// Turns ijk indices into global point number
	inline int ijk_to_p(const int i, const int j, const int k, const Simdat& sim) {
		return i * (sim.domain.znum * sim.domain.ynum) + j * sim.domain.znum + k;
	}
	// Initializes the locks. They make sure only one thread can access a point at a time
	void	SetLocks(vector<omp_lock_t>&, const Simdat&);
	// Function to easily add integration segment to nodes
	void	AddToNodes(Nodes&, const int_seg);
	// Function to combine nodes
	void	CombineNodes(Nodes&, const Nodes&);
	// Function to cleat all nodes
	void	ClearNodes(Nodes&);
	// Checks if it is inside the maximum radius
	bool	InRMax(const double, const double, const Domain&, const Settings&);
	// Calculates the time to integrate back to
	double	t0calc(const double, const Beam&, const Material&, const Settings&);
	// Gets the maximum allowable step size for a path segment
	double	GetRefTime(const double, const int, const vector<path_seg>&, const Beam&);
	// Finds the current beam location
	int_seg GetBeamLoc(const double, const int, const vector<path_seg>&, const Simdat&);
	// Indicates when all points have solidified and the full scan is over
	bool	sim_finish(const double, const Simdat&, const int);

	// Calculates scan end time
	void	Calc_AllScansEndTime(Simdat&);
	// Calculates bounds of scan path
	void	Calc_ScanBounds(Domain&, const vector<vector<path_seg>>&);
	// Calculates the nondimensional integration time
	void	Calc_NonD_dt(vector<Beam>&, const Material&);
	// Calculates the maximum radius around the domain to be considered
	void	Calc_RMax(Simdat&);

	// Rotates melt pool based on scan angle
	vector<vector<double>> rotateField(const vector<vector<double>>&, double, const int, const int);
	// Calculates min and max elements of a specific vector in the df
	double getMin(const vector<vector<double>>&, int);
	double getMax(const vector<vector<double>>&, int);
	// Calculates length, width and origin of melt pool
	array<double, 4> getLengthWidthOrigin(const vector<vector<double>>&, double, const int, const int);
	// Calculates percentage of the melt pool box that is melted
	double getPerBoxMelted(const vector<vector<double>>&, double, double, double);
	}
