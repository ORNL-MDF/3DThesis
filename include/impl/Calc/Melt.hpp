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
#include "impl/Structs/DataStructs.hpp"
#include "impl/Structs/Grid.hpp"
#include "impl/Calc/Util.hpp"

using std::vector;
using std::string;

namespace Thesis::impl{
	namespace Melt {
		// Adds relevant points to check if they are melted or not based on the beam path
		void beam_trace(vector<int>&, Grid&, const Simdat&, const double, const double);
		// Checks the neigbors of points to see if they are melted too
		void neighbor_check(vector<int>&, vector<int>& , vector<int>& , Grid&, vector<omp_lock_t>&, const Nodes&, const Simdat&, const double, const bool);
		// Calculate the depth of a melt at a specific x,y coordinate. Only used in mode_3.
		void calc_depth(vector<int>&, vector<int>&, vector<int>&, Grid&, const Nodes&, const Nodes&, const Simdat&, const double);
		// Calculate the max depth of a melt at a specific x,y coordinate. Only works in mode_3.
		void calc_depth_max(vector<int>&, vector<double>&, vector<int>&, Grid&, const Nodes&, const Simdat&);
		// Calculate meltpool statistics 
		// NOTE:: Will not work well when spatially decomposed 
		// NOTE:: Will not work well with multiple beams
		void calc_mp_info(const vector<int>&, Grid&, const Simdat&, const double);
		void local_neighbor_check(vector<int>&, vector<int>&, Grid&,const Simdat&);
	}
}