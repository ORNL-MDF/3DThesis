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
#include <vector>
#include "impl/Structs/DataStructs.hpp"

using std::vector;
using std::string;

namespace Thesis::impl {
	namespace Calc{
			// Get quadrature nodes and points for several subsequent timesteps in parallel
		void Integrate_Parallel(Nodes&, vector<Nodes>&, vector<int>&, const Simdat&, const double, const bool);
		// Get quadrature nodes and points for one time
		void Integrate_Serial(Nodes&, vector<int>&, const Simdat&, const double, const bool);

		// Adaptive Integration Scheme
		void GaussIntegrate(Nodes&, vector<int>&, const Simdat&, const double, const bool);
		// Adaptive Integration Scheme with a compression scheme between neighboring path segments. Very usefull for point rasters.
		void GaussCompressIntegrate(Nodes&, vector<int>&, const Simdat&, const double, const bool);

		// Adds Simple boundary conditions (x and y) via method of images
		void AddBCs(Nodes&, const Domain&);
	}
}