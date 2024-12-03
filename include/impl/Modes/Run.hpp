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

#include <omp.h>
#include <cmath>

#include <iostream>
#include <fstream>

#include <vector>
#include <atomic>
#include "DataStructs.h"
#include "Grid.h"

using std::vector;
using std::string;
using std::atomic;

namespace Run{
	//Chooses between modes
	void Simulate(Grid&, const Simdat&);
	
	//Choose between snapshot modes
	void Snapshots(Grid&, const Simdat&);

	//Snapshot modes
	void Snapshots_NoTracking(Grid&, const Simdat&);
	void Snapshots_Volume(Grid&, const Simdat&);
	void Snapshots_GeometryBounds(Grid&, const Simdat&);

	//Choose between solidification modes
	void Solidify(Grid&, const Simdat&);

	//Solidification Modes
	void Solidify_NoTracking(Grid&, const Simdat&); 
	void Solidify_Volume(Grid&, const Simdat&);
	void Solidify_Surface(Grid&, const Simdat&);

	//Special Mode for Stork Data
	void Stork(Grid&, const Simdat&);
}