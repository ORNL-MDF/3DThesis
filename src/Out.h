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
#include "DataStructs.h"
#include "Grid.h"

#include <cstdint>
#include <deque>

using std::vector;
using std::string;

namespace Out {
	// Writes the progress to the console
	void Progress(const Simdat&, const int);
	// Writes the progress of points to the console
	void Point_Progress(const Simdat&, const int);
	// Write the RRDF
	void RRDF_csv(const Simdat&, const vector<uint32_t>&, const vector<float>&, const vector<float>&);
	void RRDF_bin(const Simdat&, const vector<uint32_t>&, const vector<float>&, const vector<float>&);
}