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

#pragma once
#include <vector>
#include <algorithm>
#include <omp.h>
#include "DataStructs.h"

using std::vector;
using std::string;
using std::max;

namespace Util {
	// Turns an integer into a zero-padded string
	string ZeroPadNumber(const int);
	string ZeroPadNumber(const int, const int);
	// Turns ijk indices into global point number
	int		ijk_to_p(const int, const int, const int, const Simdat&);
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
}