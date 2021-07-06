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
#include <omp.h>
#include "DataStructs.h"
#include "Point.h"

using std::vector;
using std::string;

namespace Util {
	// Turns an integer into a zero-padded string
	string ZeroPadNumber(int);
	// Turns ijk indices into global point number
	int		ijk_to_p(int, int, int, Simdat&);
	// Creates a point based on ijk indices
	void	MakePoint(Point&, Simdat&, int);
	// Initializes the locks. They make sure only one thread can access a point at a time
	void	SetLocks(vector<omp_lock_t>&, Simdat&);
	// Calculates the maximum radius around the domain to be considered
	void	CalcRMax (Simdat&);
	// Checks if it is inside the maximum radius
	bool	InRMax(double, double, Simdat&);
	// Calculates the time to integrate back to
	double	t0calc(double, Simdat&);
	// Creates a vector allowing quicker (than binary search) finding of current path segment
	void	InitStartSeg(vector<int>&, vector<path_seg>, Simdat&);
	// Retrieves the initial path segment to look at
	void	GetStartSeg(Simdat&, vector<int>&, int);
	// Gets the maximum allowable step size for a path segment
	double	GetRefTime(double&, vector<path_seg>&, Simdat&, int&);
	double	GetRefTimeShape(double&, infBeam&, Simdat&, int&);
	// Finds the current beam location
	int_seg GetBeamLoc(double, vector<path_seg>&, Simdat&, int&);
	// Finds the current beam location for Inf Beams
	int_shape_seg GetBeamLocShape(double, vector<path_shape_seg>&, Simdat&, int&);
	// Estimates the actual end of the simulation
	void	EstimateEndTime(Simdat&, vector<path_seg>&);
	// Indicates when all points have solidified and the full scan is over
	bool	sim_finish(double, Simdat&, int);
}