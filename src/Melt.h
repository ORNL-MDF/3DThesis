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
#include "DataStructs.h"
#include "Grid.h"
#include "Util.h"

using std::vector;
using std::string;

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
	void calc_mp_info(const vector<int>&, const vector<int>&, Grid&, const Simdat&, const double);
}

