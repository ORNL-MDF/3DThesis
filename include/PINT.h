//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright 2019 UT-Battelle, LLC
*
* All Rights Reserved
*
* Authors: Benjamin Stump <stumpbc@ornl.gov>, Alex Plotkowski, James Ferguson, Kevin Sisco
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
#include "Point.h"
#include <deque>

namespace PINT {
	// Check if godpoint references a point, if not, make the point
	void GodCheck(std::deque<Point>&, std::vector<int>&, std::vector<omp_lock_t>&, Simdat&, int);
	// If the point melted after the recorded time, update it
	void GodToPtv(std::vector<Point>&, std::vector<int>&, std::deque<Point>&, Simdat&, std::vector<omp_lock_t>&);
	// Melt::beam_trace with the addition of godpoints AND always looking at the current time as well
	void beam_trace(std::vector<int>&, std::vector<int>&, std::deque<Point>&, std::vector<omp_lock_t>& lock, std::vector<path_seg>&, Simdat&, std::vector<int>&, int, int);
	// Melt::neighbor_check with the addition of godpoints
	void neighbor_check(std::vector<int>&, std::vector<int>&, std::vector<int>&, std::vector<int>&, std::deque<Point>&, std::vector<omp_lock_t>&, double&, std::vector<int_seg>&, Simdat&, int, int);
	// Melt::calc_depth with the addition of godpoints
	void calc_depth(std::vector<int>&, std::vector<int>&, std::vector<int>&, std::vector<int>&, std::deque<Point>&, std::vector<omp_lock_t>&, double&, std::vector<int_seg>&, std::vector<int_seg>&, std::vector<path_seg>&, Simdat&, int);
	// Help fix the scaling power [0,1] to make sure there are no repeat iterations
	void calcSpeedPow(double&, Simdat&);
	// Calculate what iteration individual threads should start on
	void calc_iterts(int&, int&, Simdat&, int, double);
}