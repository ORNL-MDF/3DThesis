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
#include <deque>
#include "DataStructs.h"
#include "Point.h"

using std::vector;
using std::string;

namespace Run{
	//Chooses between modes
	void Simulate(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);
	
	//Take a Temperature snapshot at end of the path
	void Mode_0(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);

	//Track all points at all times
	void Mode_1(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);

	//Track only meltpool at all times
	void Mode_2(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);
	//Mode_2, parallel in time, godPoints
	void Mode_2_PINT(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);

	//Track only meltpool perimeter at all times
	void Mode_3(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);
	//Mode_3, parallel in time
	void Mode_3_PINT(Point * const ptv, vector<path_seg>&, Simdat&, vector<int>&);

	//Track like 2 surfaces at the same time or something... for Matt's CA code
	/*void Mode_4(Point* const ptv, vector<path_seg>&, Simdat&, vector<int>&);*/
}