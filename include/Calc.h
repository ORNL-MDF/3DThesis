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
#include "DataStructs.h"

using std::vector;
using std::string;

namespace Calc {
	// Chooses the mode of integrations. Only compares methods if it's not for a single point solidifying
	void Integrate(vector<int_seg>&, vector<vector<int_seg>>&, vector<path_seg>&, Simdat&, vector<int>&, int, double, int, int);
	void Integrate_thread(vector<int_seg>&,vector<path_seg>&, Simdat&, double, int);

	// Adaptive Integration Scheme
	void GaussIntegrate(vector<int_seg>&, vector<path_seg>&, Simdat&, double, int);
	// Adaptive Integration Scheme with a compression scheme between neighboring path segments. Very usefull for point rasters.
	void GaussCompressIntegrate(vector<int_seg>&, vector<path_seg>&, Simdat&, double, int);
	// Adaptive Integration Scheme for use with multiple independent beams
	void GaussIntegrateInfBeams(vector<infBeam>&, Simdat&, double, int);

	// Adds Simple boundary conditions (x and y) via method of images
	void AddBCs(vector<int_seg>&, Simdat&);
	
	void UseParBeams(vector<int_seg>&, Simdat&);
}