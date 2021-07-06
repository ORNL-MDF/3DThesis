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
#include "Point.h"

using std::vector;
using std::string;

namespace Init {
	// Read for Keywords
	void	Keywords_Lv2(vector<string>&, vector<vector<string>>&, vector<vector<string>>&, string&);
	void	Keywords_Lv2(vector<string>&, vector<vector<string>>&, vector<vector<double>>&, string&);
	// Helps Set Values and Display Errors
	void	SetValues(string&, string, string, string, int);
	void	SetValues(int&,  double, int, string, int);
	void	SetValues(double&, double, double, string, int);

	// Gets file names from the input file
	void	GetFileNames(Simdat&,string);
	void	ReadSimParams(vector<path_seg>&, Simdat&);
	void	MakeDataDirectory(Simdat&);

	// Reads in all necessary simulation parameters
	void	FileRead_Material(Simdat&);
	void	FileRead_Beam(Simdat&);
	void	FileRead_Path(vector<path_seg>&, Simdat&);

	// Reads in fine tuned simulation parameters
	void	FileRead_Domain(vector<path_seg>&, Simdat&);
	void	FileRead_Settings(Simdat&);

	// Reads in utility files (miscellaneous stuff for research)
	void	FileRead_Points(Simdat&);
	void	FileRead_ParBeams(Simdat&);
	void	FileRead_InfBeams(Simdat&);

	// Initializes the positions of all points
	void	SetPoints(Point * const ptv, Simdat&);
}