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
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "DataStructs.h"

using std::vector;
using std::string;
using std::stringstream;
using std::to_string;

namespace Init {
	// Are 1st level Keywords in the file
	void	Keywords_Lv1(vector<string>&, vector<int>&, const string&);

	// Read for Keywords
	void	Keywords_Lv2(vector<string>&, vector<vector<string>>&, vector<vector<string>>&, const string&);
	void	Keywords_Lv2(vector<string>&, vector<vector<string>>&, vector<vector<double>>&, const string&);
	
	// Helps Set Values and Display Errors
	void	SetValues(string&, string, string, string, int);
	void	SetValues(int&,  string, int, string, int);
	void	SetValues(bool&, string, bool, string, int);
	void	SetValues(double&, string, double, string, int);

	// Gets file names from the input file
	void	GetFileNames(FileNames&, const string&);
	void	MakeDataDirectory(const string&);
	void	ReadSimParams(Simdat&);
	

	// Reads in simulation type information
	void	FileRead_Mode(Simdat&, const string&);
	void	FileRead_Mode_Snapshot(Simdat&, const string&);
	void	FileRead_Mode_Solidification(Simdat&, const string&);

	// Reads in all necessary simulation parameters
	void	FileRead_Material(Material&, const string&);
	void	FileRead_Beam(Beam&, const string&);
	void	FileRead_Path(vector<path_seg>&, const string&);

	// Reads in fine tuned simulation parameters
	void	FileRead_Domain(Domain&, const string&);
	void	FileRead_Output(Output&, const string&);
	void	FileRead_Settings(Settings&, const string&);

	// Reads in utility files (miscellaneous stuff for research)
	void	FileRead_Points(Domain&, const string&);

	// Set thermal diffusivity
	void	SetDiffusivity(Material&);
	// Set beam efficiency
	void	SetBeamPower(Beam&);
	// Initializes the positions of all points
	void	SetDomainParams(Domain&);

	// For reading in multiple beams
	void	FileRead_Beams(vector<Beam>&, const string&);
	// For reading in multiple paths
	void	FileRead_Paths(vector<vector<path_seg>>&, const string&);
}