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

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

#include "DataStructs.h"
#include "Init.h"
#include "Util.h"
#include "Run.h"
#include "Out.h"

#include "Grid.h"

using std::vector;
using std::string;

using std::chrono::high_resolution_clock;
using std::chrono::duration;

int main(int argc, char * argv[]) {	

	// Start initialization clock
	auto start_in = high_resolution_clock::now();

	//Command line input to program is the file that contains the simulation file names
	string	inputFile;
	if (argc <= 1) { inputFile = "TestInputs/ParamInput.txt"; }
	else { inputFile = argv[1]; }

	// Initialize struct for simulation parameters
	Simdat sim;

	// Get names of intput files
	Init::GetFileNames(sim.files, inputFile);	
	
	// Read input files and set simulation parameters
	Init::ReadSimParams(sim);

	// Initialize grid
	Grid grid(sim); 

	//// MISC INIT STUFF ////
	//if (sim.settings.use_PINT) { Util::EstimateEndTime(sim, segv); }
	//else { sim.util.approxEndTime = sim.util.scanEndTime; }
	sim.util.approxEndTime = sim.util.allScansEndTime;
	////////////////////////

	// Output initialization time
	auto stop_in = high_resolution_clock::now();
	std::cout << "Initialization time (s): " << (duration<double, std::milli>(stop_in - start_in).count())/1000.0 << "\n\n"; //(stop_in - start_in) / double(CLOCKS_PER_SEC) << "\n\n";

	// Start simulation clock
	auto start_sim = high_resolution_clock::now();

	// Run the simulation
	Run::Simulate(grid, sim);

	// Output simulation time
	auto stop_sim = high_resolution_clock::now();
	std::cout << "Execution time (s): " << (duration<double, std::milli>(stop_sim - start_sim).count())/1000.0 << "\n\n";//(stop_sim - start_sim) / double(CLOCKS_PER_SEC) << "\n\n";

	// Start output clock
	auto start_out = high_resolution_clock::now();

	if (sim.param.mode=="Solidification" && sim.param.tracking!="Stork"){ grid.Output(sim, "Solidification.Final"); }
	if (sim.output.T_hist) { grid.Output_T_hist(sim, "T.hist"); }
	if (sim.output.RDF) { grid.Output_RDF(sim, "RDF.Final"); }
	
	// Output output time
	auto stop_out = high_resolution_clock::now();
	std::cout << "Output time (s): " << (duration<double, std::milli>(stop_out - start_out).count()) / 1000.0 << "\n\n";

	return 0;
}
