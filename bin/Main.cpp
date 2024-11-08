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
#include "ThesisConfig.h"

#ifdef Thesis_ENABLE_MPI
#include "MpiStructs.h"
#include <mpi.h>
#endif

using std::vector;
using std::string;

using std::chrono::high_resolution_clock;
using std::chrono::duration;

// Read in file
// File also has MPI parameters??
//// - Make MPI file?
//// - Overlap option
// Split file based on rank (spatially)
// Run everything else normally
// Output should say MPI information

inline void run(int argc, char * argv[])
{
	// Start initialization clock
	auto start_in = high_resolution_clock::now();

	//Command line input to program is the file that contains the simulation file names
	string	inputFile;
	if (argc <= 1) { inputFile = "TestInputs/ParamInput.txt"; }
	else { inputFile = argv[1]; }

	// Initialize struct for simulation parameters
	Simdat sim;

#ifdef Thesis_ENABLE_MPI
	// Initialize MPI
	ThesisMPI mpi(MPI_COMM_WORLD);
	mpi.setPrint(sim);
#endif

	// Get names of intput files
	Init::GetFileNames(sim.files, inputFile, sim.print);	
	
	// Read input files and set simulation parameters
	Init::ReadSimParams(sim);

#ifdef Thesis_ENABLE_MPI
	// Make local bounds and set local rank
	mpi.makeLocalBounds(sim);
#endif

	// Initialize grid
	Grid grid(sim); 

	// MISC INIT STUFF
	sim.util.approxEndTime = sim.util.allScansEndTime;

	// Output initialization time
	auto stop_in = high_resolution_clock::now();
	if (sim.print)
		std::cout << "Initialization time (s): " << (duration<double, std::milli>(stop_in - start_in).count())/1000.0 << "\n\n"; //(stop_in - start_in) / double(CLOCKS_PER_SEC) << "\n\n";

	// Start simulation clock
	auto start_sim = high_resolution_clock::now();

	// Run the simulation
	Run::Simulate(grid, sim);

	// Output simulation time
	auto stop_sim = high_resolution_clock::now();
	if (sim.print)
		std::cout << "Execution time (s): " << (duration<double, std::milli>(stop_sim - start_sim).count())/1000.0 << "\n\n";//(stop_sim - start_sim) / double(CLOCKS_PER_SEC) << "\n\n";

	// Start output clock
	auto start_out = high_resolution_clock::now();

std::string rank_name = "";
#ifdef Thesis_ENABLE_MPI
if (sim.mpi)
	rank_name = "." + mpi.name;
#endif

	if (sim.param.mode=="Solidification"){ 
		grid.Output(sim, "Solidification.Final" + rank_name); 
	}
	if (sim.output.T_hist) { 
		grid.Output_T_hist(sim, "T.hist" + rank_name); 
	}
	if (sim.output.RDF) { 
		grid.Output_RDF(sim, "RDF.Final" + rank_name); 
	}
	if (sim.param.mode=="Stork"){
		// Output RRDF
		grid.Output_RRDF_csv(sim, "RRDF" + rank_name);
		//grid.Output_RRDF_bin(sim, "RRDF" + rank_name);
	}											
	
	// Output output time
	auto stop_out = high_resolution_clock::now();
	if (sim.print)
		std::cout << "Output time (s): " << (duration<double, std::milli>(stop_out - start_out).count()) / 1000.0 << "\n\n";
}

#ifdef Thesis_ENABLE_MPI
int main(int argc, char * argv[]) {	
	// Initialize MPI
    MPI_Init(&argc, &argv);
	run(argc, argv);
	// Finalize MPI
    MPI_Finalize();

	return 0;
}
#else
int main(int argc, char * argv[]) {	
	run(argc, argv);
	return 0;
}
#endif
