//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) - Transformer Resilience and Advanced Components (TRAC) Program.

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
