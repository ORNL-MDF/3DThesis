#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <mpi.h>

#include "DataStructs.h"
#include "Init.h"
#include "Util.h"
#include "Run.h"
#include "Out.h"

#include "Grid.h"
#include "Mpi_structs.hpp"

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

int main(int argc, char * argv[]) {	
	// Initialize MPI
    MPI_Init(&argc, &argv);
    {
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

		// Initialize MPI
		ThesisMPI mpi(MPI_COMM_WORLD);

		// Make local bounds
		mpi.makeLocalBounds(sim);

		// Initialize grid
		Grid grid(sim); 

		// MISC INIT STUFF
		sim.util.approxEndTime = sim.util.allScansEndTime;

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

		if (sim.param.tracking!="Stork"){
			if (sim.param.mode=="Solidification"){ 
			grid.Output(sim, "Solidification.Final." + mpi.name); 
			}
			if (sim.output.T_hist) { 
				grid.Output_T_hist(sim, "T.hist." + mpi.name); 
			}
			if (sim.output.RDF) { 
				grid.Output_RDF(sim, "RDF.Final." + mpi.name); 
			}
		}
		else{
			// Output RRDF
			grid.Output_RRDF_csv(sim, "RRDF." + mpi.name);
			//grid.Output_RRDF_bin(sim, "RRDF." + mpi.name);
		}
		
		
		// Output output time
		auto stop_out = high_resolution_clock::now();
		std::cout << "Output time (s): " << (duration<double, std::milli>(stop_out - start_out).count()) / 1000.0 << "\n\n";
	}
	// Finalize MPI
    MPI_Finalize();

	return 0;
}
