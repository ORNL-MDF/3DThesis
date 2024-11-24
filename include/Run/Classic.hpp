#pragma once

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

#include "impl/Modes/ClassicTypes.hpp"

#include "impl/Structs/DataStructs.hpp"
#include "impl/IO/Init.hpp"
#include "impl/Calc/Util.hpp"
#include "impl/IO/Out.hpp"
#include "impl/Structs/Grid.hpp"
#include "ThesisConfig.h"

#ifdef Thesis_ENABLE_MPI
#include "impl/Structs/MpiStructs.h"
#include <mpi.h>
#endif

namespace Thesis::Run{

    inline void Classic(int argc, char * argv[])
    {
        // Namespace declarations
        using namespace impl;
        using std::vector;
        using std::string;
        using std::chrono::high_resolution_clock;
        using std::chrono::duration;
        
        // Start initialization clock
        auto start_in = high_resolution_clock::now();

        //Command line input to program is the file that contains the simulation file names
        if (argc <= 1) { throw std::runtime_error( "Input file argument required: ./3DThesis ParamInput.txt" ); }
            string inputFile = argv[1];

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
        Modes::Simulate(grid, sim);

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
            if (sim.param.tracking=="Stork"){
                // Output RRDF
                //grid.Output_RRDF_csv(sim, "RRDF" + rank_name);
                grid.Output_RRDF_bin(sim, "RRDF" + rank_name);
            }	
            else{
                grid.Output(sim, "Solidification.Final" + rank_name); 
            }
        }
        if (sim.output.T_hist) { 
            grid.Output_T_hist(sim, "T.hist" + rank_name); 
        }
        if (sim.output.RDF) { 
            grid.Output_RDF(sim, "RDF.Final" + rank_name); 
        }
										
        
        // Output output time
        auto stop_out = high_resolution_clock::now();
        if (sim.print)
            std::cout << "Output time (s): " << (duration<double, std::milli>(stop_out - start_out).count()) / 1000.0 << "\n\n";
    }

}