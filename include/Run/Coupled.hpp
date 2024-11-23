#pragma once

#include <Stork_Core.hpp>

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
#include "MpiStructs.h"
#include <mpi.h>
#endif

namespace Thesis::Run{

    template<typename FloatType>
    Stork::Structs::RDF_Dual<FloatType> Coupled_RDF(int argc, char * argv[])
    {
        // Namespace declarations
        using namespace impl;
        using Stork::Structs::RDF_Dual;
        using Stork::Structs::RDF;
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

            // Throw errors if not the correct modes
        if (sim.param.mode!="Solidification"){
            throw std::runtime_error("Simulation Type Must Be: 'Solidification'");
        }
        if (!sim.output.RDF){
            throw std::runtime_error("'RDF' Output Must Be Turned On");
        }
        if (sim.param.mode=="Stork"){
            throw std::runtime_error("Function 'Coupled_RDF' does not work with the Stork Tracking mode.");
        }

        // Run the Simulation
        Modes::Simulate(grid, sim);

        // Collect results
        vector<FloatType> RDF_data;
        const size_t numPoints = sim.domain.pnum;
        #pragma omp parallel num_threads(sim.settings.thnum)
        {
            vector<FloatType> th_data;
            #pragma omp for schedule(static)
            for (int p = 0; p < numPoints; p++) {
                if (grid.get_output_flag(p)) {
                    const size_t pointNumEvents = grid.get_RDF_tm(p).size();
                    for (size_t e = 0; e < pointNumEvents; e++) {
                        th_data.push_back(grid.get_x(p));
                        th_data.push_back(grid.get_y(p));
                        th_data.push_back(grid.get_z(p));
                        th_data.push_back(grid.get_RDF_tm(p)[e]);
                        th_data.push_back(grid.get_RDF_tl(p)[e]);
                        th_data.push_back(grid.get_RDF_cr(p)[e]);
                    }
                }
            }
            #pragma omp critical
            {
                RDF_data.insert(RDF_data.end(), th_data.begin(), th_data.end());
            }
        }
        
        // Make views
        RDF_Dual<FloatType> DualRDF;           
        DualRDF.SetSizes(RDF_data.size() / 6);
        DualRDF.template MakeViews<host_space>();     

        // Create a temporary unmanaged Kokkos::View from RDF_data
        using floatType_hostView = Kokkos::View<FloatType*, layout, host_space>;
        floatType_hostView unmanagedView_RDF(RDF_data.data(), RDF_data.size());

        // Perform the copy into the managed Kokkos view
        Kokkos::deep_copy(DualRDF.hostRDF, unmanagedView_RDF);

        // Return Kokkos data
        return DualRDF;
    }

    template<typename FloatType>
    Stork::Structs::SRDF_Dual<FloatType>  Coupled_SRDF(int argc, char * argv[])
    {
        // Namespace declarations
        using namespace impl;
        using Stork::Structs::SRDF_Dual;
        using Stork::Structs::SRDF;
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

        // Throw errors if not the correct modes
        if (sim.param.mode!="Solidification"){
            throw std::runtime_error("Simulation Type Must Be: 'Solidification'");
        }
        if (sim.param.tracking!="Stork"){
            throw std::runtime_error("Function 'Coupled_SRDF' only works with the Stork Tracking mode.");
        }

        // Run the Simulation
        Modes::Simulate(grid, sim);

        // Set references to data
        vector<uint32_t>& idxs = grid.get_RRDF_idxs();
        vector<double>& ts = grid.get_RRDF_ts();
        vector<double>& Ts = grid.get_RRDF_Ts();

        // Create a pointer to hold either the original or converted data
        FloatType* ts_data = nullptr;
        FloatType* Ts_data = nullptr;

        // Temporary storage for conversion (only created if needed)
        vector<FloatType> ts_converted;
        vector<FloatType> Ts_converted;

        if constexpr (std::is_same<FloatType, double>::value) 
        {
            // If FloatType is double, use the original data directly
            ts_data = ts.data();
            Ts_data = Ts.data();
        } 
        else 
        {
            // Otherwise, create temporary storage and convert
            ts_converted.resize(ts.size());
            Ts_converted.resize(Ts.size());

            using std::transform;
            transform(ts.begin(), ts.end(), ts_converted.begin(), [](double value) { return static_cast<FloatType>(value); });
            transform(Ts.begin(), Ts.end(), Ts_converted.begin(), [](double value) { return static_cast<FloatType>(value); });

            ts_data = ts_converted.data();
            Ts_data = Ts_converted.data();
        }

        // Initialize Kokoks views
        SRDF_Dual<FloatType> DualSRDF;
        SRDF<FloatType, host_space>& hostSRDF = DualSRDF.hostSRDF;
        DualSRDF.size = ts.size()/2;
        DualSRDF.Init_Host_Views_Header();
        DualSRDF.Init_Host_Views_Data();

        // Header extents
        hostSRDF.i_num() = static_cast<uint32_t>(sim.domain.xnum);
        hostSRDF.j_num() = static_cast<uint32_t>(sim.domain.ynum);
        hostSRDF.k_num() = static_cast<uint32_t>(sim.domain.znum);

        // Header floats
        hostSRDF.x0() = static_cast<FloatType>(sim.domain.xmin);
        hostSRDF.y0() = static_cast<FloatType>(sim.domain.ymin);
        hostSRDF.z0() = static_cast<FloatType>(sim.domain.zmin);
        hostSRDF.res() = static_cast<FloatType>(sim.domain.xres);
        hostSRDF.T_liq() = static_cast<FloatType>(sim.material.T_liq);

        // Make unmanaged views
        using uint32_hostView = Kokkos::View<uint32_t*, layout, host_space>;
        using floatType_hostView = Kokkos::View<FloatType*, layout, host_space>;
        uint32_hostView unmanagedView_cellNum(idxs.data(), hostSRDF.cellNum_view.size());
        floatType_hostView unmanagedView_t(ts_data, ts.size());
        floatType_hostView unmanagedView_T(Ts_data, Ts.size());

        // Copy to managed views
        Kokkos::deep_copy(hostSRDF.cellNum_view,unmanagedView_cellNum);
        Kokkos::deep_copy(hostSRDF.times_view,unmanagedView_t);
        Kokkos::deep_copy(hostSRDF.thermals_view,unmanagedView_T);

        // Return Kokkos data
        return DualSRDF;
    }

}