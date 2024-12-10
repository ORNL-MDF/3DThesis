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
        using Stork::Structs::RDF_Header;
        using Stork::Structs::RDF_Data;
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
        sim.print = false;

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

        // Single pass to count and populate events
        size_t numEvents = 0;
        const size_t numPoints = sim.domain.pnum;

        // First, count total events
        for (size_t n = 0; n < numPoints; n++) {
            if (grid.get_output_flag(n)) {
                numEvents += grid.get_RDF_tm(n).size();
            }
        }

        // Make structure with precise size
        RDF_Dual<FloatType> RDF;
        RDF.template MakeViews<host_space>(numEvents);
        RDF_Header<FloatType, host_space>& header = RDF.host_header;
        RDF_Data<FloatType, host_space>& data = RDF.host_data;

        // TODO::MPI DECOMPOSITION AFFECTS GLOBAL i0, etc.
        // Set up header
        header.global_i0() = static_cast<uint32_t>(0);
        header.global_j0() = static_cast<uint32_t>(0);
        header.global_k0() = static_cast<uint32_t>(0);
        header.local_inum() = static_cast<uint32_t>(sim.domain.xnum);
        header.local_jnum() = static_cast<uint32_t>(sim.domain.ynum);
        header.local_knum() = static_cast<uint32_t>(sim.domain.znum);
        header.global_x0() = static_cast<FloatType>(sim.domain.xmin);
        header.global_y0() = static_cast<FloatType>(sim.domain.ymin);
        header.global_z0() = static_cast<FloatType>(sim.domain.zmin);
        header.gridResolution() = static_cast<FloatType>(sim.domain.xres);

        // Single pass population with index tracking
        size_t currentEventIndex = 0;
        for (size_t n = 0; n < numPoints; n++) {
            if (grid.get_output_flag(n)) {
                const auto& rdf_tm = grid.get_RDF_tm(n);
                const auto& rdf_tl = grid.get_RDF_tl(n);
                const auto& rdf_cr = grid.get_RDF_cr(n);
                const size_t eSize = rdf_tm.size();
                const uint32_t p = RDF.template ijk_to_p<host_space>(grid.get_i(n), grid.get_j(n), grid.get_k(n));
                for (size_t e = 0; e < eSize; e++) {
                    data.p(currentEventIndex) = p;
                    data.tm(currentEventIndex) = rdf_tm[e];
                    data.tl(currentEventIndex) = rdf_tl[e];
                    data.cr(currentEventIndex) = rdf_cr[e];
                    currentEventIndex++;
                }
            }
        }

        // Return Kokkos data
        return RDF;
    }

    template<typename FloatType>
    Stork::Structs::SRDF_Dual<FloatType>  Coupled_SRDF(int argc, char * argv[])
    {
        // Namespace declarations
        using namespace impl;
        using Stork::Structs::SRDF_Dual;
        using Stork::Structs::SRDF_Header;
        using Stork::Structs::SRDF_Data;
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
        sim.print = false;

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

        // Get number of events
        const uint32_t numEvents = ts.size()/2;

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
        SRDF_Dual<FloatType> SRDF;
        SRDF.template MakeViews<host_space>(numEvents);
        SRDF_Header<FloatType, host_space>& header = SRDF.host_header;
        SRDF_Data<FloatType, host_space>& data = SRDF.host_data;

        // TODO::MPI DECOMPOSITION AFFECTS GLOBAL i0, etc.
        // Set up header
        header.global_i0() = static_cast<uint32_t>(0);
        header.global_j0() = static_cast<uint32_t>(0);
        header.global_k0() = static_cast<uint32_t>(0);
        header.local_inum() = static_cast<uint32_t>(sim.domain.xnum);
        header.local_jnum() = static_cast<uint32_t>(sim.domain.ynum);
        header.local_knum() = static_cast<uint32_t>(sim.domain.znum);
        header.global_x0() = static_cast<FloatType>(sim.domain.xmin);
        header.global_y0() = static_cast<FloatType>(sim.domain.ymin);
        header.global_z0() = static_cast<FloatType>(sim.domain.zmin);
        header.gridResolution() = static_cast<FloatType>(sim.domain.xres);
        header.T_critical() = static_cast<FloatType>(sim.material.T_liq);

        // Make unmanaged views
        using uint32_hostView = Kokkos::View<uint32_t*, layout, host_space>;
        using floatType_hostView = Kokkos::View<FloatType*, layout, host_space>;
        uint32_hostView unmanagedView_cellIdxs(idxs.data(), idxs.size());
        floatType_hostView unmanagedView_t(ts_data, ts.size());
        floatType_hostView unmanagedView_T(Ts_data, Ts.size());

        // Convert i,j,k -> p and store value
        Kokkos::parallel_for(
        "STORK - <i,j,k> to <p>",
        Kokkos::RangePolicy<host_space>(0, numEvents),
        KOKKOS_LAMBDA(const uint32_t n)
            {   
                const uint32_t& i = unmanagedView_cellIdxs(3*n+0);
                const uint32_t& j = unmanagedView_cellIdxs(3*n+1);
                const uint32_t& k = unmanagedView_cellIdxs(3*n+2);
                data.p(n) = SRDF.template ijk_to_p<host_space>(i,j,k);
            }
        );

        // Copy rest to to managed views
        Kokkos::deep_copy(data.times_view, unmanagedView_t);
        Kokkos::deep_copy(data.thermals_view, unmanagedView_T);

        // Return Kokkos data
        return SRDF;
    }

}