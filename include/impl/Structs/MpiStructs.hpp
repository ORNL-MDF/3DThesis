#pragma once

#include <Stork_Core.hpp>

#include <mpi.h>
#include <string>

#include "impl/Structs/DataStructs.hpp"
#include "impl/Calc/Util.hpp"

using std::string;

namespace Thesis::impl
{

    template<typename FloatType>
    class ThesisMPI{
    private:
        // Communication Stuff
        MPI_Comm comm;

        // MPI info
        int rank;
        int nproc;

        // Locality info
        int i_min, i_max;
        int j_min, j_max;

    public:
        ThesisMPI(MPI_Comm inComm){
            // Set Comm information
            comm = inComm;

            // Get MPI Info
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &nproc);

            // Set mpi name
            name = Util::ZeroPadNumber(rank, 1+nproc/10);

            // TODO::CUSTOM DECOMPOSITION
            // Get dims for 2D decomposition
            MPI_Dims_create(nproc, 2, dims);

            // Create coms
            int periods[2] = {0,0}; // non-periodic boundaries
            MPI_Comm cart_comm;
            MPI_Cart_create(comm, 2, dims, periods, 0, &cart_comm); 

            // Set coords to coordinate of local process
            MPI_Cart_coords(cart_comm, rank, 2, coords);  
        }

        // Dimensionality of decomposition
        int dims[2] = {0,0};    // Dimensionality of decomposition
        int coords[2] = {0,0};   // Own coordinates

        // Stork Header
        Stork::Structs::RegularGrid_Header<FloatType, host_space> header;

        // Name addition
        string name;

        int size() { return nproc; }

        void setPrint(Simdat& sim) {
            // Update local rank printing
            sim.print = rank == 0;
            sim.mpi = size() > 1;
        }

        // Make x-y bounds for local domain
        void makeLocalBounds(Simdat& sim){
            // Global Decomposition
            const int I = dims[0];
            const int J = dims[1];
            
            // Local Decomposition Number
            const int i = coords[0];
            const int j = coords[1];

            // Find local domain bounds (no overlap for else)
            if (sim.param.tracking=="Stork" || sim.settings.mpi_overlap){ 
                // Find local domain bounds (overlap for stork)
                // Note:: Doesn't just allow for all points to be filled but overlaps the interpolated points by (so 1 to 1 data handoff is valid)
                // i_min = ((sim.domain.xnum-1)*i)/I;
                // i_max = ((sim.domain.xnum-1)*(i+1))/I + (i!=(I-1));

                // j_min = ((sim.domain.ynum-1)*j)/J;
                // j_max = ((sim.domain.ynum-1)*(j+1))/J + (j!=(J-1));
                // TODO::DEBUG
                sim.param.radiusCheck = 3;
                if (rank==0){
                    i_min = 0;
                    i_max = (sim.domain.xnum-1)-50;

                    j_min = 0;
                    j_max = (sim.domain.ynum-1);
                }
                else{
                    i_min = (sim.domain.xnum-1)-51;
                    i_max = (sim.domain.xnum-1);

                    j_min = 0;
                    j_max = (sim.domain.ynum-1);
                }
            }
            else{
                // Find local domain bounds (no overlap for else)
                i_min = ((sim.domain.xnum-1)*i)/I + (i!=0);
                i_max = ((sim.domain.xnum-1)*(i+1))/I;

                j_min = ((sim.domain.ynum-1)*j)/J + (j!=0);
                j_max = ((sim.domain.ynum-1)*(j+1))/J;
            }

            // Set header info
            header.global_i0() = i_min;
            header.global_j0() = j_min;
            header.global_x0() = sim.domain.xmin;
            header.global_y0() = sim.domain.ymin;
            header.local_inum() = 1+(i_max-i_min);
            header.local_jnum() = 1+(j_max-j_min);

            // Find local domain bounds (no)
            const double xmin = sim.domain.xmin + sim.domain.xres*i_min;
            const double xmax = sim.domain.xmin + sim.domain.xres*i_max;
            
            const double ymin = sim.domain.ymin + sim.domain.yres*j_min;
            const double ymax = sim.domain.ymin + sim.domain.yres*j_max;

            // Now adjust local domain
            sim.domain.xmin = xmin; sim.domain.xmax = xmax;
            sim.domain.xnum = 1 + int(0.5 + (sim.domain.xmax - sim.domain.xmin) / sim.domain.xres);
            sim.domain.xmax = sim.domain.xmin + (sim.domain.xnum - 1) * sim.domain.xres;
            
            sim.domain.ymin = ymin; sim.domain.ymax = ymax;
            sim.domain.ynum = 1 + int(0.5 + (sim.domain.ymax - sim.domain.ymin) / sim.domain.yres);
            sim.domain.ymax = sim.domain.ymin + (sim.domain.ynum - 1) * sim.domain.yres;

            sim.domain.pnum = sim.domain.xnum*sim.domain.ynum*sim.domain.znum;
        }
    };
}
