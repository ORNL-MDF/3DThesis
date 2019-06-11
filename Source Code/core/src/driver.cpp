//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright 2019 UT-Battelle, LLC
*
* All Rights Reserved
*
* Authors: Benjamin Stump <stumpbc@ornl.gov>, Alex Plotkowski, James Ferguson, Kevin Sisco
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

#include "driver.h"
#include <iostream>
#include <vector>
#include <thread>

using namespace SAHTM;

void Driver::run(std::string sim_file, std::string path_file)
{
    Simdat sim(sim_file);
    util::MakeDirs(sim);
    //Create a vector of points
    std::vector<Point> ptv(sim.pnum);

    util::SetPositions(ptv, sim);

    for (int p = 0; p < sim.pnum; p++){ ptv[p].Initialize(sim); }

    sim.a = sim.kon / (sim.rho * sim.cps);	//Calculate the thermal diffusivity

    // Auto-Dtau Stuff
    sim.nond_dt = sim.ax*sim.ax / sim.a;
    sim.xres = 1e+10; //Default to large number
    sim.yres = 1e+10; //Default to large number
    sim.zres = 1e+10; //Default to large number
    if (sim.imax != 1) { sim.xres = (sim.xmax - sim.xmin) / (sim.imax - 1); }
    if (sim.jmax != 1) { sim.yres = (sim.ymax - sim.ymin) / (sim.jmax - 1); }
    if (sim.kmax != 1) { sim.zres = (sim.zmax - sim.zmin) / (sim.kmax - 1); }

    //Read in path information
    sim.max_speed = 0.0001;
    std::vector<PathSeg> segv = util::ReadPath(path_file, sim);

    // Better Search Algorythm Setup
    std::vector<int> seg_num = util::StartSeg(segv, sim);
    int seg_sum_size = (int)seg_num.size();

    //*****************************************************************
    //*****************************************************************
    //Time loop
    int tflag = 1, itert = 0, endflag = 1;

    double t2 = 0.0;
    std::vector<int> liq_pts;

    while (tflag){
        //std::cout << "Time step: " << itert << "\t\t";
        //std::cout << "% of Path: " << 100*t2 / segv[segv.size()-1].seg_time << "%" << "\n";
		if(usingProgress)
			progress_->progress(t2 / segv[segv.size()-1].seg_time);
		else
		{
			std::cout << "Time step: " << itert << "\t\t";
			std::cout << "% of Path: " << 100*t2 / segv[segv.size()-1].seg_time << "%" << "\n";
		}
        //Does something else just for getting integration info
        if (itert < seg_sum_size) {
            sim.start_seg = seg_num[itert];
        }
        else {
            sim.start_seg = seg_num[seg_sum_size-1];
        }
        if (!sim.start_seg) {
            sim.start_seg = 1;
        }
        std::vector<int_seg> isegv;
        util::GetIntegrationInfo(t2, segv, sim, isegv,0);

        //Mode 1 calculates all information at all points at all time steps
        if (sim.mode == 1){

            #pragma omp parallel for num_threads(sim.thnum) schedule(dynamic)
            for (int p = 0; p < sim.pnum; p++){
                ptv[p].Tlast = ptv[p].T;
                ptv[p].Temp_Calc_Pre_Path(t2, isegv, sim, 1, 0);
            }


            #pragma omp parallel for num_threads(sim.thnum) schedule(dynamic)
            for (int p = 0; p < sim.pnum; p++){
                if ((ptv[p].T < sim.Tsol) && (ptv[p].Tlast > sim.Tsol)){ ptv[p].CalcGo_2(t2, segv, sim); }
            }
        }

        //Mode 2 searches for only the melt melt pool
        if (sim.mode == 2){
            //Set T_calc_flag at all points to indicate that they have not yet been calculated
            util::SetFlags(ptv);

            //See which points from previous time step have solidified and do appropriate calculations
            std::vector<int> last_liq_pts = liq_pts;
            liq_pts.clear();
            liq_pts = util::CalcSolidifiedPoints(ptv, last_liq_pts, t2, segv, sim, isegv);

            //Start search from points that are known to be liquid
            std::vector<int> test_pts = liq_pts;

            //Trace beam path between time steps and add relevant points to the test vector
            if (endflag) {
                if (t2 > segv[sim.np - 1].seg_time) {
                    endflag = 0;
                }
                else if (itert) {
                    util::BeamTrace(seg_num,itert, ptv, segv, sim, test_pts);
                }
            }


            int melt_pool_flag = 1;		//Flag, gets set to zero when there are no more new points to check
            while (melt_pool_flag){		//Iterative loop expanding from previously identified points to find melt pool boundary
                std::vector<int> test_tmp;
                std::vector<int> test_check;

                //identify neighbors of liquid points
                util::NeighborCheck(ptv, sim, test_pts, test_check);

                //Perform temperature calculations of necessary points
                util::VectorTempCalc(t2, ptv, sim, test_check, isegv);

                //Add points above liquidus into melt pool vector
                for (size_t it = 0; it < test_check.size(); it++){
                    if (ptv[test_check[it]].T > sim.Tsol){
                        test_tmp.push_back(test_check[it]);
                        liq_pts.push_back(test_check[it]);
                    }
                }

                test_pts.clear();
                test_pts = test_tmp;
                if (test_pts.size() == 0){ melt_pool_flag = 0; }
            }
        }

        //If we are writing an output file create a thread to do it.
        if(itert && itert % sim.out_freq == 0)
        {
            std::thread t1(util::Output, ptv, itert, 0, sim);
            t1.detach(); //detech thread so rest of program continues.
        }

        t2 += sim.dt;							//Increment time
        itert++;								//Increment time step counter
        tflag = util::SimFinished(ptv, t2, sim, segv);	//Check if simulation is finished
    }
    //*****************************************************************
    //*****************************************************************

    //Output final file
    util::Output(ptv, itert, 1, sim);
}
