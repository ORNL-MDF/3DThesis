//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright(c) 2019 Benjamin Stump, Alex Plotkowski, James Ferguson, Kevin Sisco <stumpbc@ornl.gov>
* All rights reserved.
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

#include "util.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <io.h>

#ifdef _WIN32 // conditional needed to prevent error when compiling for unix systems
#include <direct.h>
#else // conditional for unix systems; should work *most of the time
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

using namespace SAHTM;

// Takes in an integer and outputs a string with that integer plus leading zeros
std::string ZeroPadNumber(int num)
{
    std::ostringstream ss;
    ss << std::setw(7) << std::setfill('0') << num;
    return ss.str();
}

void util::MakeDirs(Simdat& sim)
{
    //Make Data directory if it does not already exist
    std::string strPath = sim.out_path + "/Data";
    #if defined(_WIN32)
        _mkdir(strPath.c_str());
    #else
        mkdir(strPath.c_str(), 0777); // notice that 777 is different than 0777
    #endif
    //Make SubData directory if it does not already exist
    std::string strSubPath = sim.out_path + "/Data/" +sim.sim_name+ "(Data)";
    #if defined(_WIN32)
        _mkdir(strSubPath.c_str());
    #else
        mkdir(strSubPath.c_str(), 0777); // notice that 777 is different than 0777
    #endif
}

// Sets positions of points on a regular three dimensional grid
void util::SetPositions(std::vector<Point>& ptv, Simdat& sim){
    double xp, yp, zp;
    int p = 0;
    for (int i = 0; i < sim.imax; i++){
        if (sim.imax == 1){ xp = sim.xmax; }
        else{ xp = sim.xmin + ((float)i * ((sim.xmax - sim.xmin) / ((float)sim.imax - 1))); }
        for (int j = 0; j < sim.jmax; j++){
            if (sim.jmax == 1){ yp = sim.ymax; }
            else{ yp = sim.ymin + ((float)j * ((sim.ymax - sim.ymin) / ((float)sim.jmax - 1))); }
            for (int k = 0; k < sim.kmax; k++){
                if (sim.kmax == 1){ zp = sim.zmax; }
                else{ zp = sim.zmin + ((float)k * ((sim.zmax - sim.zmin) / ((float)sim.kmax - 1))); }
                ptv[p].i = i;
                ptv[p].j = j;
                ptv[p].k = k;
                ptv[p].x = xp;
                ptv[p].y = yp;
                ptv[p].z = zp;
                p++;
            }
        }
    }
}

//Map three-dimensional points onto one-dimensional vector
int util::GetP(int i, int j, int k, Simdat& sim){
    return i * (sim.kmax * sim.jmax) + j * sim.kmax + k;
}

//Jamie Function for start location at each dt
std::vector<int> util::StartSeg(std::vector<PathSeg>& segv, Simdat& sim) {
    std::vector<int> segnum;
    double f_time = segv[sim.np - 1].seg_time;
    int seg_n = 0;
    double a = 0.0;
    while (a <= f_time) {
        while ((segv[seg_n].seg_time<(a - (1e-9))) && (segv[seg_n].seg_time<(a + (1e-9)))) { // Toleratnce of 1ns
            seg_n++;
        }
        segnum.push_back(seg_n);
        a += sim.dt;
    }
    return segnum;
}

//Function to read the path file and convert to a vector of path segments
std::vector<PathSeg>  util::ReadPath(std::string path_file, Simdat& sim){
    std::vector<PathSeg> segv;
    //int mode;
    //double xp, yp, zp, pm, par;
    double convert = 0.001;	//Currently hard coding a conversion from mm to m

    std::ifstream pathfile;
	pathfile.exceptions(ifstream::badbit);
    std::string line;
    try{
        pathfile.open(path_file, std::ios::in);

        //Set initial position as 0, 0, 0
        int smode = 1;
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double sqmod = 0.0;
        double sparam = 0.0;
        segv.push_back(PathSeg(smode,sx,sy,sz,sqmod,sparam));

        //Read in path information from file
        while (getline(pathfile, line))
        {
            pathfile >> smode >> sx >> sy >> sz >> sqmod >> sparam;
            sx *= convert;
            sy *= convert;
            sz *= convert;
            segv.push_back(PathSeg(smode,sx,sy,sz,sqmod,sparam));
        }
    }
    catch (const std::ifstream::failure&) {
        std::cout << "Exception opening/reading path file\n";
    }
    catch (const std::exception &e) {
        std::cout << "caught an exception " << e.what() << '\n';
    }
    pathfile.close();

    sim.np = (int)segv.size();
    std::cout << sim.np << " path segments were successfully imported \n";
    //Right now, the above always returns an exception if the last line of the file is blank,
    //but will run fine with the path segments that were correctly imported, so not really an error

    //Calculate path times
    double dt_seg, dist, dx, dy, dz;
    segv[0].seg_time = 0;
    for (int p = 1; p < sim.np; p++){
        if (segv[p].smode){	//For spot mode
            segv[p].seg_time = segv[p - 1].seg_time + segv[p].sparam;
        }
        else{						//For line mode
            dx = segv[p].sx - segv[p - 1].sx;
            dy = segv[p].sy - segv[p - 1].sy;
            dz = segv[p].sz - segv[p - 1].sz;
            dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));
            dt_seg = dist / segv[p].sparam;
            segv[p].seg_time = segv[p - 1].seg_time + dt_seg;
            if (segv[p].sparam > sim.max_speed) { sim.max_speed = segv[p].sparam; }
        }
    }
    return segv;
}

double GetRefTime(double& t2, std::vector<PathSeg>& segv, Simdat& sim, int& seg) {

    double tau = (t2 - segv[seg].seg_time);
    if (tau < 0) { tau = 0; }
    double SOI = sqrt((sim.ax*sim.ax + 12 * sim.a*tau)*log(sim.tmod) / 3); //The "sphere of influence" of the heat source

    double ref_t=0;
    double dx, dy, dz, dt_cur;
    dt_cur = segv[seg].seg_time - segv[seg - 1].seg_time;
    if (segv[seg].smode) {	//Location calculation for spot mode
        ref_t = dt_cur;
        if (ref_t==0){
            ref_t = 0.0001;
        }
    }
    else {							//Location calculation for line mode
        dx = std::fabs(segv[seg].sx - segv[seg - 1].sx);
        dy = std::fabs(segv[seg].sy - segv[seg - 1].sy);
        dz = std::fabs(segv[seg].sz - segv[seg - 1].sz);
        double ref_t_min = dt_cur/fmax(fmax(dx/sim.xres,dy/sim.yres), dz/sim.zres); //Makes sure to, at a minimum, integrate twice between each point (to maintain symmetry)
        double ref_t_tmod = dt_cur*(SOI/fmax(fmax(dx,dy),dz));
        ref_t = fmax(ref_t_min, ref_t_tmod);
    }
    if (ref_t < 0) {ref_t = 0;}
    //if (ref_t > sim.dt) { ref_t = sim.dt; } //Put in mainly for cross sections
    return ref_t;
}

// Uses a binary search to calculate the location of the beam at
// an arbitrary time
int_seg GetBeamLoc(double time, std::vector<PathSeg>& segv, int& ref_seg){

    int flag = 1;
    double dx, dy, dz, tcur, dt_cur;
    int_seg current_seg;
    //int ref_seg = sim.start_seg;
    int seg = 0;

    //Skips binary search, starts at known segment
    seg = ref_seg;
    while (flag) {
        if (time < segv[ref_seg - 1].seg_time){
            ref_seg--;
        }
        else {
            seg = ref_seg;
            flag = 0;
        }
    }

    if (segv[seg].smode) {	//Location calculation for spot mode
        current_seg.xb = segv[seg].sx;
        current_seg.yb = segv[seg].sy;
        current_seg.zb = segv[seg].sz;
    }
    else {							//Location calculation for line mode
        dx = segv[seg].sx - segv[seg - 1].sx;
        dy = segv[seg].sy - segv[seg - 1].sy;
        dz = segv[seg].sz - segv[seg - 1].sz;
        tcur = time - segv[seg - 1].seg_time;
        dt_cur = segv[seg].seg_time - segv[seg - 1].seg_time;
        current_seg.xb = segv[seg - 1].sx + (tcur / dt_cur)*dx;
        current_seg.yb = segv[seg - 1].sy + (tcur / dt_cur)*dy;
        current_seg.zb = segv[seg - 1].sz + (tcur / dt_cur)*dz;
    }
    current_seg.qmod = segv[seg].sqmod;
    return current_seg;
}

// Determine if simulation has finished
int  util::SimFinished(std::vector<Point>& ptv, double t2, Simdat& sim, std::vector<PathSeg>& segv){
    int tflag = 1;
    if (t2 > segv[sim.np - 1].seg_time) {
        tflag = 0;
        size_t p_size = ptv.size();
        for (size_t p = 0; p < p_size; p++) {
            if (ptv[p].T> sim.Tsol) { tflag = 1; }
        }
    }
    return tflag;
}

// Claculate t1, the lower bound of numerical integration
double  util::t0calc(double t2, Simdat& sim){
    double t0;
    double safe_t = sim.nond_dt/12.0*(pow(sim.t_hist, (2.0 / 3.0)) - 1);
    if (t2 < safe_t){ t0 = 0; }
    else { t0 = t2-safe_t; }
    return t0;
}

// Writes output files

void util::Output(std::vector<Point> ptv, int itert, int out, Simdat sim){

    if (itert && (itert % sim.out_freq == 0 || out == 1)){
        std::ofstream datafile;
        datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        std::string time_step = ZeroPadNumber(itert);
        std::string out_file = sim.out_path + "/Data/" + sim.sim_name + "(Data)/" + sim.sim_name + "." + time_step + ".csv";
        try{
            datafile.open(out_file);
            datafile << "x,y,z,T,G,V\n";		//File header
            //Write data file
            for (int p = 0; p < sim.pnum; p++){
                if (ptv[p].output_flag || (sim.mode == 1)){
                    datafile << ptv[p].Print();
//                    datafile << ptv[p].x << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
//                    datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << "\n";
                }
            }
        }
        catch (const std::ofstream::failure&){ std::cout << "Exception writing data file, check that output directory exists\n"; }
        datafile.close();
    }

}

//Set T_calc_flag at all points to indicate that they have not yet been calculated
void util::SetFlags(std::vector<Point>& ptv){

    for (size_t p = 0; p < ptv.size(); p++){
        ptv[p].T_calc_flag = false;
        ptv[p].s_flag = false;
    }
}

//This function call temperature calculation for all points that were liquid in the previous time step
//Solidification parameters are calcualted for any points that have solidified since then
//The points that remain liquid are returned to main
std::vector<int> util::CalcSolidifiedPoints(std::vector<Point>& ptv, std::vector<int>& last_liq_pts, double t2, std::vector<PathSeg>& segv, Simdat& sim, std::vector<int_seg>& isegv){

    for (size_t p = 0; p < ptv.size(); p++){ ptv[p].Tlast = ptv[p].T; }
    #pragma omp parallel for num_threads(sim.thnum) schedule(dynamic)
    for (int it = 0; it < last_liq_pts.size(); it++){
        ptv[last_liq_pts[it]].Temp_Calc_Pre_Path(t2, isegv, sim, 1, 0);
        if (ptv[last_liq_pts[it]].T < sim.Tsol){ ptv[last_liq_pts[it]].CalcGo_2(t2, segv, sim); }	//Check liquid points to see if they have solidified
    }

    std::vector<int> liq_pts;

    //Add any points that are still liquid from the previous time step
    for (size_t it = 0; it < last_liq_pts.size(); it++){ if (ptv[last_liq_pts[it]].T > sim.Tsol) liq_pts.push_back(last_liq_pts[it]); }

    return liq_pts;
}

//Traces the path of the beam between time steps and adds corresponding grid points to vector for testing
void util::BeamTrace(std::vector<int>& seg_num, int& itert, std::vector<Point>& ptv, std::vector<PathSeg>& segv, Simdat& sim, std::vector<int>& test_pts){
    int seg_temp_now = seg_num[itert];
    int seg_temp_prev = seg_num[itert-1]; //make the last start_seg
    int x_grid_num=0, y_grid_num=0, z_grid_num=0;
    int x_flat=0, y_flat=0, z_flat=0;
    int point_num = 0;
    //int x_grid_floor, x_grid_ceil, y_grid_floor, y_grid_ceil, z_grid_floor, z_grid_ceil;
    for (int i = seg_temp_prev-1; i <= seg_temp_now; i++) {
        if (i < 0) { i = 0; }
        if (sim.imax - 1) { x_grid_num = (int)((segv[i].sx - sim.xmin) / sim.xres); x_flat=1; }
        if (sim.jmax - 1) { y_grid_num = (int)((segv[i].sy - sim.ymin) / sim.yres); y_flat=1; }
        if (sim.kmax - 1) { z_grid_num = (int)((segv[i].sz - sim.zmin) / sim.zres); z_flat=1; }

        if (x_grid_num < 0) { x_grid_num = 0; x_flat = 0;}
        else if (x_grid_num >= (sim.imax - 1)) { x_grid_num = sim.imax - 1; x_flat = 0;}
        if (y_grid_num < 0) { y_grid_num = 0; y_flat = 0;}
        else if (y_grid_num >= (sim.jmax - 1)) { y_grid_num = sim.jmax - 1; y_flat = 0;}
        if (z_grid_num < 0) { z_grid_num = 0; z_flat = 0;}
        else if (z_grid_num >= (sim.kmax - 1)) { z_grid_num = sim.kmax - 1; z_flat = 0;}

        std::vector<int> point_nums;
        for (int a = 0; a <= z_flat; a++) {
            for (int b = 0; b <= y_flat; b++) {
                for (int c = 0; c <= x_flat; c++) {
                    point_num=(z_grid_num+a) + sim.kmax*(y_grid_num+b) + sim.kmax*sim.jmax*(x_grid_num+c);
                    if (!ptv[point_num].T_calc_flag && !ptv[point_num].check_flag) {
                        test_pts.push_back(point_num);
                        ptv[point_num].check_flag = true;
                    }
                }
            }
        }
    }
}


//Neighbor check adds relevant neighbors to test_check vector for later tempeature calculation
void util::NeighborCheck(std::vector<Point>& ptv, Simdat& sim, std::vector<int>& test_pts, std::vector<int>& test_check){

    for (size_t p = 0; p < ptv.size(); p++){ ptv[p].check_flag = false; }
    //Find neighbors in test_pts for checking
    for (size_t it = 0; it < test_pts.size(); it++){
        //Get i, j, k location of current point, then construct array of neighbors
        int i = ptv[test_pts[it]].i;
        int j = ptv[test_pts[it]].j;
        int k = ptv[test_pts[it]].k;

        int neighborhood = 1;  //How many points to check to each size (1 means 3x3x3 cube and since it's the center, check 26pts...2 would be 5x5x5, etc)
        int n = neighborhood;

        std::vector <int> ijkminmax;
        ijkminmax.push_back(n - i);
        ijkminmax.push_back(i+1 + n - sim.imax);
        ijkminmax.push_back(n - j);
        ijkminmax.push_back(j+1 + n - sim.jmax);
        ijkminmax.push_back(n - k);
        ijkminmax.push_back(k+1 + n - sim.kmax);
        for (size_t temp = 0; temp < ijkminmax.size(); temp++) {
            if (ijkminmax[temp] < 0) { ijkminmax[temp] = 0; }
        }

        std::vector<int> nbs;
        for (int di = -n + ijkminmax[0]; di <= (n - ijkminmax[1]); di++) {
            for (int dj = -n + ijkminmax[2]; dj <= (n - ijkminmax[3]); dj++) {
                for (int dk = -n + ijkminmax[4]; dk <= (n - ijkminmax[5]); dk++) {
                    nbs.push_back(GetP(i + di, j + dj, k + dk, sim));
                }
            }
        }

        for (size_t nb = 0; nb < nbs.size(); nb++){
            if (nbs[nb] >= 0 && nbs[nb] < sim.pnum){
                if (!ptv[nbs[nb]].T_calc_flag && !ptv[nbs[nb]].check_flag && nbs[nb]<sim.pnum){
                    ptv[nbs[nb]].check_flag = true;
                    test_check.push_back(nbs[nb]);	//Flag for temperature calculation
                }
            }
        }
    }
}

void util::VectorTempCalc(double t2, std::vector<Point>& ptv, Simdat& sim, std::vector<int>& test_check, std::vector<int_seg>& isegv){
    //Perform temperature calculations of requestd points
    #pragma omp parallel for num_threads(sim.thnum) schedule(dynamic)
    for (int it = 0; it < test_check.size(); it++){
        ptv[test_check[it]].Temp_Calc_Pre_Path(t2, isegv, sim, 1, 0);
    }
}

void util::GetIntegrationInfo(double t2, std::vector<PathSeg>& segv, Simdat& sim, std::vector<int_seg>& isegv, int is_sol) {
    int seg_temp = sim.start_seg; //Stores the starting seg to change it back

    while ((t2 > segv[seg_temp].seg_time) && (t2 < segv[segv.size() - 1].seg_time)) {seg_temp++;}
    while (t2 < segv[seg_temp - 1].seg_time) {seg_temp--;}

    Point pt_temp;

    double t0 = t0calc(t2, sim); //the time you want to integrate back too
    int maxOrder = 16; //if you want more than 16, you'd have to add them
    double stepmod = 2; // 2 is default

    double nodes[30] = {
        -0.57735027,  0.57735027,
        -0.86113631, -0.33998104,  0.33998104,  0.86113631,
        -0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
        -0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
    };
    double weights[30] = {
        1.0, 1.0,
        0.34785485, 0.65214515, 0.65214515, 0.34785485,
        0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
        0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
    };


    double ref_time;
    if (t2 > segv[sim.np - 1].seg_time) { ref_time = 1e9; }
    else { ref_time = GetRefTime(t2, segv, sim, seg_temp); }

    double ref_speed = sim.max_speed; //usually max velocity
    double nond_v = (sim.ax / sim.a)*ref_speed;
    int n = (int)ceil(nond_v / 5.0); //incorporates nondimensional velocity into the minimum timestep

    double curStep = sim.nond_dt / n;
    if (is_sol) { curStep *= sim.az; }
    double curStep_use = curStep;
    double curStep_tot = 0;
    int curOrder = maxOrder;
    int tflag = 1;
    double tt = t2;

    int_seg current_beam = GetBeamLoc(t2, segv, seg_temp); //Make 1st segment at the exact time...for instantaneous heat source additon to laplacian
    current_beam.taui = t2;
    current_beam.dtau = 0.0;
    isegv.push_back(current_beam);

    while (tflag) {
        int nflag = 1;
        int i = 1;

        while (nflag) {
            int gflag = 1;
            if (curStep > (ref_time * curOrder)) {
                curStep_use = ref_time * curOrder;
            }
            else {
                curStep_use = curStep;
            }
            double t1_temp = tt - curStep_use;
            double next_time = segv[seg_temp-1].seg_time;
            if (tt > (segv[sim.np - 1].seg_time)) {
                next_time = segv[sim.np - 1].seg_time;
            }
            if (t1_temp < next_time) {	//If we are at the end of a segment, hit the end of it and set the program to jump to the next segment next time
                t1_temp = next_time;
                if (next_time != 0.0) {
                    gflag = 0;
                }
                else {
                    nflag = 0;
                    tflag = 0;
                }
            }
            else if (t1_temp < t0) { //Only will hit if auto.dtau is turned off, then it will just end the loop  (old code and don't want to mes with it)
                t1_temp = t0;
                nflag = 0;
                tflag = 0;
            }
            if (t1_temp < segv[sim.np - 1].seg_time) {
                for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
                    double tp = 0.5 * ((tt - t1_temp)*nodes[a] + (tt + t1_temp));
                        int_seg current_beam = GetBeamLoc(tp, segv, seg_temp);
                        current_beam.taui = tp;
                        current_beam.dtau = 0.5 * (tt - t1_temp) * weights[a];
                        if (current_beam.qmod>0.0) { isegv.push_back(current_beam);}
                }
            }
            if (gflag) {	//If we are not switching segments, do the normal thing (just increment time backwards)
                tt -= curStep_use;
                curStep_tot += curStep_use;
                i++;
                if (i >= n) {
                    nflag = 0;
                }
            }
            else {	//If we are switching segments, set start time to start of next segment and increment the start segment down
                curStep_tot += (tt - next_time);
                tt = next_time;
                if (tt<segv[sim.np - 1].seg_time){
                    seg_temp--;
                    ref_time = GetRefTime(t2, segv, sim, seg_temp); //Gets reference time (dtau_min) for the next segment
                }

                gflag = 1;
            }
            if (curStep_tot >= (curStep*n)) {
                curStep *= stepmod;
                if (curOrder != 2) {
                    curOrder = (curOrder / 2);
                }
            }
        }
    }
}
