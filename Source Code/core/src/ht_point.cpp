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

#include <vector>
#include <cmath>
#include <sstream>

#include "ht_point.h"
#include "util.h"

#define PI 3.14159

using namespace std;
using namespace SAHTM;

//Constructor and destructor
Point::Point(void){}
Point::~Point(void){}

// Set initial parameters
void Point::Initialize(Simdat& sim){
    T = sim.Tinit;

    //Uses very small value instead of 0 for
    //variables that will be displayed in log scale
    Gtemp = 1e-100;
    G = 1e-100;
    V = 1e-100;
    Gmag = 1e-100;
    kh = 0;
    khmax = 0;
    theta = 0;
    Gx = Gy = Gz = 0;
    output_flag = 0;
}

void Point::CalcGo_2(double t2, vector<PathSeg>& segv, Simdat& sim) {
    int flag = 1;
    int iter = 0;

    double T2 = T;
    double T1 = Tlast;
    double T0 = sim.Tsol;
    double T_temp = 0.0;
    double T_err = 0.0;

    //double err_limit = sim.dttest;
    //int max_iter = sim.max_iter;
    double err_limit = 1.0/1000.0;
    double max_iter = 10;

    double time2 = t2;
    double time1 = t2 - sim.dt;
    double time0 = 0;
    double m = 0.0;

    vector<int_seg> isegv;

    while (flag) { //Continue until you get a temperature really close to the solidification temp
        m = (T2 - T1) / (time2 - time1);
        time0 = time1 + ((T0 - T1) / m);
        util::GetIntegrationInfo(time0, segv, sim, isegv,0);
        T_temp = Temp_Calc_Pre_Path(time0, isegv, sim, 0, 0);
        if (T_temp>T0) {
            time1 = time0;
            T1 = T_temp;
        }
        else {
            time2 = time0;
            T2 = T_temp;
        }
        T_err = fabs(1 - T_temp / T0);
        iter++;
        if ((T_err < err_limit)||(iter>=max_iter)) {
            flag = 0;
        }
        isegv.clear();
    }

    util::GetIntegrationInfo(time0, segv, sim, isegv, 1);
    T_temp = Temp_Calc_Pre_Path(time0, isegv, sim, 0, 1);

    if (Gtemp > (1e-50)) {
        G = Gtemp;

        double dTdt = sim.a*laplace +dT_cur;
        /*if (fabs(dT_cur / (sim.a*laplace)) > (1e-9)) { dTdt = (T2 - T1) / sim.dt; }
        else{ dTdt = sim.a*laplace + dT_cur; }*/

        dTdt_sol = dTdt;
        V = -dTdt_sol / G;

        Gx = Gx_temp;
        Gy = Gy_temp;
        Gz = Gz_temp;

        Gxu = Gx / G;
        Gyu = Gy / G;
        Gzu = Gz / G;

    }

    //sim.start_seg = seg_temp;
    output_flag = 1;
}


double Point::Temp_Calc_Pre_Path(double t2, vector<int_seg>& isegv, Simdat& sim, int set_T, int is_sol){
    double tau;
    double phix, phiy, phiz, psi, beta, phi, dTp;
    double dT = 0, dx, dy, dz, ct;
    Gx_temp = 0; Gy_temp = 0; Gz_temp = 0;
    laplace = 0;
    dT_cur = 0;

    //*************************************************************************************************
    //*************************************************************************************************
    //Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
    //Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
    beta = pow(3 / PI, 1.5) * sim.q / (sim.rho * sim.cps);	//Beam power pre-factor
    for (size_t i = 0; i < isegv.size(); i++){
        tau = t2 - isegv[i].taui;
        dx = x - isegv[i].xb;
        dy = y - isegv[i].yb;
        dz = z - isegv[i].zb;

        ct = 12 * sim.a * tau;
        phix = (sim.ax * sim.ax + ct);
        phiy = (sim.ay * sim.ay + ct);
        phiz = (sim.az * sim.az + ct);
        psi = 1 / (sqrt(phix*phiy*phiz));
        phi = exp(-(3 * dx * dx / phix) - (3 * dy * dy / phiy) - (3 * dz * dz / phiz));

        dTp = isegv[i].dtau * isegv[i].qmod * beta * phi * psi;
        dT += dTp;

        if (is_sol) {
            Gx_temp += dTp * (-6 * dx / phix);	//x-gradient
            Gy_temp += dTp * (-6 * dy / phiy);	//y-gradient
            Gz_temp += dTp * (-6 * dz / phiz);	//z-gradient

            laplace += (6 * dTp*(6 * (dx / phix*dx / phix + dy / phiy*dy / phiy + dz / phiz*dz / phiz) - 1.0 / phix - 1.0 / phiy - 1.0 / phiz)); //laplacian
            if (i == 0) {
                dT_cur = (isegv[i].qmod * beta * phi * psi);
            }// / (sim.rho*sim.cps); } //Add current heat source
        }
    }


    //*************************************************************************************************
    //*************************************************************************************************

    Gtemp = sqrt(Gx_temp*Gx_temp + Gy_temp*Gy_temp + Gz_temp*Gz_temp);

    if (Gtemp < 1e-50 || !isfinite(Gtemp)){ Gtemp = 1e-100; }
    if (set_T) {
        T = dT + sim.Tinit;
        T_calc_flag = 1;
    }
    if (T > sim.Tsol){
        G = 1e-100;
        V = 1e-100;
        output_flag = 1;
    }
    return dT + sim.Tinit;
}

std::string Point::Print()
{
    std::string str =
            std::to_string(x) + "," +
            std::to_string(y) + "," +
            std::to_string(z) + "," +
            std::to_string(T) + "," +
            std::to_string(G) + "," +
            std::to_string(V) + "\n";

    return str;
}
