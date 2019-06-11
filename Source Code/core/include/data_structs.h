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

#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <string>
#include <fstream>
#include <iostream>

struct int_seg{
    double xb, yb, zb, taui, dtau, qmod;
};

class PathSeg {
public:
    PathSeg(){}
    PathSeg(double x, double y, double param) : PathSeg(1, x, y, 0.0, 1.0, param){}
    PathSeg(int mode, double x, double y, double z, double qmod, double param)
        : smode(mode), sx(x), sy(y), sz(z), sqmod(qmod), sparam(param), seg_time(0.0) {}

    int smode;
    double sx, sy, sz;		//Segment start coordinates
    double sqmod;		//Segmod power modulation
    double sparam;
    double seg_time;
};

class Simdat {
public:

    Simdat(){}

    Simdat(std::string path)
    {
        readSimFile(path);
    }

    void readSimFile(std::string path)
    {
        std::string line;
        std::ifstream simfile;
        simfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        int count = 0;
        try{
            simfile.open(path, std::ios::in);
            simfile >> sim_name; std::getline(simfile, line); count++;
            simfile >> out_path; std::getline(simfile, line); count++;
            simfile >> thnum; getline(simfile, line); count++;	 		//Number of threads to use in parallel
            simfile >> imax >> jmax >> kmax; getline(simfile, line); count++;	//Number of points in each direction
            simfile >> xmin >> xmax; getline(simfile, line); count++;	//x-direction limits
            simfile >> ymin >> ymax; getline(simfile, line); count++;   //y-direction limits
            simfile >> zmin >> zmax; getline(simfile, line); count++;   //z-direction limits
            simfile >> mode; getline(simfile, line); count++;   //Simulation mode (1 - normal, 2 - GV only)
            simfile >> dt; getline(simfile, line); count++;	//Time step
            simfile >> out_freq; getline(simfile, line); count++;	//Output frequency in number of time steps
            //material part
            simfile >> Tinit >> Tsol; getline(simfile, line); count++; //Pre-heat and Solidification Temperatures
            simfile >> cps; getline(simfile, line);	count++;//Specific heat
            simfile >> kon; getline(simfile, line);	count++;//Thermal conductivity
            simfile >> rho; getline(simfile, line);	count++;//Density
            //beam size
            simfile >> ax >> ay >> az; getline(simfile, line); count++;//Beam widths
            simfile >> eff; getline(simfile, line); count++;		//Efficiency of absorption
            simfile >> q; getline(simfile, line); count++;		//Power
            q = q * 2 * eff;					//Apply efficiency and boundary conditions factors
            simfile.close();
        }
        catch (const std::ifstream::failure& e){
			std::cout << e.what() << std::endl;
			std::cout << "Exception opening/reading sim file at " << count << std::endl;

        }
        tmod = 1.0001; //Sets tmod to 1.0001
        t_hist = 1e+9; //Sets thist to 1e+9
        pnum = imax*jmax*kmax;	//Claculate total number of points
    }


    std::string sim_name; //Name for output data files
    std::string out_path; //Directory of output data files
    //Material properties
    double kon, rho, cps, cpl, Lf, Tliq, Tsol, Tinit;
    double a;

    //Parameter for columnar to equiaxed transition
    double cet_a, cet_n, cet_N0;

    //Path parameters
    int np;

    //Beam parameters
    double ax, ay, az;
    double eff, q;

    //Simulation parameters
    int thnum;
    int imax, jmax, kmax, pnum;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int mode, out_freq;
    double dt, tmod, t_hist;
    int ptmode;

    //Pre-flight parameters
    int	PF_file_flag;
    int Rnum, tnum;
    double Rmax, R_time, PFdt, spot_time;

    //Integration mode
    int int_mode;

    //Other params
    int start_seg;
    double nond_dt;
    double max_speed;

    //Resolutions
    double xres, yres, zres;
};

#endif // DATASTRUCTS_H
