//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

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

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <cmath>

#include "Out.h"
#include "DataStructs.h"
#include "Point.h"

void Out::Progress(Simdat& sim, int itert) {
	if (int(10.0 * itert * sim.param.dt / sim.util.scanEndTime) != sim.util.prog_print_last) {
		sim.util.prog_print_last = int(10.0 * itert * sim.param.dt / sim.util.scanEndTime);
		if (sim.util.prog_print_last <= 10) {
			std::cout << "Time step: " << itert << "\t\t";
			std::cout << "% of Path: " << 10 * sim.util.prog_print_last << "%" << "\n";
		}
		else{
			std::cout << "Time step: " << itert << "\t\t";
			std::cout << "Cooling... \n";
		}
	}
	return;
}

void Out::Point_Progress(Simdat& sim, int p) {
	if (10 * p / sim.param.pnum != sim.util.prog_print_last) {
		sim.util.prog_print_last = 10 * p / sim.param.pnum;
		std::cout << "% of Points: " << 10 * sim.util.prog_print_last << "%" << "\n";
	}
	
	return;
}

void Out::Write_csv(Point * const ptv, Simdat& sim, string name, int out_mode) {
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	string out_file = "Data/" + sim.file.name + "/" + sim.file.name + "." + name + ".csv";
	try {
		datafile.open(out_file.c_str());
		if (out_mode == 0) {
			datafile << "x,y,z,T\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "\n";
				}
			}
		}
		else if (out_mode == 1) {
			datafile << "x,y,z,T,G,V\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << "\n";
				}
			}
		}
		else if (out_mode == 2) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << ",";
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << "\n";
				}
			}
		}
		else if (out_mode == 3) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac,Gx,Gy,Gz\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << "," ;
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << ",";
					datafile << ptv[p].get_Gxu() << "," << ptv[p].get_Gyu() << "," << ptv[p].get_Gzu() << "\n";
				}
			}
		}
		else if (out_mode == 4) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac,Gx,Gy,Gz,H,Hx,Hy,Hz\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << ",";
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << ",";
					datafile << ptv[p].get_Gxu() << "," << ptv[p].get_Gyu() << "," << ptv[p].get_Gzu() << ",";
					datafile << ptv[p].get_H() << "," << ptv[p].get_Hxu() << "," << ptv[p].get_Hyu() << "," << ptv[p].get_Hzu() << "\n";
				}
			}
		}
		else {
			datafile << "x,y,z,tl,cr\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag() && ptv[p].get_t_last_liq() > 0.0) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << "," << ptv[p].get_t_last_liq() << "," << ptv[p].get_dTdt_sol() << "\n";
				}
			}
			/*datafile << "x,y,z,tl,ts\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag() && ptv[p].get_t_last_liq() > 0.0) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << "," << ptv[p].get_t_last_liq() << "," << ptv[p].get_t_last_sol() << "\n";
				}
			}*/
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
	return;
}

void Out::Write_csv_PINT(std::deque<Point>& ptv, Simdat& sim, string name, int out_mode) {
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	string out_file = "Data/" + sim.file.name + "/" + sim.file.name + "." + name + ".csv";
	try {
		datafile.open(out_file.c_str());
		if (out_mode == 0) {
			datafile << "x,y,z,T\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "\n";
				}
			}
		}
		else if (out_mode == 1) {
			datafile << "x,y,z,T,G,V\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << "\n";
				}
			}
		}
		else if (out_mode == 2) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << ",";
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << "\n";
				}
			}
		}
		else if (out_mode == 3) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac,Gx,Gy,Gz\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << ",";
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << ",";
					datafile << ptv[p].get_Gxu() << "," << ptv[p].get_Gyu() << "," << ptv[p].get_Gzu() << "\n";
				}
			}
		}
		else if (out_mode == 4) {
			datafile << "x,y,z,T,G,V,dTdt,eq_frac,Gx,Gy,Gz,H,Hx,Hy,Hz\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag()) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
					datafile << ptv[p].get_T() << "," << ptv[p].get_G() << "," << ptv[p].get_V() << ",";
					datafile << -ptv[p].get_dTdt_sol() << "," << ptv[p].get_eq_frac() << ",";
					datafile << ptv[p].get_Gxu() << "," << ptv[p].get_Gyu() << "," << ptv[p].get_Gzu() << ",";
					datafile << ptv[p].get_H() << "," << ptv[p].get_Hxu() << "," << ptv[p].get_Hyu() << "," << ptv[p].get_Hzu() << "\n";
				}
			}
		}
		else {
			datafile << "x,y,z,tl,cr\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag() && ptv[p].get_t_last_liq() > 0.0) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << "," << ptv[p].get_t_last_liq() << "," << ptv[p].get_dTdt_sol() << "\n";
				}
			}
			/*datafile << "x,y,z,tl,ts\n";
			for (int p = 0; p < sim.param.pnum; p++) {
				if (ptv[p].get_output_flag() && ptv[p].get_t_last_liq() > 0.0) {
					datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << "," << ptv[p].get_t_last_liq() << "," << ptv[p].get_t_last_sol() << "\n";
				}
			}*/
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
	return;
}

void Out::Write_T_hist(Point * const ptv, Simdat& sim, string name) {
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	string out_file = "Data/" + sim.file.name + "/" + sim.file.name + "." + name + ".csv";
	int max_len = 0;
	for (int p = 0; p < sim.param.pnum; p++) {
		int len = ptv[p].get_T_hist().size();
		if (len > max_len) {max_len = len;}
	}
	try {
		datafile.open(out_file.c_str());
		datafile << "x,y,z,";
		for (int t_i = 0; t_i < max_len; t_i++) {datafile << "T(" << t_i * sim.param.dt << " s),";}
		datafile << "\n";
		for (int p = 0; p < sim.param.pnum; p++) {
			if (ptv[p].get_output_flag()) {
				datafile << ptv[p].get_x() << "," << ptv[p].get_y() << "," << ptv[p].get_z() << ",";
				vector<double> T_hist = ptv[p].get_T_hist();
				vector<int> T_hist_iter = ptv[p].get_T_hist_iter();
				int iter_temp = 0;
				for (int iter = 0; iter < max_len; iter++) { 
					if (T_hist_iter.size()==iter_temp) { datafile << T_hist[iter_temp - 1] << ","; }
					else if (T_hist_iter[iter_temp] == iter) { datafile << T_hist[iter_temp] << ","; iter_temp++; }
					else if (iter_temp) { datafile << T_hist[iter_temp - 1] << ","; }
					else {datafile << sim.mat.Tinit << ","; }
				}
				datafile << "\n";
			}
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
	return;
}