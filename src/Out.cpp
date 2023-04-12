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

void Out::Progress(const Simdat& sim, const int itert) {
	static int prog_print_last = 0;
	int prog_now = int(10.0 * itert * sim.param.dt / sim.util.allScansEndTime);
	if (prog_now != prog_print_last) {
		prog_print_last = prog_now;
		if (prog_print_last <= 10) {
			std::cout << "Time step: " << itert << "\t\t";
			std::cout << "% of Path: " << 10 * prog_print_last << "%" << "\n";
		}
		else{
			std::cout << "Time step: " << itert << "\t\t";
			std::cout << "Cooling... \n";
		}
	}
	return;
}

void Out::Point_Progress(const Simdat& sim, const int p) {
	static int prog_print_last = 0;
	int prog_now = 10 * p / sim.domain.pnum;
	if (prog_now != prog_print_last) {
		prog_print_last = prog_now;
		std::cout << "% of Points: " << 10 * prog_print_last << "%" << "\n";
	}
	
	return;
}

void Out::Write_csv_temp(Grid& grid, const Simdat& sim, const string name) {
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	string out_file = "Data/" + sim.files.name + "/" + sim.files.name + "." + name + ".csv";
	try {
		datafile.open(out_file.c_str());

		if (sim.output.x) {datafile << "x,"; }
		if (sim.output.y) { datafile << "y,"; }
		if (sim.output.z) { datafile << "z,"; }
		if (sim.output.T) { datafile << "T,"; }

		if (sim.output.G) { datafile << "G,"; }
		if (sim.output.Gx) { datafile << "Gx,"; }
		if (sim.output.Gy) { datafile << "Gy,"; }
		if (sim.output.Gz) { datafile << "Gz,"; }
		if (sim.output.V) { datafile << "V,"; }
		if (sim.output.dTdt) { datafile << "dTdt,"; }
		if (sim.output.eqFrac) { datafile << "eqFrac,"; }

		if (sim.output.H) { datafile << "H,"; }
		if (sim.output.Hx) { datafile << "Hx,"; }
		if (sim.output.Hy) { datafile << "Hy,"; }
		if (sim.output.Hz) { datafile << "Hz,"; }

		datafile << "\n";

		for (int p = 0; p < sim.domain.pnum; p++) {			
			if (grid.get_output_flag(p)){ 
				if (sim.output.x) { datafile << grid.get_x(p) << ","; }
				if (sim.output.y) { datafile << grid.get_y(p) << ","; }
				if (sim.output.z) { datafile << grid.get_z(p) << ","; }
				if (sim.output.T) { datafile << grid.get_T(p) << ","; }

				if (sim.output.G) { datafile << grid.get_G(p) << ","; }
				if (sim.output.Gx) { datafile << grid.get_Gx(p) << ","; }
				if (sim.output.Gy) { datafile << grid.get_Gy(p) << ","; }
				if (sim.output.Gz) { datafile << grid.get_Gz(p) << ","; }
				if (sim.output.V) { datafile << grid.get_V(p) << ","; }
				if (sim.output.dTdt) { datafile << grid.get_dTdt(p) << ","; }
				if (sim.output.eqFrac) { datafile << grid.get_eqFrac(p) << ","; }

				if (sim.output.H) { datafile << grid.get_H(p) << ","; }
				if (sim.output.Hx) { datafile << grid.get_Hx(p) << ","; }
				if (sim.output.Hy) { datafile << grid.get_Hy(p) << ","; }
				if (sim.output.Hz) { datafile << grid.get_Hz(p) << ","; }

				datafile << "\n";
			}
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
	return;
}

void Out::Write_csv(Grid& grid, const Simdat& sim, const string name) {
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	string out_file = "Data/" + sim.files.name + "/" + sim.files.name + "." + name + ".csv";
	try {
		datafile.open(out_file.c_str());
		datafile << "x,y,z,T,G,V\n";
		for (int p = 0; p < sim.domain.pnum; p++) {
			if (grid.get_output_flag(p)) {
				datafile << grid.get_x(p) << "," << grid.get_y(p) << "," << grid.get_z(p) << ",";
				datafile << grid.get_T(p) << "," << grid.get_G(p) << "," << grid.get_V(p) << "\n";
			}
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
	return;
}
//
//void Out::Write_T_hist(Grid& grid, Simdat& sim, string name) {
//	std::ofstream datafile;
//	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
//	string out_file = "Data/" + sim.files.name + "/" + sim.files.name + "." + name + ".csv";
//	int max_len = 0;
//	for (int p = 0; p < sim.domain.pnum; p++) {
//		int len = grid[p].get_T_hist().size();
//		if (len > max_len) {max_len = len;}
//	}
//	try {
//		datafile.open(out_file.c_str());
//		datafile << "x,y,z,";
//		for (int t_i = 0; t_i < max_len; t_i++) {datafile << "T(" << t_i * sim.param.dt << " s),";}
//		datafile << "\n";
//		for (int p = 0; p < sim.domain.pnum; p++) {
//			if (grid[p].get_output_flag()) {
//				datafile << grid[p].get_x() << "," << grid[p].get_y() << "," << grid[p].get_z() << ",";
//				vector<double> T_hist = grid[p].get_T_hist();
//				vector<int> T_hist_iter = grid[p].get_T_hist_iter();
//				int iter_temp = 0;
//				for (int iter = 0; iter < max_len; iter++) { 
//					if (T_hist_iter.size()==iter_temp) { datafile << T_hist[iter_temp - 1] << ","; }
//					else if (T_hist_iter[iter_temp] == iter) { datafile << T_hist[iter_temp] << ","; iter_temp++; }
//					else if (iter_temp) { datafile << T_hist[iter_temp - 1] << ","; }
//					else {datafile << sim.material.T_init << ","; }
//				}
//				datafile << "\n";
//			}
//		}
//	}
//	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
//	datafile.close();
//	return;
//}
//
//void Out::Write_Reduced(Grid& grid, Simdat& sim, string name) {
//	std::ofstream datafile;
//	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
//	string out_file = "Data/" + sim.files.name + "/" + sim.files.name + "." + name + ".csv";
//	try {
//		datafile.open(out_file.c_str());
//		datafile << "x,y,z,tm,tl,cr\n";
//		for (int p = 0; p < sim.domain.pnum; p++) {
//			if (grid.get_output_flag(p)){
//				//const int numEvents = grid.get_
//				/*for (int e = 0; e < grid[p].get_numEvents()[1]; e++) {
//					datafile << grid.get_x(p) << "," << grid.get_y(p) << "," << grid.get_z(p) << ",";
//					datafile << grid[p].get_t_melt(e) << "," << grid[p].get_t_liq(e) << "," << grid[p].get_cr(e) << "\n";
//				}*/
//			}
//		}
//	}
//	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
//	datafile.close();
//	return;
//}