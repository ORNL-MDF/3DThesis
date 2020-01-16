//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) � Transformer Resilience and Advanced Components (TRAC) Program.

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
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <climits>

#include "DataStructs.h"
#include "Point.h"
#include "Init.h"
#include "Util.h"

#ifdef _WIN32 // conditional needed to prevent error when compiling for unix systems
#include <direct.h>
#else // conditional for unix systems; should work *most of the time
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

void	Init::Keywords_Lv2(std::vector<std::string>& mainWords, std::vector<std::vector<std::string>>& subWords, std::vector<std::vector<std::string>>& strings, std::string& input_file) {
	
	for (int i = 0; i < subWords.size(); i++) {
		for (int j = 0; j < subWords[i].size(); j++) {
			strings[i].push_back("");
		}
	}
	
	std::string line;
	std::ifstream readFile;
	int num_read = 0;
	readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		readFile.open(input_file.c_str(), std::ios::in);
		while (true) {
			std::string tempWord = "";
			std::string tempWord2 = "";
			int mainWord = -1;
			std::string subWord = "";
			readFile >> tempWord; getline(readFile, line);
			for (int i = 0; i < mainWords.size(); i++) {
				if (tempWord == mainWords[i]) {
					mainWord = i;
					break;
				}
			}
			if (mainWord >= 0) {
				readFile >> tempWord; getline(readFile, line);
				if (tempWord == "{") {
					while (true) {
						readFile >> tempWord;
						if (tempWord != "}") {
							readFile >> tempWord2; getline(readFile, line);
							for (int i = 0; i < subWords[mainWord].size(); i++) {
								if (tempWord == subWords[mainWord][i]) {
									strings[mainWord][i] = tempWord2;
									num_read++;
								}
							}
						}
						else { getline(readFile, line); break; }
					}
				}
			}
		}
	}
	catch (const std::ifstream::failure&) { if (num_read) { std::cout << "Finished Reading " << input_file << "\n"; readFile.close();} }
}
void	Init::Keywords_Lv2(std::vector<std::string>& mainWords, std::vector<std::vector<std::string>>& subWords, std::vector<std::vector<double>>& values, std::string& input_file) {

	for (int i = 0; i < subWords.size(); i++) {
		for (int j = 0; j < subWords[i].size(); j++) {
			values[i].push_back(DBL_MAX);
		}
	}

	std::string line;
	std::ifstream readFile;
	int num_read = 0;
	readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		readFile.open(input_file.c_str(), std::ios::in);
		while (true) {
			std::string tempWord = "";
			double tempValue = DBL_MAX;
			int mainWord = -1;
			std::string subWord = "";
			readFile >> tempWord; getline(readFile, line);
			for (int i = 0; i < mainWords.size(); i++) {
				if (tempWord == mainWords[i]) {
					mainWord = i;
					break;
				}
			}
			if (mainWord >= 0) {
				readFile >> tempWord; getline(readFile, line);
				if (tempWord == "{") {
					while (true) {
						readFile >> tempWord;
						if (tempWord != "}") {
							readFile >> tempValue; getline(readFile, line);
							for (int i = 0; i < subWords[mainWord].size(); i++) {
								if (tempWord == subWords[mainWord][i]) {
									values[mainWord][i] = tempValue;
									num_read++;
								}
							}
						}
						else { getline(readFile, line); break; }
					}
				}
			}
		}
	}
	catch (const std::ifstream::failure&) { if (num_read){ std::cout << "Finished Reading "<< input_file << "\n"; readFile.close();}}
	
}

void	Init::SetValues(std::string& simValue, std::string input, std::string simDefault, std::string name, int err) {
	if (input == "") {
		simValue = simDefault; 
		if (err==1) {
			std::cout << "Error: " << name << " not read. Setting to default value " << simDefault << "\n";
		} 
		if (err == 2) {
			std::cout << "Fatal Error: " << name << " not read\n"; 
			//system("pause");
		}
	}
	else { simValue = input; }
	return;
}
void	Init::SetValues(int& simValue, double input, int simDefault, std::string name, int err) {
	if (input == DBL_MAX) {
		simValue = simDefault;
		if (err == 1) {
			std::cout << "Error: " << name << " not read. Setting to default value " << simDefault << "\n";
		}
		if (err == 2) {
			std::cout << "Fatal Error: " << name << " not read\n";
			//system("pause");
		}
	}
	else { simValue = int(input); }
	return;
}
void	Init::SetValues(double& simValue, double input, double simDefault, std::string name, int err) {
	if (input == DBL_MAX) {
		simValue = simDefault;
		if (err == 1) {
			std::cout << "Error: " << name << " not read. Setting to default value " << simDefault << "\n";
		}
		if (err == 2) {
			std::cout << "Fatal Error: " << name << " not read\n";
			//system("pause");
		}
	}
	else { simValue = input; }
	return;
}

void	Init::GetFileNames(Simdat& sim, std::string input_file) {
	std::vector<std::string> mainWords;
	mainWords.push_back("Simulation");
	mainWords.push_back("Options");
	mainWords.push_back("Utility");
	
	std::vector<std::vector<std::string>> subWords(mainWords.size());
	subWords[0].push_back("Name");
	subWords[0].push_back("Material");
	subWords[0].push_back("Beam");
	subWords[0].push_back("Path");

	subWords[1].push_back("Domain");
	subWords[1].push_back("Settings");

	subWords[2].push_back("ParBeams");

	std::vector<std::vector<std::string>> files(mainWords.size());

	Init::Keywords_Lv2(mainWords, subWords, files, input_file);

	Init::SetValues(sim.file.name, files[0][0], "N/A", "Simulation Name", 2);
	Init::SetValues(sim.file.mat, files[0][1], "", "Material File", 2);
	Init::SetValues(sim.file.beam, files[0][2], "", "Beam File", 2);
	Init::SetValues(sim.file.path, files[0][3], "", "Path File", 2);

	Init::SetValues(sim.file.domain, files[1][0], "", "Domain File", 0);
	Init::SetValues(sim.file.settings, files[1][1], "", "Settings File", 0);

	Init::SetValues(sim.file.parBeams, files[2][0], "", "ParBeams File", 0);

	return;
}
void	Init::ReadSimParams(std::vector<path_seg>& segv, Simdat& sim) {

	Init::MakeDataDirectory(sim);

	Init::FileRead_Path(segv, sim);
	Init::FileRead_Beam(sim);
	Init::FileRead_Material(sim);
	
	Init::FileRead_Domain(segv, sim);
	Init::FileRead_Settings(sim);

	Init::FileRead_ParBeams(sim);

	return;
}
void	Init::MakeDataDirectory(Simdat& sim) {
	//Make Data directory if it does not already exist
	std::string strPath = "Data";
#if defined(_WIN32)
	_mkdir(strPath.c_str());
#else 
	mkdir(strPath.c_str(), 0777); // notice that 777 is different than 0777
#endif

//Make SubData directory if it does not already exist
	std::string strSubPath = "Data/" + sim.file.name;
#if defined(_WIN32)
	_mkdir(strSubPath.c_str());
#else 
	mkdir(strSubPath.c_str(), 0777); // notice that 777 is different than 0777
#endif
}

void	Init::FileRead_Material(Simdat& sim) {
	std::vector<std::string> mainWords;
	mainWords.push_back("Constants");
	mainWords.push_back("CET");

	std::vector<std::vector<std::string>> subWords(mainWords.size());
	subWords[0].push_back("T_0");
	subWords[0].push_back("T_L");
	subWords[0].push_back("k");
	subWords[0].push_back("c");
	subWords[0].push_back("p");

	subWords[1].push_back("N0");
	subWords[1].push_back("n");
	subWords[1].push_back("a");

	std::vector<std::vector<double>> values(mainWords.size());

	Init::Keywords_Lv2(mainWords, subWords, values, sim.file.mat);

	Init::SetValues(sim.mat.Tinit, values[0][0], 1273.0, "Initial Temperature", 1);
	Init::SetValues(sim.mat.T_liq, values[0][1], 1610.0, "Liquidus Temperature", 1);
	Init::SetValues(sim.mat.kon, values[0][2], 26.6, "Thermal Conductivity", 1);
	Init::SetValues(sim.mat.cps, values[0][3], 600.00, "Specific Heat", 1);
	Init::SetValues(sim.mat.rho, values[0][4], 7451.0, "Density", 1);

	Init::SetValues(sim.mat.cet_N0, values[1][0], DBL_MAX, "CET: N0", 0);
	Init::SetValues(sim.mat.cet_n, values[1][1], DBL_MAX, "CET: n", 0);
	Init::SetValues(sim.mat.cet_a, values[1][2], DBL_MAX, "CET: a", 0);

	sim.mat.a = sim.mat.kon / (sim.mat.rho * sim.mat.cps);	//Calculate the thermal diffusivity
	sim.util.nond_dt = sim.beam.ax * sim.beam.ax / sim.mat.a; // Auto-Dtau Stuff

	sim.mat.calc_CET = 1;
	if (sim.mat.cet_N0 == DBL_MAX || sim.mat.cet_n == DBL_MAX || sim.mat.cet_a == DBL_MAX) { sim.mat.calc_CET = 0; }

	return;
}
void	Init::FileRead_Beam(Simdat& sim) {
	std::vector<std::string> mainWords;
	mainWords.push_back("Shape");
	mainWords.push_back("Intensity");

	std::vector<std::vector<std::string>> subWords(mainWords.size());
	subWords[0].push_back("Width_X");
	subWords[0].push_back("Width_Y");
	subWords[0].push_back("Depth_Z");

	subWords[1].push_back("Power");
	subWords[1].push_back("Efficiency");

	std::vector<std::vector<double>> values(mainWords.size());

	Init::Keywords_Lv2(mainWords, subWords, values, sim.file.beam);

	Init::SetValues(sim.beam.ax, values[0][0], 10.0e-6, "X Width", 1);
	Init::SetValues(sim.beam.ay, values[0][1], 10.0e-6, "Y Width", 1);
	Init::SetValues(sim.beam.az, values[0][2], 1.0e-6, "Z Depth", 1);

	Init::SetValues(sim.beam.q, values[1][0], 1200, "Power", 1);
	Init::SetValues(sim.beam.eff, values[1][1], 1.0, "Efficiency", 1);

	if (sim.beam.eff == DBL_MAX) { sim.beam.eff = 1.0; }
	sim.beam.q = sim.beam.q * sim.beam.eff * 2.0;

	return;
}
void	Init::FileRead_Path(std::vector<path_seg>& segv, Simdat& sim) {
	sim.param.xmin = DBL_MAX;
	sim.param.xmax = DBL_MIN;
	sim.param.ymin = DBL_MAX;
	sim.param.ymax = DBL_MIN;
	sim.param.zmin = -0.001;
	sim.param.zmax = 0;
	
	//Currently hard coding a conversion from mm to m
	double convert = 1e-3;

	std::ifstream pathfile;
	pathfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	std::string line;

	try {
		pathfile.open(sim.file.path.c_str(), std::ios::in);
		path_seg temp;

		//Set initial position as 0, 0, 0
		temp.smode = 1;
		temp.sx = 0.0;
		temp.sy = 0.0;
		temp.sz = 0.0;
		temp.sqmod = 0.0;
		temp.sparam = 0.0;
		segv.push_back(temp);

		//Read in path information from file
		while (getline(pathfile, line))
		{
			temp.smode = 1;
			temp.sx = 0.0;
			temp.sy = 0.0;
			temp.sz = 0.0;
			temp.sqmod = 0.0;
			temp.sparam = 0.0;

			pathfile >> temp.smode >> temp.sx >> temp.sy >> temp.sz >> temp.sqmod >> temp.sparam;
			temp.sx *= convert;
			temp.sy *= convert;
			temp.sz *= convert;
			segv.push_back(temp);
			if (temp.sqmod > 0) {
				if (temp.sx < sim.param.xmin) { sim.param.xmin = temp.sx; }
				if (temp.sx > sim.param.xmax) { sim.param.xmax = temp.sx; }
				if (temp.sy < sim.param.ymin) { sim.param.ymin = temp.sy; }
				if (temp.sy > sim.param.ymax) { sim.param.ymax = temp.sy; }
				if (temp.sz < sim.param.zmin) { sim.param.zmin = temp.sz; }
				if (temp.sz > sim.param.zmax) { sim.param.zmax = temp.sz; }
			}
		}
	}
	catch (const std::ifstream::failure&) {if (segv.size()) { std::cout << "Finished Reading " << sim.file.path << "\n"; pathfile.close(); }}
	sim.util.np = segv.size();

	//Calculate path times
	double dt_seg, dist, dx, dy, dz;
	segv[0].seg_time = 0;
	for (int p = 1; p < segv.size(); p++) {
		if (segv[p].smode) {	//For spot mode
			segv[p].seg_time = segv[p - 1].seg_time + segv[p].sparam;
		}
		else {						//For line mode			
			dx = segv[p].sx - segv[p - 1].sx;
			dy = segv[p].sy - segv[p - 1].sy;
			dz = segv[p].sz - segv[p - 1].sz;
			dist = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
			dt_seg = dist / segv[p].sparam;
			segv[p].seg_time = segv[p - 1].seg_time + dt_seg;
		}
	}


	sim.util.scanEndTime = segv[segv.size() - 1].seg_time;

	return;
}

void	Init::FileRead_Domain(std::vector<path_seg>& segv, Simdat& sim) {
	std::vector<std::string> mainWords;
	mainWords.push_back("X");
	mainWords.push_back("Y");
	mainWords.push_back("Z");
	mainWords.push_back("BoundaryConditions");

	std::vector<std::vector<std::string>> subWords(mainWords.size());
	subWords[0].push_back("Min");
	subWords[0].push_back("Max");
	subWords[0].push_back("Num");
	subWords[0].push_back("Res");

	subWords[1].push_back("Min");
	subWords[1].push_back("Max");
	subWords[1].push_back("Num");
	subWords[1].push_back("Res");

	subWords[2].push_back("Min");
	subWords[2].push_back("Max");
	subWords[2].push_back("Num");
	subWords[2].push_back("Res");

	subWords[3].push_back("X_min");
	subWords[3].push_back("X_max");
	subWords[3].push_back("Y_min");
	subWords[3].push_back("Y_max");

	std::vector<std::vector<double>> values(mainWords.size());

	Init::Keywords_Lv2(mainWords, subWords, values, sim.file.domain);

	Init::SetValues(sim.param.xmin, values[0][0], sim.param.xmin - 1e-3, "X min", 0);
	Init::SetValues(sim.param.xmax, values[0][1], sim.param.xmax + 1e-3, "X max", 0);
	Init::SetValues(sim.param.xnum, values[0][2], INT_MAX, "X num", 0);
	Init::SetValues(sim.param.xres, values[0][3], 50e-6, "X res", 0);

	Init::SetValues(sim.param.ymin, values[1][0], sim.param.ymin - 1e-3, "Y min", 0);
	Init::SetValues(sim.param.ymax, values[1][1], sim.param.ymax + 1e-3, "Y max", 0);
	Init::SetValues(sim.param.ynum, values[1][2], INT_MAX, "Y num", 0);
	Init::SetValues(sim.param.yres, values[1][3], 50e-6, "Y res", 0);

	Init::SetValues(sim.param.zmin, values[2][0], sim.param.zmin, "Z min", 0);
	Init::SetValues(sim.param.zmax, values[2][1], sim.param.zmax, "Z max", 0);
	Init::SetValues(sim.param.znum, values[2][2], INT_MAX, "Z num", 0);
	Init::SetValues(sim.param.zres, values[2][3], 50e-6, "Z res", 0);

	Init::SetValues(sim.param.BC_xmin, values[3][0], DBL_MAX, "BC X min", 0);
	Init::SetValues(sim.param.BC_xmax, values[3][1], DBL_MAX, "BC X max", 0);
	Init::SetValues(sim.param.BC_ymin, values[3][2], DBL_MAX, "BC Y min", 0);
	Init::SetValues(sim.param.BC_ymax, values[3][3], DBL_MAX, "BC Y max", 0);

	if (sim.param.xnum == INT_MAX) {
		sim.param.xnum = 1 + int((sim.param.xmax - sim.param.xmin) / sim.param.xres);
		sim.param.xmax = sim.param.xmin + (sim.param.xnum - 1) * sim.param.xres;
	}
	else if (sim.param.xnum != 1) {
		sim.param.xres = (sim.param.xmax - sim.param.xmin) / (sim.param.xnum - 1);
	}
	else {
		sim.param.xres == DBL_MAX;
	}

	if (sim.param.ynum == INT_MAX) {
		sim.param.ynum = 1 + int((sim.param.ymax - sim.param.ymin) / sim.param.yres);
		sim.param.ymax = sim.param.ymin + (sim.param.ynum - 1) * sim.param.yres;
	}
	else if (sim.param.ynum != 1) {
		sim.param.yres = (sim.param.ymax - sim.param.ymin) / (sim.param.ynum - 1);
	}
	else {
		sim.param.yres == DBL_MAX;
	}

	if (sim.param.znum == INT_MAX) {
		sim.param.znum = 1 + int((sim.param.zmax - sim.param.zmin) / sim.param.zres);
		sim.param.zmax = sim.param.zmin + (sim.param.znum - 1) * sim.param.zres;
	}
	else if (sim.param.ynum != 1) {
		sim.param.zres = (sim.param.zmax - sim.param.zmin) / (sim.param.znum - 1);
	}
	else {
		sim.param.zres == DBL_MAX;
	}

	if (sim.param.BC_xmin != DBL_MAX || sim.param.BC_xmax != DBL_MAX || sim.param.BC_ymin != DBL_MAX || sim.param.BC_ymax != DBL_MAX) { sim.setting.use_BCs = 1; }
	else { sim.setting.use_BCs = 0; }

	sim.param.pnum = sim.param.xnum * sim.param.ynum * sim.param.znum;

	return;
}
void	Init::FileRead_Settings(Simdat& sim) {
	std::vector<std::string> mainWords;
	mainWords.push_back("Simulation");
	mainWords.push_back("Output");
	mainWords.push_back("Temperature");
	mainWords.push_back("Path");
	mainWords.push_back("Neighbors");

	std::vector<std::vector<std::string>> subWords(mainWords.size());
	subWords[0].push_back("TimeStep");
	subWords[0].push_back("Mode");
	subWords[0].push_back("MaxThreads");
	subWords[0].push_back("PINT");

	subWords[1].push_back("Mode");
	subWords[1].push_back("Interval");

	subWords[2].push_back("Sol_Tol");
	subWords[2].push_back("Sol_Iter");
	subWords[2].push_back("Cutoff_Peak");
	subWords[2].push_back("Cutoff_T0TL");

	subWords[3].push_back("Buffer");
	subWords[3].push_back("Compression");

	subWords[4].push_back("Neighborhood");

	std::vector<std::vector<double>> values(mainWords.size());

	Init::Keywords_Lv2(mainWords, subWords, values, sim.file.settings);

	Init::SetValues(sim.param.dt, values[0][0], sim.util.nond_dt,"Timestep",0);
	Init::SetValues(sim.param.mode, values[0][1], 3, "Simulation Mode", 0);
	Init::SetValues(sim.setting.thnum, values[0][2], omp_get_max_threads()/2, "Number of Threads", 0);
	Init::SetValues(sim.param.use_PINT, values[0][3], 0, "Parallel in Time Mode", 0);

	Init::SetValues(sim.setting.out_mode, values[1][0], 1, "Output Mode", 0);
	Init::SetValues(sim.param.out_freq, values[1][1], INT_MAX, "Output Interval", 0);

	Init::SetValues(sim.setting.dttest, values[2][0], 1e-3, "Solidification Tolerance", 0);
	Init::SetValues(sim.setting.max_iter, values[2][1], 10, "Solidification Iterations", 0);
	Init::SetValues(sim.setting.t_hist, values[2][2], 1e-9, "Peak Temperature Cutoff Ratio", 0);
	Init::SetValues(sim.setting.p_hist, values[2][3], 1e-2, "Solidfication to Initial Temperature Cutoff Ratio", 0);

	Init::SetValues(sim.setting.r_max, values[3][0], -1.0, "Domain Buffer", 0);
	Init::SetValues(sim.setting.compress, values[3][1], 0, "Path Compression", 0);

	Init::SetValues(sim.setting.neighborhood, values[4][0], 1, "Point Neighborhood", 0);
	
	Util::CalcRMax(sim);

	return;
}

void	Init::FileRead_ParBeams(Simdat& sim){

	std::ifstream readFile;
	readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	std::string line;
	int num_read = 0;
	try {
		readFile.open(sim.file.parBeams.c_str(), std::ios::in);
		//getline(pathfile, line);
		parBeam temp;

		//Set initial position as 0, 0, 0
		temp.Xr = 0.0;
		temp.Yr = 0.0;
		temp.Pmod = 1.0;
		sim.beams.push_back(temp);

		//Read in path information from file
		while (getline(readFile, line))
		{

			readFile >> temp.Xr >> temp.Yr >> temp.Pmod;

			if (num_read){ sim.beams.push_back(temp); }
			else { sim.beams[0] = temp; }
			num_read++;
		}
	}
	catch (const std::ifstream::failure&) { if (num_read) { std::cout << "Finished Reading " << sim.file.parBeams << "\n"; readFile.close(); } }

	if (num_read) { sim.setting.parBeams = 1; }

	return;
}

void	Init::SetPoints(std::vector<Point>& ptv, Simdat& sim) {
	ptv.resize(sim.param.pnum);

	#pragma omp parallel num_threads(sim.setting.thnum)
	{
		double xp, yp, zp;
		#pragma omp for schedule(static)
		for (int i = 0; i < sim.param.xnum; i++) {
			if (sim.param.xnum == 1) { xp = sim.param.xmax; }
			else { xp = sim.param.xmin + ((float)i * ((sim.param.xmax - sim.param.xmin) / ((float)sim.param.xnum - 1))); }
			for (int j = 0; j < sim.param.ynum; j++) {
				if (sim.param.ynum == 1) { yp = sim.param.ymax; }
				else { yp = sim.param.ymin + ((float)j * ((sim.param.ymax - sim.param.ymin) / ((float)sim.param.ynum - 1))); }
				for (int k = 0; k < sim.param.znum; k++) {
					if (sim.param.znum == 1) { zp = sim.param.zmax; }
					else { zp = sim.param.zmin + ((float)k * ((sim.param.zmax - sim.param.zmin) / ((float)sim.param.znum - 1))); }
					int p = Util::ijk_to_p(i, j, k, sim);
					ptv[p].set_i(i);
					ptv[p].set_j(j);
					ptv[p].set_k(k);
					ptv[p].set_xloc(xp);
					ptv[p].set_yloc(yp);
					ptv[p].set_zloc(zp);
					ptv[p].Initialize();
					p++;
				}
			}
		}
	}
}