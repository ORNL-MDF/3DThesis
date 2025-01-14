/****************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC                                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of 3dThesis. 3dThesis is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <cfloat>
#include <climits>

#include "impl/Structs/DataStructs.hpp"
#include "impl/IO/Init.hpp"
#include "impl/Calc/Util.hpp"

#ifdef _WIN32 // conditional needed to prevent error when compiling for unix systems
#include <direct.h>
#else // conditional for unix systems; should work *most of the time
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

namespace Thesis::impl
{
	void	Init::Keywords_Lv1(vector<string>& mainWords, vector<int>& isIn, const string& input_file) {
		string line;
		std::ifstream readFile;
		readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		try {
			readFile.open(input_file.c_str(), std::ios::in);
			while (true) {
				string tempWord = "";	
				string tempWord2 = "";
				int mainWord = -1;
				string subWord = "";
				readFile >> tempWord; getline(readFile, line);
				for (int i = 0; i < mainWords.size(); i++) {
					if (tempWord == mainWords[i]) {
						mainWord = i;
						isIn[i] += 1;
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
							}
							else { getline(readFile, line); break; }
						}
					}
				}
			}
		}
		catch (const std::ifstream::failure&) {}
		catch (const std::exception& e) {}
	}

	void Init::catchPrint(const string& input_file, const int num_read, const bool print){
		if (print) {
			if (num_read) { std::cout << "Finished Reading " << input_file << std::endl;}
			else { std::cout << "Failed to Open " << input_file << std::endl;}
		}
	}

	void	Init::Keywords_Lv2(vector<string>& mainWords, vector<vector<string>>& subWords, vector<vector<string>>& strings, const string& input_file, const bool print) {
		
		for (int i = 0; i < subWords.size(); i++) {
			for (int j = 0; j < subWords[i].size(); j++) {
				strings[i].push_back("");
			}
		}
		
		string line;
		std::ifstream readFile;
		int num_read = 0;
		readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		try {
			readFile.open(input_file.c_str(), std::ios::in);
			while (true) {
				string tempWord = "";
				string tempWord2 = "";
				int mainWord = -1;
				string subWord = "";
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
		catch (const std::ifstream::failure&) {
			catchPrint(input_file, num_read, print);
			readFile.close();
		}
		catch (const std::exception& e) {
			catchPrint(input_file, num_read, print);
			readFile.close();
		}
	}
	void	Init::Keywords_Lv2(vector<string>& mainWords, vector<vector<string>>& subWords, vector<vector<double>>& values, const string& input_file, const bool print) {

		for (int i = 0; i < subWords.size(); i++) {
			for (int j = 0; j < subWords[i].size(); j++) {
				values[i].push_back(DBL_MAX);
			}
		}

		string line;
		std::ifstream readFile;
		int num_read = 0;
		readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		try {
			readFile.open(input_file.c_str(), std::ios::in);
			while (true) {
				string tempWord = "";
				double tempValue = DBL_MAX;
				int mainWord = -1;
				string subWord = "";
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
		catch (const std::ifstream::failure&) { 
			catchPrint(input_file, num_read, print);
			readFile.close();
		}	
		catch (const std::exception& e) {
			catchPrint(input_file, num_read, print);
			readFile.close();
		}
	}

	void	Init::SetValues(string& simValue, string input, string simDefault, string name, int err, const bool print) {
		if (input == string("")) {
			simValue = simDefault; 
			setPrint(simDefault, name, err, print);
		}
		else { simValue = input; }
		return;
	}
	void	Init::SetValues(int& simValue, string input, int simDefault, string name, int err, const bool print) {
		if (input == string("")) {
			simValue = simDefault;
			setPrint(simDefault, name, err, print);
		}
		else { 
			simValue = std::stoi(input);
		}
		return;
	}
	void	Init::SetValues(bool& simValue, string input, bool simDefault, string name, int err, const bool print) {
		if (input == string("")) {
			simValue = simDefault;
			setPrint(simDefault, name, err, print);
		}
		else { simValue = bool(std::stoi(input)); }
		return;
	}
	void	Init::SetValues(double& simValue, string input, double simDefault, string name, int err, const bool print) {
		if (input == string("")) {
			simValue = simDefault;
			Init::setPrint(simDefault, name, err, print);
		}
		else { simValue = std::stod(input); }
		return;
	}

	void    Init::checkAsterisks(const std::string s, const std::string del, const std::string file, const bool print) {
		// If there are multiple '*', there is a problem
		if (s.find(del) != string::npos) {
			if (print)
				std::cout << "Too many '*' in: " << file << "\n";
			exit(1);
		}
	}

	void	Init::GetFileNames(FileNames& files, const string& file, const bool print) {
		
		size_t pos = file.find_last_of("/");
		files.dataDir = "Data";
		if (pos != string::npos) { files.dataDir = file.substr(0, pos) + "/Data";}
		Init::MakeDataDirectory(files.dataDir);
		
		vector<string> mainWords;
		mainWords.push_back("Simulation");
		mainWords.push_back("Options");
		//mainWords.push_back("Utility");
		
		vector<vector<string>> subWords(mainWords.size());
		subWords[0].push_back("Name");
		subWords[0].push_back("Mode");
		subWords[0].push_back("Material");
		subWords[0].push_back("Beam");
		subWords[0].push_back("Path");

		subWords[1].push_back("Domain");
		subWords[1].push_back("Output");
		subWords[1].push_back("Settings");

		//subWords[2].push_back("Points");
		//subWords[2].push_back("ParBeams");
		//subWords[2].push_back("InfBeams");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(files.name, values[0][0], "N/A", "Simulation Name", 2, print);
		Init::SetValues(files.mode, values[0][1], "", "Mode File", 2, print);
		Init::SetValues(files.material, values[0][2], "", "Material File", 2, print);
		Init::SetValues(files.beam, values[0][3], "", "Beam File", 2, print);
		Init::SetValues(files.path, values[0][4], "", "Path File", 2, print);


		Init::SetValues(files.domain, values[1][0], "", "Domain File", 0, print);
		Init::SetValues(files.output, values[1][1], "", "Output File", 0, print);
		Init::SetValues(files.settings, values[1][2], "", "Settings File", 0, print);

		//Init::SetValues(sim.files.points, files[2][0], "", "Points File", 0, print);
		//Init::SetValues(sim.files.parBeams, files[2][1], "", "ParBeams File", 0, print);
		//Init::SetValues(sim.files.infBeams, files[2][2], "", "InfBeams File", 0, print);

		return;
	}
	void	Init::MakeDataDirectory(const string& file) {
		//Make Data directory if it does not already exist
		#if defined(_WIN32)
			_mkdir(file.c_str());
		#else 
			mkdir(file.c_str(), 0777); // notice that 777 is different than 0777
		#endif
	}

	void	Init::ReadSimParams(Simdat& sim) {

		Init::FileRead_Beams(sim.beams, sim.files.beam, sim.print);				// Read beam files
		Init::FileRead_Material(sim.material, sim.files.material, sim.print); Util::Calc_NonD_dt(sim.beams, sim.material);		// Read material properties // Calculate nonDimensional integration time for all beams
		Init::FileRead_Paths(sim.paths, sim.files.path, sim.print); Util::Calc_AllScansEndTime(sim);	// Read path files	// Calculate Important Simulation Parameters
		
		Init::FileRead_Mode(sim, sim.files.mode);	// Initialize simulation mode

		Init::FileRead_Domain(sim.domain, sim.files.domain, sim.print); Util::Calc_ScanBounds(sim.domain, sim.paths); Init::SetDomainParams(sim.domain);	// Read domain // Look at bounds of scan path and make domain bounds if bounds are unspecified // Set domain parameters
		Init::FileRead_Output(sim.output, sim.files.output, sim.print);	// Read what to output
		Init::FileRead_Settings(sim.settings, sim.files.settings, sim.print); Util::Calc_RMax(sim);	// Read fine tuned settings // Calculate maximum radius to care about

		return;
	}

	void	Init::FileRead_Mode(Simdat& sim, const string& file) {
		vector<string> mainWords;
		mainWords.push_back("Snapshots");
		mainWords.push_back("Solidification");

		vector<int> hasMainWords(mainWords.size());
		Keywords_Lv1(mainWords, hasMainWords, file);

		int sum = 0;
		int loc = -1;
		for (int i = 0; i < hasMainWords.size(); i++) {
			sum += hasMainWords[i];
			if (hasMainWords[i] > 0) {loc = i;}
		}

		if (sum > 1) {
			if (sim.print)
				std::cout << "Fatal Error: Too many modes in file " << file << std::endl;
			exit(1);
		}
		if (loc < 0) {
			if (sim.print)
				std::cout << "Fatal Error: No modes in file " << file << std::endl;
			exit(1);
		}

		sim.param.mode = "";
		switch (loc) {
			case 0:
				FileRead_Mode_Snapshot(sim, file);
				break;
			case 1:
				FileRead_Mode_Solidification(sim, file);
				break;
		}
	}
	void	Init::FileRead_Mode_Snapshot(Simdat& sim, const string& file) {
		sim.param.mode = "Snapshots";
		
		vector<string> mainWords;
		mainWords.push_back("Snapshots");

		vector<vector<string>> subWords(mainWords.size());
		subWords[0].push_back("Times");
		subWords[0].push_back("ScanFracs");
		subWords[0].push_back("Tracking");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, sim.print);

		if (values[0][0].size() && values[0][1].size()) {
			if (sim.print)
				std::cout << "Fatal Error: Too many time types for snapshot in " << file << std::endl;
			exit(1);
		}

		// For Times
		if (values[0][0].size()) {
			string s = values[0][0];
			string tmp;
			stringstream ss(s);
			while (getline(ss, tmp, ',')) {
				sim.param.SnapshotTimes.push_back(std::stod(tmp));
			}
		}

		// For ScanFracs
		if (values[0][1].size()) {
			string s = values[0][1];
			string tmp;
			stringstream ss(s);
			while (getline(ss, tmp, ',')) {
				sim.param.SnapshotTimes.push_back(std::stod(tmp)/100.0*sim.util.allScansEndTime);
			}
		}

		Init::SetValues(sim.param.tracking, values[0][2], string("None"), "Tracking", 0, sim.print);

	}
	void	Init::FileRead_Mode_Solidification(Simdat& sim, const string& file) {
		sim.param.mode = "Solidification";
		
		vector<string> mainWords;
		mainWords.push_back("Solidification");

		vector<vector<string>> subWords(mainWords.size());

		subWords[0].push_back("Tracking");			// None, Volume, Perimeter
		subWords[0].push_back("Timestep");
		subWords[0].push_back("OutputFrequency");
		subWords[0].push_back("CheckRadius");
		subWords[0].push_back("Secondary");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, sim.print);

		Init::SetValues(sim.param.tracking, values[0][0], string("None"), "Tracking", 0, sim.print);
		Init::SetValues(sim.param.dt, values[0][1], 1e-5, "Timestep", 0, sim.print);
		Init::SetValues(sim.param.out_freq, values[0][2], INT_MAX, "Output Frequency", 0, sim.print);
		Init::SetValues(sim.param.radiusCheck, values[0][3], -1.0, "Check Radius", 0, sim.print);
		Init::SetValues(sim.param.secondary, values[0][4], 0, "Secondary Solidfication", 0, sim.print);
	}
	
	void	Init::FileRead_Material(Material& material, const string& file, const bool print) {
		vector<string> mainWords;
		mainWords.push_back("Constants");
		mainWords.push_back("CET");

		vector<vector<string>> subWords(mainWords.size());
		subWords[0].push_back("T_0");
		subWords[0].push_back("T_L");
		subWords[0].push_back("k");
		subWords[0].push_back("c");
		subWords[0].push_back("p");

		subWords[1].push_back("N0");
		subWords[1].push_back("n");
		subWords[1].push_back("a");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(material.T_init, values[0][0], 1273.0, "Initial Temperature", 1, print);
		Init::SetValues(material.T_liq, values[0][1], 1610.0, "Liquidus Temperature", 1, print);
		Init::SetValues(material.kon, values[0][2], 26.6, "Thermal Conductivity", 1, print);
		Init::SetValues(material.cps, values[0][3], 600.00, "Specific Heat", 1, print);
		Init::SetValues(material.rho, values[0][4], 7451.0, "Density", 1, print);

		Init::SetValues(material.cet_N0, values[1][0], DBL_MAX, "CET: N0", 0, print);
		Init::SetValues(material.cet_n, values[1][1], DBL_MAX, "CET: n", 0, print);
		Init::SetValues(material.cet_a, values[1][2], DBL_MAX, "CET: a", 0, print);

		Init::SetDiffusivity(material);

		return;
	}

	void	Init::FileRead_Beams(vector<Beam>& beams, const string& file, const bool print) {

		// Find location of *
		string s = file;
		string del = "*";
		size_t pos = s.find(del);

		// If no * found, read file as is
		if (pos == string::npos) {
			Beam beam;
			Init::FileRead_Beam(beam, file, print);
			beams.push_back(beam);
		}
		// If there is a * found, start reading from ('*' = 1)
		else {
			// Find what comes before (s1) and after (s2) the '*'
			const string s1 = s.substr(0, pos);
			s.erase(0, pos + del.length());
			const string s2 = s.substr(0, string::npos);

			checkAsterisks(s, del, file, print);
			
			// Start reading from '*' = 1
			uint16_t beamNum = 1;
			while (true) {
				std::ifstream readFile;
				readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
				const string wildFile = s1 + to_string(beamNum) + s2;
				try {
					readFile.open(wildFile.c_str(), std::ios::in);
				}
				catch (const std::ifstream::failure&) {
					break;
				}
				Beam beam;
				Init::FileRead_Beam(beam, wildFile, print);
				beams.push_back(beam);
				beamNum++; 
			}
		}

		return;
	}
	void	Init::FileRead_Beam(Beam& beam, const string& file, const bool print) {
		vector<string> mainWords;
		mainWords.push_back("Shape");
		mainWords.push_back("Intensity");

		vector<vector<string>> subWords(mainWords.size());
		subWords[0].push_back("Width_X");
		subWords[0].push_back("Width_Y");
		subWords[0].push_back("Depth_Z");

		subWords[1].push_back("Power");
		subWords[1].push_back("Efficiency");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(beam.ax, values[0][0], 10.0e-6, "X Width", 1, print);
		Init::SetValues(beam.ay, values[0][1], 10.0e-6, "Y Width", 1, print);
		Init::SetValues(beam.az, values[0][2], 1.0e-6, "Z Depth", 1, print);

		Init::SetValues(beam.q, values[1][0], 1200, "Power", 1, print);
		Init::SetValues(beam.eff, values[1][1], 1.0, "Efficiency", 1, print);

		Init::SetBeamPower(beam);

		return;
	}

	void	Init::FileRead_Paths(vector<vector<path_seg>>& paths, const string& file, const bool print) {

		// Find location of *
		string s = file;
		string del = "*";
		size_t pos = s.find(del);

		// If no * found, read file as is
		if (pos == string::npos) {
			vector<path_seg> path;
			Init::FileRead_Path(path, file, print);
			paths.push_back(path);
		}
		// If there is a * found, start reading from ('*' = 1)
		else {
			// Find what comes before (s1) and after (s2) the '*'
			const string s1 = s.substr(0, pos);
			s.erase(0, pos + del.length());
			const string s2 = s.substr(0, string::npos);

			checkAsterisks(s, del, file, print);

			// Start reading from '*' = 1
			uint16_t pathNum = 1;
			while (true) {
				std::ifstream readFile;
				readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
				const string wildFile = s1 + to_string(pathNum) + s2;
				try {
					readFile.open(wildFile.c_str(), std::ios::in);
				}
				catch (const std::ifstream::failure&) {
					break;
				}
				vector<path_seg> path;
				Init::FileRead_Path(path, wildFile, print);
				paths.push_back(path);
				pathNum++;
			}
		}

		return;
	}
	void	Init::FileRead_Path(vector<path_seg>& path, const string& file, const bool print) {	
		//Currently hard coding a conversion from mm to m
		double convert = 1e-3;

		std::ifstream pathfile;
		pathfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		string line;

		try {
			pathfile.open(file.c_str(), std::ios::in);
			path_seg seg;

			//Set initial position as 0, 0, 0
			seg.smode = 1;
			seg.sx = 0.0;
			seg.sy = 0.0;
			seg.sz = 0.0;
			seg.sqmod = 0.0;
			seg.sparam = 0.0;
			path.push_back(seg);

			//Read in path information from file
			while (getline(pathfile, line))
			{
				seg.smode = 1;
				seg.sx = 0.0;
				seg.sy = 0.0;
				seg.sz = 0.0;
				seg.sqmod = 0.0;
				seg.sparam = 0.0;

				pathfile >> seg.smode >> seg.sx >> seg.sy >> seg.sz >> seg.sqmod >> seg.sparam;
				seg.sx *= convert;
				seg.sy *= convert;
				seg.sz *= convert;
				path.push_back(seg);
			}
		}
		catch (const std::ifstream::failure&) {
			catchPrint(file, path.size(), print);
			pathfile.close();
		}
		catch (const std::exception& e) {
			catchPrint(file, path.size(), print);
			pathfile.close();
		}

		//Calculate path times
		double dt_seg, dist, dx, dy, dz;
		path[0].seg_time = 0;
		for (int seg = 1; seg < path.size(); seg++) {
			if (path[seg].smode) {	//For spot mode
				path[seg].seg_time = path[seg - 1].seg_time + path[seg].sparam;
			}
			else {						//For line mode			
				dx = path[seg].sx - path[seg - 1].sx;
				dy = path[seg].sy - path[seg - 1].sy;
				dz = path[seg].sz - path[seg - 1].sz;
				dist = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
				dt_seg = dist / path[seg].sparam;
				path[seg].seg_time = path[seg - 1].seg_time + dt_seg;
			}
		}

		return;
	}

	void	Init::FileRead_Domain(Domain& domain, const string& file, const bool print) {
		vector<string> mainWords;
		mainWords.push_back("X");
		mainWords.push_back("Y");
		mainWords.push_back("Z");
		mainWords.push_back("BoundaryConditions");
		mainWords.push_back("Custom");

		vector<vector<string>> subWords(mainWords.size());
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
		subWords[3].push_back("Z_min");
		subWords[3].push_back("Reflections");

		subWords[4].push_back("File");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(domain.xmin, values[0][0], -DBL_MAX, "X min", 0, print);
		Init::SetValues(domain.xmax, values[0][1], DBL_MAX, "X max", 0, print);
		Init::SetValues(domain.xnum, values[0][2], INT_MAX, "X num", 0, print);
		Init::SetValues(domain.xres, values[0][3], 50e-6, "X res", 0, print);

		Init::SetValues(domain.ymin, values[1][0], -DBL_MAX, "Y min", 0, print);
		Init::SetValues(domain.ymax, values[1][1], DBL_MAX, "Y max", 0, print);
		Init::SetValues(domain.ynum, values[1][2], INT_MAX, "Y num", 0, print);
		Init::SetValues(domain.yres, values[1][3], 50e-6, "Y res", 0, print);

		Init::SetValues(domain.zmin, values[2][0], -DBL_MAX, "Z min", 0, print);
		Init::SetValues(domain.zmax, values[2][1], DBL_MAX, "Z max", 0, print);
		Init::SetValues(domain.znum, values[2][2], INT_MAX, "Z num", 0, print);
		Init::SetValues(domain.zres, values[2][3], 50e-6, "Z res", 0, print);

		Init::SetValues(domain.BC_xmin, values[3][0], DBL_MAX, "BC X min", 0, print);
		Init::SetValues(domain.BC_xmax, values[3][1], DBL_MAX, "BC X max", 0, print);
		Init::SetValues(domain.BC_ymin, values[3][2], DBL_MAX, "BC Y min", 0, print);
		Init::SetValues(domain.BC_ymax, values[3][3], DBL_MAX, "BC Y max", 0, print);
		Init::SetValues(domain.BC_zmin, values[3][4], DBL_MAX, "BC Z min", 0, print);
		Init::SetValues(domain.BC_reflections, values[3][5], INT_MAX, "BC reflections", 0, print);

		Init::SetValues(domain.pointsFile, values[4][0], string(""), "Point File", 0, print);

		// For boundary conditions
		if (domain.BC_xmin != DBL_MAX || domain.BC_xmax != DBL_MAX || domain.BC_ymin != DBL_MAX || domain.BC_ymax != DBL_MAX || domain.BC_zmin != DBL_MAX) {
			domain.use_BCs = 1;
		}
		else {
			domain.use_BCs = 0;
		}

		// For custom point files
		if (domain.pointsFile == string("")) {
			domain.customPoints = false;
		}
		else {
			domain.customPoints = true;
			FileRead_Points(domain, domain.pointsFile, print);
		}


		return;
	}

	void	Init::FileRead_Output(Output& output, const string& file, const bool print) {
		vector<string> mainWords;
		mainWords.push_back("Grid");
		mainWords.push_back("Temperature");
		mainWords.push_back("Solidification");
		mainWords.push_back("Solidification+");

		vector<vector<string>> subWords(mainWords.size());
		subWords[0].push_back("x");
		subWords[0].push_back("y");
		subWords[0].push_back("z");

		subWords[1].push_back("T");
		subWords[1].push_back("T_hist");

		subWords[2].push_back("tSol");
		subWords[2].push_back("G");
		subWords[2].push_back("Gx");
		subWords[2].push_back("Gy");
		subWords[2].push_back("Gz");
		subWords[2].push_back("V");
		subWords[2].push_back("dTdt");
		subWords[2].push_back("eqFrac");
		subWords[2].push_back("depth");
		subWords[2].push_back("numMelt");
		subWords[2].push_back("RDF");
		subWords[2].push_back("MP_Stats");

		subWords[3].push_back("H");
		subWords[3].push_back("Hx");
		subWords[3].push_back("Hy");
		subWords[3].push_back("Hz");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(output.x, values[0][0], true, "output-x", 1, print);
		Init::SetValues(output.y, values[0][1], true, "output-y", 1, print);
		Init::SetValues(output.z, values[0][2], true, "output-z", 1, print);

		Init::SetValues(output.T, values[1][0], true, "output-T", 1, print);
		Init::SetValues(output.T_hist, values[1][1], false, "output-T_hist", 1, print);

		Init::SetValues(output.tSol, values[2][0], false, "output-tSol", 1, print);
		Init::SetValues(output.G, values[2][1], false, "output-G", 1, print);
		Init::SetValues(output.Gx, values[2][2], false, "output-Gx", 1, print);
		Init::SetValues(output.Gy, values[2][3], false, "output-Gy", 1, print);
		Init::SetValues(output.Gz, values[2][4], false, "output-Gz", 1, print);
		Init::SetValues(output.V, values[2][5], false, "output-V", 1, print);
		Init::SetValues(output.dTdt, values[2][6], false, "output-dTdt", 1, print);
		Init::SetValues(output.eqFrac, values[2][7], false, "output-eqFrac", 1, print);
		Init::SetValues(output.depth, values[2][8], false, "output-depth", 1, print);
		Init::SetValues(output.numMelt, values[2][9], false, "output-numMelt", 1, print);
		Init::SetValues(output.RDF, values[2][10], false, "output-RDF", 1, print);
		Init::SetValues(output.mp_stats, values[2][11], false, "output-mpStats", 1, print);

		Init::SetValues(output.H, values[3][0], false, "output-H", 1, print);
		Init::SetValues(output.Hx, values[3][1], false, "output-Hx", 1, print);
		Init::SetValues(output.Hy, values[3][2], false, "output-Hy", 1, print);
		Init::SetValues(output.Hz, values[3][3], false, "output-Hz", 1, print);
		return;
	}
	void	Init::FileRead_Settings(Settings& settings, const string& file, const bool print) {

		vector<string> mainWords;

		mainWords.push_back("Temperature");
		mainWords.push_back("Path");
		mainWords.push_back("Compute");
		mainWords.push_back("MPI");

		vector<vector<string>> subWords(mainWords.size());
		
		subWords[0].push_back("Sol_Tol");
		subWords[0].push_back("Sol_Iter");
		subWords[0].push_back("Cutoff_Peak");
		subWords[0].push_back("Cutoff_T0TL");
		
		subWords[1].push_back("Buffer");
		subWords[1].push_back("Compression");

		subWords[2].push_back("MaxThreads");
		subWords[2].push_back("PINT");

		subWords[3].push_back("Overlap");

		vector<vector<string>> values(mainWords.size());

		Init::Keywords_Lv2(mainWords, subWords, values, file, print);

		Init::SetValues(settings.dttest, values[0][0], 1e-3, "Solidification Tolerance", 0, print);
		Init::SetValues(settings.max_iter, values[0][1], 10, "Solidification Iterations", 0, print);
		Init::SetValues(settings.t_hist, values[0][2], 1e-9, "Peak Temperature Cutoff Ratio", 0, print);
		Init::SetValues(settings.p_hist, values[0][3], 1e-2, "Solidfication to Initial Temperature Cutoff Ratio", 0, print);

		Init::SetValues(settings.r_max, values[1][0], -1.0, "Domain Buffer", 0, print);
		Init::SetValues(settings.compress, values[1][1], 0, "Path Compression", 0, print);

		Init::SetValues(settings.thnum, values[2][0], omp_get_max_threads()/2, "Number of Threads", 0, print);
		Init::SetValues(settings.use_PINT, values[2][1], 0, "Parallel in Time Mode", 0, print);

		Init::SetValues(settings.mpi_overlap, values[3][0], 0, "Mpi Overlap", 0, print);

		return;
	}

	void	Init::FileRead_Points(Domain& domain, const string& file, const bool print) {
		std::ifstream readFile;
		readFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		string line;
		int num_read = 0;
		try {
			readFile.open(file.c_str(), std::ios::in);
			coord temp;

			//Set initial position as 0, 0, 0
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;

			//Read in path information from file
			while (getline(readFile, line))
			{
				readFile >> temp.x >> temp.y >> temp.z;
				temp.x /= 1000.0; temp.y /= 1000.0; temp.z /= 1000.0;
				domain.points.push_back(temp); 
				num_read++;
			}
		}
		catch (const std::ifstream::failure&) {
			catchPrint(file, num_read, print);
			readFile.close();
		}
		catch (const std::exception& e) {
			catchPrint(file, num_read, print);
			readFile.close();
		}

		domain.pnum = num_read;
		if (num_read) { 
			if (print) std::cout << domain.pnum << " points found in " << file << std::endl;
		}
		else { 
			if (print) std::cout << "ERROR: No points found in " << file << std::endl;
		}

		return;
	}

	void	Init::SetDiffusivity(Material& material) {
		material.a = material.kon / (material.rho * material.cps);
	}
	void	Init::SetBeamPower(Beam& beam){
		if (beam.eff == DBL_MAX) { beam.eff = 1.0; }
		beam.q = beam.q * beam.eff * 2.0;
	}
	void	Init::SetDomainParams(Domain& domain) {

		// If using a custom point file, skip
		if (domain.customPoints) { return; }

		// Otherwise, set (x,y,z) parameters (num, res, etc.)
		if (domain.xnum == INT_MAX) {
			domain.xnum = 1 + int(0.5 + (domain.xmax - domain.xmin) / domain.xres);
			domain.xmax = domain.xmin + (domain.xnum - 1) * domain.xres;
		}
		else if (domain.xnum != 1) {
			domain.xres = (domain.xmax - domain.xmin) / (domain.xnum - 1);
		}
		else {
			domain.xres = DBL_MAX;
		}

		if (domain.ynum == INT_MAX) {
			domain.ynum = 1 + int(0.5 + (domain.ymax -domain.ymin) / domain.yres);
			domain.ymax = domain.ymin + (domain.ynum - 1) * domain.yres;
		}
		else if (domain.ynum != 1) {
			domain.yres = (domain.ymax - domain.ymin) / (domain.ynum - 1);
		}
		else {
			domain.yres = DBL_MAX;
		}

		if (domain.znum == INT_MAX) {
			domain.znum = 1 + int(0.5 + (domain.zmax - domain.zmin) / domain.zres);
			domain.zmax = domain.zmin + (domain.znum - 1) * domain.zres;
		}
		else if (domain.ynum != 1) {
			domain.zres = (domain.zmax - domain.zmin) / (domain.znum - 1);
		}
		else {
			domain.zres = DBL_MAX;
		}

		domain.pnum = domain.xnum * domain.ynum * domain.znum;
	}
}
