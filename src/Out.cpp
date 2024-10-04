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

void Out::RRDF_csv(const Simdat& sim, const vector<uint32_t>& idxs, const vector<float>& ts, const vector<float>& Ts){
	
	// Make fileName
	const string csvFile = sim.files.dataDir + "/ThesisData.csv";

	// Get header info
	const uint32_t size = ts.size()/2;
	const uint32_t extent[3] = {static_cast<uint32_t>(sim.domain.xnum), static_cast<uint32_t>(sim.domain.ynum), static_cast<uint32_t>(sim.domain.znum)};
	const float floats[5] = {static_cast<float>(sim.domain.xmin), static_cast<float>(sim.domain.ymin), static_cast<float>(sim.domain.zmin), static_cast<float>(sim.domain.xres), static_cast<float>(sim.material.T_liq)};

	// Output CSV
	std::ofstream datafile;
	datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try {
		datafile.open(csvFile.c_str());
		datafile << "x,y,z,t_prev,t_cur";
        datafile << ",T000_prev,T001_prev,T010_prev,T011_prev,T100_prev,T101_prev,T110_prev,T111_prev";
        datafile << ",T000_cur,T001_cur,T010_cur,T011_cur,T100_cur,T101_cur,T110_cur,T111_cur\n";
		for (uint32_t p = 0; p < size; p++) {
            datafile << idxs[3*p+0] << "," << idxs[3*p+1] << "," << idxs[3*p+2] << "," << ts[2*p+0] << "," << ts[2*p+1];
            for (int n=0;n<8;n++){
                datafile << "," << Ts[16*p+n];
            }
            for (int n=0;n<8;n++){
                datafile << "," << Ts[16*p+8+n];
            }
            datafile << "\n";
		}
	}
	catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
	datafile.close();
}

void Out::RRDF_bin(const Simdat& sim, const vector<uint32_t>& idxs, const vector<float>& ts, const vector<float>& Ts){
	
	// Make file name
	const string binFile = sim.files.dataDir + "/ThesisData.bin";
	
	// Get header info
	const uint32_t size = ts.size()/2;
	const uint32_t extent[3] = {static_cast<uint32_t>(sim.domain.xnum), static_cast<uint32_t>(sim.domain.ynum), static_cast<uint32_t>(sim.domain.znum)};
	const float floats[5] = {static_cast<float>(sim.domain.xmin), static_cast<float>(sim.domain.ymin), static_cast<float>(sim.domain.zmin), static_cast<float>(sim.domain.xres), static_cast<float>(sim.material.T_liq)};	
	
	// Open up binary file
	std::ofstream os(binFile, std::ios::binary);
	
	// Write header to binary file	
	os.write((const char*)&size, sizeof(uint32_t));
	os.write((const char*)&extent, 3 * sizeof(uint32_t));
	os.write((const char*)&floats, 5 * sizeof(float));
	
	// Write vectors to binary file
	os.write((const char*)&idxs[0], 3 * size * sizeof(uint32_t));
	os.write((const char*)&ts[0], 2 * size * sizeof(float));
	os.write((const char*)&Ts[0], 16 * size * sizeof(float));
	
	// Close binary file
	os.close();

	std::cout << "Sizes: " << size << "\t" << idxs.size() << "\t" << ts.size() << "\t" << Ts.size() << std::endl;
}