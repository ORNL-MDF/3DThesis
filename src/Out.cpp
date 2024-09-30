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
#include "Mpi_structs.hpp"

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