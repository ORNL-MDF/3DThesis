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

#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>

#include "DataStructs.h"
#include "Point.h"
#include "Init.h"
#include "Util.h"
#include "Run.h"
#include "Out.h"

using std::vector;
using std::string;

int main(int argc, char * argv[]) {	
	
	int start_s = clock();
	
	//Command line input to program is the file that contains the simulation file names
	string	input;
	if (argc <= 1) { input = "TestInputs/ParamInput.txt"; }
	else { input = argv[1]; }	

	Simdat sim;
	Init::GetFileNames(sim,input);	//Get names of intput files

	vector<path_seg> segv;
	Init::ReadSimParams(segv, sim);		//Read simulation parameters from input files

	Point* const ptv = new Point[sim.param.pnum];
	Init::SetPoints(ptv,sim);		//Set up points

	vector<int> seg_num;
	Util::InitStartSeg(seg_num, segv, sim); 	// Better Search Algorythm Setup

	Run::Simulate(ptv, segv, sim, seg_num);

	int stop_s = clock();
	std::cout << "Execution time (s): " << (stop_s - start_s) / double(CLOCKS_PER_SEC) << "\n\n";
	
	if (sim.param.mode){ Out::Write_csv(ptv, sim, "Final", sim.setting.out_mode); }
	else { Out::Write_csv(ptv, sim, "Snapshot", 0); }

	if (sim.setting.T_hist) { Out::Write_T_hist(ptv, sim, "TemperatureHistory"); }
	//system("pause");

	return 0;
}
