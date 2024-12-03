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

#pragma once
#include "impl/Structs/DataStructs.hpp"
#include "impl/Calc/Util.hpp"
#include "impl/Calc/Calc.hpp"

#include <iostream>
#include <fstream>

#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <functional> 

using std::bind;
using std::function;
using std::vector;
using std::string;
using std::to_string;
using std::placeholders::_1;

////// Modes (exclusive)
//// Snapshot
// Calculates temperature at all points at end of scan path
//
//// Temp History
// Tracks temperature history for select points
//
//// Solidificiation
// Calculates solidification conditions
//
//// Solidification+
// Calculates primary and secondary solidifiaction conditions
//
//// Reduced Data Format
// Calculates melting and solidification events
//

////// Heat Source 
//// Many Beams
// Usual Gaussian (BeamX)
// Single Path File Per Gaussian (PathX) 

namespace Thesis::impl
{
	class Grid {
	private:
		double xmin = -DBL_MAX;
		double xmax = DBL_MAX;
		double ymin = -DBL_MAX;
		double ymax = DBL_MAX;
		double zmin = -DBL_MAX;
		double zmax = DBL_MAX;
		
		uint16_t xnum = 0;
		uint16_t ynum = 0;
		uint16_t znum = 0;
		
		uint16_t* i = NULL; // index - always needed
		uint16_t* j = NULL; // index - always needed
		uint16_t* k = NULL; // index - always needed

		bool* T_calc_flag = NULL; // make sure temperature only calculated once - needed from any MP tracking mode
		bool* output_flag = NULL; // should a point be output - needed for MP tracking modes

		double* x = NULL; // coordinate - always needed
		double* y = NULL; // coordinate - always needed
		double* z = NULL; // coordinate - always needed

		double* T = NULL; // temperature - always needed
		double* T_last = NULL; // last temperature - needed for solidification times

		double* tSol = NULL;
		double* G = NULL;  // solidfication gradient - store if output
		double* V = NULL;  // solidfication velocity - store if output
		double* Gx = NULL; // solidification direction - store if output
		double* Gy = NULL; // solidification direction - store if output
		double* Gz = NULL; // solidification direction - store if output
		double* dTdt = NULL; // cooling rate - store if output
		double* eqFrac = NULL; // equiaxed fraction from CET - store if output
		uint16_t* numMelt = NULL;

		double* H = NULL; // secondary solidification magnitude - store if output
		double* Hx = NULL; // secondary solidification direction - store if output
		double* Hy = NULL; // secondary solidification direction - store if output
		double* Hz = NULL; // secondary solidification direction - store if output

		double* depth = NULL; // max depth under point - store if output

		vector<double>* T_hist = NULL; // temperature history - store if output
		vector<double>* t_hist = NULL; // iteration history - store if output

		vector<double>* RDF_tm = NULL; // melting time for reduced data format - store if output
		vector<double>* RDF_tl = NULL; // liquid time for reduced data format - store if output
		vector<double>* RDF_cr = NULL; // cooling rate for reduced data format - store if output 

		double* MP_Width = NULL; // maximum meltpool width - store if output
		double* MP_Length = NULL; // maximum meltpool length - store if output
		double* MP_Depth = NULL; // maximum meltpool depth - store if output

		vector<uint32_t> RRDF_idxs; // indices for doubly reduced data format - initilaize size if used
		vector<double> RRDF_ts;  // times for doubly reduced data format - initilaize size if used
		vector<double> RRDF_Ts; // temperatures for doubly reduced data format - initilaize size if used

		vector<string> outputNames; // what all to output ("x", "y", "z", etc.)
		vector<function<double(const int)>> outputFuncs; // the functions for outputting the variables

	public:
		Grid() {};
		Grid(Simdat& sim) {
			
			xmin = sim.domain.xmin; xmax = sim.domain.xmax; xnum = sim.domain.xnum;
			ymin = sim.domain.ymin; ymax = sim.domain.ymax; ynum = sim.domain.ynum;
			zmin = sim.domain.zmin; zmax = sim.domain.zmax; znum = sim.domain.znum;

			const int pnum = sim.domain.pnum;

			i = new uint16_t[pnum]();
			j = new uint16_t[pnum]();
			k = new uint16_t[pnum]();
			
			if (sim.domain.customPoints) {
				x = new double[pnum]();
				y = new double[pnum]();
				z = new double[pnum]();
			}

			if (sim.output.x) {
				outputNames.push_back("x");
				outputFuncs.push_back(bind(&Grid::get_x, this, _1));
			}
			
			if (sim.output.y) {
				outputNames.push_back("y");
				outputFuncs.push_back(bind(&Grid::get_y, this, _1));
			}

			if (sim.output.z) {
				outputNames.push_back("z");
				outputFuncs.push_back(bind(&Grid::get_z, this, _1));
			}

			T = new double[pnum]();
			if (sim.output.T) {
				outputNames.push_back("T");
				outputFuncs.push_back(bind(&Grid::get_T, this, _1));
			}
			if (sim.output.T_hist) {
				T_hist = new vector<double>[pnum];
				t_hist = new vector<double>[pnum];
			}

			InitializeGridPoints(sim);

			if (sim.param.mode == "Solidification" || sim.param.tracking != "None") {
				T_calc_flag = new bool[pnum]();
				output_flag = new bool[pnum]();
			}
			
			if (sim.param.mode == "Stork"){
				RRDF_idxs.reserve(pnum);
				RRDF_ts.reserve(2*pnum);
				RRDF_Ts.reserve(16*pnum);
			}

			if (sim.param.mode == "Solidification") {
				T_last = new double[pnum]();
				if (sim.output.tSol) {
					sim.util.do_sol = true;
					tSol = new double[pnum]();
					outputNames.push_back("tSol");
					outputFuncs.push_back(bind(&Grid::get_tSol, this, _1));
				}
				if (sim.output.G) { 
					sim.util.do_sol = true;
					G = new double[pnum]();
					outputNames.push_back("G");
					outputFuncs.push_back(bind(&Grid::get_G, this, _1));
				}
				if (sim.output.V) { 
					sim.util.do_sol = true;
					V = new double[pnum]();
					outputNames.push_back("V");
					outputFuncs.push_back(bind(&Grid::get_V, this, _1));
				}
				if (sim.output.Gx) { 
					sim.util.do_sol = true;
					Gx = new double[pnum]();
					outputNames.push_back("Gx");
					outputFuncs.push_back(bind(&Grid::get_Gx, this, _1));
				}
				if (sim.output.Gy) { 
					sim.util.do_sol = true;
					Gy = new double[pnum]();
					outputNames.push_back("Gy");
					outputFuncs.push_back(bind(&Grid::get_Gy, this, _1));
				}
				if (sim.output.Gz) { 
					sim.util.do_sol = true;
					Gz = new double[pnum]();
					outputNames.push_back("Gz");
					outputFuncs.push_back(bind(&Grid::get_Gz, this, _1));
				}
				if (sim.output.dTdt) { 
					sim.util.do_sol = true;
					dTdt = new double[pnum]();
					outputNames.push_back("dTdt");
					outputFuncs.push_back(bind(&Grid::get_dTdt, this, _1));
				}
				if (sim.output.eqFrac) { 
					sim.util.do_sol = true;
					eqFrac = new double[pnum]();
					outputNames.push_back("eqFrac");
					outputFuncs.push_back(bind(&Grid::get_eqFrac, this, _1));
				}
				if (sim.output.depth && sim.param.tracking == "Surface") {
					depth = new double[pnum]();
					outputNames.push_back("depth");
					outputFuncs.push_back(bind(&Grid::get_depth, this, _1));
				}
				if (sim.output.numMelt) {
					if (!sim.output.RDF) { numMelt = new uint16_t[pnum](); }
					outputNames.push_back("numMelt");
					outputFuncs.push_back(bind(&Grid::get_numMelt, this, _1));
				}
				if (sim.output.RDF) {
					sim.util.do_sol = true;
					RDF_tm = new vector<double>[pnum]();
					RDF_tl = new vector<double>[pnum]();
					RDF_cr = new vector<double>[pnum]();
				}
				if (sim.output.mp_stats){
					MP_Width = new double[pnum]();
					MP_Length = new double[pnum]();
					MP_Depth = new double[pnum]();
					
					outputNames.push_back("MP_width");
					outputNames.push_back("MP_length");
					outputNames.push_back("MP_depth");

					outputFuncs.push_back(bind(&Grid::get_mpWidth, this, _1));
					outputFuncs.push_back(bind(&Grid::get_mpLength, this, _1));
					outputFuncs.push_back(bind(&Grid::get_mpDepth, this, _1));	
				}
			}

			if (sim.param.mode == "Solidification" && sim.param.secondary == true) {
				if (sim.output.H) { 
					H = new double[pnum]; 
					outputNames.push_back("H");
					outputFuncs.push_back(bind(&Grid::get_H, this, _1));
				}
				if (sim.output.Hx) { 
					Hx = new double[pnum];
					outputNames.push_back("Hx");
					outputFuncs.push_back(bind(&Grid::get_Hx, this, _1));
				}
				if (sim.output.Hy) { 
					Hy = new double[pnum]; 
					outputNames.push_back("Hy");
					outputFuncs.push_back(bind(&Grid::get_Hy, this, _1));
				}
				if (sim.output.Hz) { 
					Hz = new double[pnum]; 
					outputNames.push_back("Hz");
					outputFuncs.push_back(bind(&Grid::get_Hz, this, _1));
				}
			}

		}
		~Grid() {
			if (i != NULL) { delete[] i; i = NULL; }
			if (j != NULL) { delete[] j; j = NULL;}
			if (k != NULL) { delete[] k; k = NULL;}
			if (T_calc_flag != NULL) { delete[] T_calc_flag; T_calc_flag = NULL;}
			if (output_flag != NULL) { delete[] output_flag; output_flag = NULL;}
			if (x != NULL) { delete[] x; x = NULL;}
			if (y != NULL) { delete[] y; y = NULL;}
			if (z != NULL) { delete[] z; z = NULL;}
			if (T != NULL) { delete[] T; T = NULL;}
			if (T_last != NULL) { delete[] T_last; T_last = NULL;}
			if (tSol != NULL) { delete[] tSol; tSol = NULL;}
			if (G != NULL) { delete[] G; G = NULL;}
			if (V != NULL) { delete[] V; V = NULL;}
			if (Gx != NULL) { delete[] Gx; Gx = NULL;}
			if (Gy != NULL) { delete[] Gy; Gy = NULL;}
			if (Gz != NULL) { delete[] Gz; Gz = NULL;}
			if (dTdt != NULL) { delete[] dTdt; dTdt = NULL;}
			if (eqFrac != NULL) { delete[] eqFrac; eqFrac = NULL;}
			if (numMelt != NULL) { delete[] numMelt; numMelt = NULL;}
			if (H != NULL) { delete[] H; H = NULL;}
			if (Hx != NULL) { delete[] Hx; Hx = NULL;}
			if (Hy != NULL) { delete[] Hy; Hy = NULL;}
			if (Hz != NULL) { delete[] Hz; Hz = NULL;}
			if (depth != NULL) { delete[] depth; depth = NULL;}
			if (T_hist != NULL) { delete[] T_hist; T_hist = NULL;}
			if (t_hist != NULL) { delete[] t_hist; t_hist = NULL;}
			if (RDF_tm != NULL) { delete[] RDF_tm; RDF_tm = NULL;}
			if (RDF_tl != NULL) { delete[] RDF_tl; RDF_tl = NULL;}
			if (RDF_cr != NULL) { delete[] RDF_cr; RDF_cr = NULL;}
			if (MP_Width != NULL) { delete[] MP_Width; MP_Width = NULL;}
			if (MP_Length != NULL) { delete[] MP_Length; MP_Length = NULL;}
			if (MP_Depth != NULL) { delete[] MP_Depth; MP_Depth = NULL;}
		};
		
		void InitializeGridPoints(const Simdat&);
		void Output(const Simdat&, const string);
		vector<vector<double>> Output_Table(const Simdat& sim);
		void Output_T_hist(const Simdat&, const string);
		void Output_RDF(const Simdat&, const string);
		void Output_RRDF_csv(const Simdat&, const string);
		void Output_RRDF_bin(const Simdat&, const string);

		double Calc_T(const double, const Nodes&, const Simdat&, const bool, const int);				// Returns temperature 
		void Solidify(vector<int>&, const double, const Simdat&, const int);												// Solidifies point
		double Calc_Solidification_time(vector<int>&, const double, const Simdat&, const int);							// Returns time of solidification
		vector<vector<double>> Calc_Solidficiaton_Primary(const double, const Nodes&, const int);	// Returns {Gx,Gy,Gz,Laplace,dT_t}
		vector<vector<double>> Calc_Solidficiaton_Secondary(const double, const Nodes&, const int);
		void Set_Solidficiaton_Primary(const vector<double>&, const Simdat&, const int);
		void Set_Solidficiaton_Secondary(const vector<double>&, const vector<double>&, const Simdat&, const int);

		// Gets //
		inline int get_i(const int p) { return i[p]; }
		inline int get_j(const int p) { return j[p]; }
		inline int get_k(const int p) { return k[p]; }

		bool get_T_calc_flag(const int p) { return T_calc_flag[p]; }
		bool get_output_flag(const int p) { 
			bool out = true;
			if (output_flag != NULL){ out = output_flag[p]; }
			return out;
		}
		
		double get_x(const int p) {
			double xcoord = xmax;
			if (x != NULL) { xcoord = x[p]; }
			else if (xnum - 1) { xcoord = xmin + (double(i[p]) * ((xmax - xmin) / (xnum - 1)));  }
			return xcoord;
		}
		double get_y(const int p) { 
			double ycoord = ymax;
			if (y != NULL) { ycoord = y[p]; }
			else if (ynum - 1) { ycoord = ymin + (double(j[p]) * ((ymax - ymin) / (ynum - 1)));  }
			return ycoord;
		}
		double get_z(const int p) {
			double zcoord = zmax;
			if (z != NULL) { zcoord = z[p]; }
			else if (znum - 1) { zcoord = zmin + (double(k[p]) * ((zmax - zmin) / (znum - 1))); }
			return zcoord; 
		}

		double get_T(const int p) { return T[p]; }
		double get_T_last(const int p) { return T_last[p]; }

		double get_tSol(const int p) { return tSol[p]; }
		double get_G(const int p) { return G[p]; }
		double get_Gx(const int p) { return Gx[p]; }
		double get_Gy(const int p) { return Gy[p]; }
		double get_Gz(const int p) { return Gz[p]; }
		double get_V(const int p) { return V[p]; }
		double get_dTdt(const int p) { return dTdt[p]; }
		double get_eqFrac(const int p) { return eqFrac[p]; }
		double get_depth(const int p) { 
			return (depth[p] * ((zmax - zmin) / (znum - 1))); 
		}
		uint16_t get_numMelt(const int p) { 
			if (numMelt != NULL) {
				return numMelt[p];
			}
			else {
				return get_RDF_tm(p).size();
			}
		}

		double get_H(const int p) { return G[p]; }
		double get_Hx(const int p) { return Hx[p]; }
		double get_Hy(const int p) { return Hy[p]; }
		double get_Hz(const int p) { return Hz[p]; }

		vector<double> get_T_hist(const int p) { return T_hist[p]; }
		vector<double> get_t_hist(const int p) { return t_hist[p]; }

		vector<double> get_RDF_tm(const int p) { return RDF_tm[p]; }
		vector<double> get_RDF_tl(const int p) { return RDF_tl[p]; }
		vector<double> get_RDF_cr(const int p) { return RDF_cr[p]; }

		double get_mpLength(const int p) { return MP_Length[p]; }
		double get_mpWidth(const int p) { return MP_Width[p]; }
		double get_mpDepth(const int p) { return MP_Depth[p]; }

		vector<uint32_t>& get_RRDF_idxs() { return RRDF_idxs; }
		vector<double>& get_RRDF_ts() { return RRDF_ts; }
		vector<double>& get_RRDF_Ts() { return RRDF_Ts; }

		// Sets //
		void set_T_calc_flag(const bool b, const int p) { if (T_calc_flag != NULL) { T_calc_flag[p] = b; } }
		void set_output_flag(const bool b, const int p) { if (output_flag != NULL) { output_flag[p] = b; } }

		void set_T(const double d, const int p) { if (T != NULL) { T[p] = d; } }
		void add_T_hist(const double d, const int p) { if (T_hist != NULL) { T_hist[p].push_back(d); } }
		void add_t_hist(const double d, const int p) { if (t_hist != NULL) { t_hist[p].push_back(d); } }
		
		void set_T_last(const int p) {if (T_last != NULL) { T_last[p] = T[p]; }}
		void set_T_last(const double d, const int p) { if (T_last != NULL) { T_last[p] = d; } }

		void set_tSol(const double d, const int p) { if (tSol != NULL) { tSol[p] = d; } }
		void set_G(const double d, const int p) { if (G != NULL) { G[p] = d; } }
		void set_Gx(const double d, const int p) { if (Gx != NULL) { Gx[p] = d; } }
		void set_Gy(const double d, const int p) { if (Gy != NULL) { Gy[p] = d; } }
		void set_Gz(const double d, const int p) { if (Gz != NULL) { Gz[p] = d; } }
		void set_V(const double d, const int p) { if (V != NULL) { V[p] = d; } }
		void set_dTdt(const double d, const int p) { if (dTdt != NULL) { dTdt[p] = d; } }
		void set_eqFrac(const double d, const int p) { if (eqFrac != NULL) { eqFrac[p] = d; } }
		void add_numMelt(const int p) { if (numMelt != NULL) { numMelt[p] += 1; } }

		void set_H(const double d, const int p) { if (H != NULL) { H[p] = d; } }
		void set_Hx(const double d, const int p) { if (Hx != NULL) { Hx[p] = d; } }
		void set_Hy(const double d, const int p) { if (Hy != NULL) { Hy[p] = d; } }
		void set_Hz(const double d, const int p) { if (Hz != NULL) { Hz[p] = d; } }

		void set_depth(const double d, const int p) { if (depth != NULL) { depth[p] = d; } }

		void add_RDF_tm(const double d, const int p) { 
			if (RDF_tm != NULL) {
				if (RDF_tm[p].size() == RDF_tl[p].size()) {RDF_tm[p].push_back(d);}
			}
		}
		void add_RDF_tl(const double d, const int p) { 
			if (RDF_tl != NULL) { RDF_tl[p].push_back(d); } 
		}
		void add_RDF_cr(const double d, const int p) {
			if (RDF_cr != NULL) { RDF_cr[p].push_back(d); } 
		}

		void set_mpLength(const double length, const int p) { MP_Length[p] = std::max(length, MP_Length[p]); }
		void set_mpWidth(const double width, const int p) { MP_Width[p] = std::max(width, MP_Width[p]); }
		void set_mpDepth(const double depth, const int p) { MP_Depth[p] = std::max(depth, MP_Depth[p]); }
	};
}