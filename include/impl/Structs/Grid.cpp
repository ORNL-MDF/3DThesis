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

#include "impl/Structs/Grid.hpp"
#include "impl/IO/Out.hpp"
#include "impl/Calc/Util.hpp"

namespace Thesis::impl
{
	void Grid::InitializeGridPoints(const Simdat& sim) {
		if (sim.domain.customPoints) {
			#pragma omp for schedule(static)
			for (int p = 0; p < sim.domain.pnum; p++) {
				i[p] = p;
				j[p] = 0;
				k[p] = 0;

				x[p] = sim.domain.points[p].x;
				y[p] = sim.domain.points[p].y;
				z[p] = sim.domain.points[p].z;
			}
		}
		else {
			#pragma omp parallel num_threads(sim.settings.thnum)
			{
				double xp, yp, zp;
				int p;
				#pragma omp for schedule(static)
				for (int ii = 0; ii < sim.domain.xnum; ii++) {
					if (sim.domain.xnum == 1) { xp = sim.domain.xmax; }
					else { xp = sim.domain.xmin + ((double)ii * ((sim.domain.xmax - sim.domain.xmin) / ((double)sim.domain.xnum - 1))); }
					for (int jj = 0; jj < sim.domain.ynum; jj++) {
						if (sim.domain.ynum == 1) { yp = sim.domain.ymax; }
						else { yp = sim.domain.ymin + ((double)jj * ((sim.domain.ymax - sim.domain.ymin) / ((double)sim.domain.ynum - 1))); }
						for (int kk = 0; kk < sim.domain.znum; kk++) {
							if (sim.domain.znum == 1) { zp = sim.domain.zmax; }
							else { zp = sim.domain.zmin + ((double)kk * ((sim.domain.zmax - sim.domain.zmin) / ((double)sim.domain.znum - 1))); }
							p = Util::ijk_to_p(ii, jj, kk, sim);		
							i[p] = ii; j[p] = jj; k[p] = kk;
							//x[p] = xp; y[p] = yp; z[p] = zp;
						}
					}
				}
			}
		}
	}

	vector<vector<double>> Grid::Output_Table(const Simdat& sim) {
		vector<vector<double>> data;
		vector<double> row;
		// If, for some reason, nothing is being output: return
		const int numOut = outputNames.size();
		if (numOut==0) { return data;}
		for (int p = 0; p < sim.domain.pnum; p++){
			if (get_output_flag(p)){
				for (int i = 0; i < numOut; i++){
					row.push_back(outputFuncs[i](p));
				}
			}
			data.push_back(row);
			row.clear();
		}
		return data;
	}

	void Grid::Output(const Simdat& sim, const string name){
		const int numOut = outputNames.size();
		if (numOut==0) { return; }

		std::ofstream datafile;
		datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		string out_file = sim.files.dataDir + "/" + sim.files.name + "." + name + ".csv";
		try {
			datafile.open(out_file.c_str());

			datafile << outputNames[0];
			for (int i = 1; i < numOut; i++) {datafile << "," << outputNames[i];}
			datafile << "\n";
			
			for (int p = 0; p < sim.domain.pnum; p++) {
				if (get_output_flag(p)) {
					datafile << outputFuncs[0](p);
					for (int i = 1; i < numOut; i++) { datafile << "," << outputFuncs[i](p); }
					datafile << "\n";
				}
			}

		}
		catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
		datafile.close();
		return;

	}

	void Grid::Output_T_hist(const Simdat& sim, const string name) {
		// If T_hist isn't ouput: return
		if (T_hist == NULL) { return; }

		// If the points have the same times
		if (sim.param.tracking=="None"){
			// Create datafile with list of points and their coordinates
			std::ofstream datafile;
			datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
			string out_file = sim.files.dataDir + "/" + sim.files.name + name + ".csv";
			try {
				datafile.open(out_file.c_str());
				datafile << "x,y,z";
				const int max_tnum = get_t_hist(0).size();
				const int t_width = 1+static_cast<int>(1+log10(max_tnum));
				// for (int tnum=0;tnum<get_t_hist(0).size();tnum++){"t"+Util::ZeroPadNumber(tnum,t_width);}
				for (int tnum=0;tnum<max_tnum;tnum++){
					datafile<< ","+std::to_string(get_t_hist(0)[tnum]);
				}
				// Newline
				datafile << "\n";
				// Output point information
				for (int p = 0; p < sim.domain.pnum; p++) {
					// x-y-z info
					datafile << get_x(p) << "," << get_y(p) << "," << get_z(p);
					// Temperatures 
					for (int tnum=0;tnum<max_tnum;tnum++){
						datafile<< ","+std::to_string(get_T_hist(p)[tnum]);
					}
					// Newline
					datafile <<  "\n";
				}
			}
			catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
			datafile.close();
		}
		// If the points have different times
		else{
			// Create datafile with list of points and their coordinates
			std::ofstream datafile;
			datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
			string out_file = sim.files.dataDir + "/" + sim.files.name + ".Points.csv";
			try {
				datafile.open(out_file.c_str());
				datafile << "p,x,y,z\n";
				for (int p = 0; p < sim.domain.pnum; p++) {
					if (get_output_flag(p)) {
						datafile << p << "," << get_x(p) << "," << get_y(p) << "," << get_z(p) << "\n";
					}
				}
			}
			catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
			datafile.close();

			// Write files for each point
			for (int p = 0; p < sim.domain.pnum; p++) {
				if (get_output_flag(p)) {
					out_file = sim.files.dataDir + "/" + sim.files.name + ".T_hist." + to_string(p) + ".csv";
					try {
						datafile.open(out_file.c_str());
						datafile << "t,T\n";
						for (int e = 0; e < get_t_hist(p).size(); e++) {
							{
								datafile << get_t_hist(p)[e] << "," << get_T_hist(p)[e] << "\n";
							}
						}
					}
					catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
					datafile.close();
				}
			}
		}

	}

	void Grid::Output_RDF(const Simdat& sim, const string name) {
		
		if (RDF_tm == NULL) { return; }

		std::ofstream datafile;
		datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		string out_file = sim.files.dataDir + "/" + sim.files.name + "." + name + ".csv";
		try {
			datafile.open(out_file.c_str());

			datafile << "x,y,z,tm,tl,cr\n";

			for (int p = 0; p < sim.domain.pnum; p++) {
				if (get_output_flag(p)) {
					for (int e = 0; e < get_RDF_tm(p).size(); e++) {
						datafile << get_x(p) << "," << get_y(p) << "," << get_z(p) << ",";
						datafile << get_RDF_tm(p)[e] << "," << get_RDF_tl(p)[e] << "," << get_RDF_cr(p)[e] << "\n";
					}
				}
			}

		}
		catch (const std::ofstream::failure& e) { std::cout << "Exception writing data file, check that Data directory exists\n"; }
		datafile.close();
		return;
	}

	void Grid::Output_RRDF_csv(const Simdat& sim, const string name){
		
		// Make fileName
		const string csvFile = sim.files.dataDir + "/" + sim.files.name + "." + name + ".csv";

		// Set references to data
		const vector<uint32_t>& idxs = RRDF_idxs;
		const vector<double>& ts = RRDF_ts;
		const vector<double>& Ts = RRDF_Ts;

		// Get header info
		const uint32_t size = ts.size()/2;
		const uint32_t extent[3] = {static_cast<uint32_t>(sim.domain.xnum), static_cast<uint32_t>(sim.domain.ynum), static_cast<uint32_t>(sim.domain.znum)};
		const double doubles[5] = {static_cast<double>(sim.domain.xmin), static_cast<double>(sim.domain.ymin), static_cast<double>(sim.domain.zmin), static_cast<double>(sim.domain.xres), static_cast<double>(sim.material.T_liq)};

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

	void Grid::Output_RRDF_bin(const Simdat& sim, const string name){
		
		// Make file name
		const string binFile = sim.files.dataDir + "/" + sim.files.name + "." + name + ".mdf";
		
		// Set references to data
		const vector<uint32_t>& idxs = RRDF_idxs;
		const vector<double>& ts = RRDF_ts;
		const vector<double>& Ts = RRDF_Ts;

		// Get header info
		const uint32_t size = ts.size()/2;
		const uint32_t extent[3] = {static_cast<uint32_t>(sim.domain.xnum), static_cast<uint32_t>(sim.domain.ynum), static_cast<uint32_t>(sim.domain.znum)};
		const double doubles[5] = {static_cast<double>(sim.domain.xmin), static_cast<double>(sim.domain.ymin), static_cast<double>(sim.domain.zmin), static_cast<double>(sim.domain.xres), static_cast<double>(sim.material.T_liq)};	
		
		// Open up binary file
		std::ofstream os(binFile, std::ios::binary);
		
		// Write header to binary file	
		os.write((const char*)&size, sizeof(uint32_t));
		os.write((const char*)&extent, 3 * sizeof(uint32_t));
		os.write((const char*)&doubles, 5 * sizeof(double));
		
		// Write vectors to binary file
		os.write((const char*)&idxs[0], 3 * size * sizeof(uint32_t));
		os.write((const char*)&ts[0], 2 * size * sizeof(double));
		os.write((const char*)&Ts[0], 16 * size * sizeof(double));
		
		// Close binary file
		os.close();

		std::cout << "Sizes: " << size << "\t" << idxs.size() << "\t" << ts.size() << "\t" << Ts.size() << std::endl;
	}

	double Grid::Calc_T(const double t, const Nodes& nodes, const Simdat& sim, const bool set, const int p) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/

		const double xp = get_x(p);
		const double yp = get_y(p);
		const double zp = get_z(p);
		
		double dT = 0;
		for (size_t iter = 0; iter < nodes.size; iter++) {
			
			const double dx = xp - nodes.xb[iter];
			const double dy = yp - nodes.yb[iter];
			const double dz = zp - nodes.zb[iter];

			const double phi = exp(-3.0 * ((dx * dx * nodes.phix[iter]) + (dy * dy * nodes.phiy[iter]) + (dz * dz * nodes.phiz[iter])) + nodes.expmod[iter]);
			const double dT_seg = nodes.dtau[iter] * phi;
			
			dT += dT_seg;
		}

		const double T_temp = sim.material.T_init + dT;
		
		if (set) { 
			set_T(T_temp, p);
			add_T_hist(T_temp, p);
			add_t_hist(t, p);
			set_T_calc_flag(true, p);
			if (T_temp >= sim.material.T_liq) {
				set_G(DBL_MIN, p);
				set_V(DBL_MIN, p);
				set_eqFrac(-1.0, p);
				set_output_flag(1, p);
			}
		}

		return T_temp;
	}

	void Grid::Solidify(vector<int>& start_seg, const double t, const Simdat& sim, const int p) {
		add_numMelt(p);
		if (sim.util.do_sol == false) { return; }

		const double t_sol = Calc_Solidification_time(start_seg, t, sim, p);
		set_tSol(t_sol, p);
		add_RDF_tl(t_sol, p);


		Nodes nodes;
		Calc::Integrate_Serial(nodes, start_seg, sim, t_sol, true);

		if (!sim.param.secondary) {
			const vector<vector<double>> params = Calc_Solidficiaton_Primary(t_sol, nodes, p);
			const vector<double> primaryParams = params[0];
			Set_Solidficiaton_Primary(primaryParams, sim, p);
		}
		else {
			const vector<vector<double>>	params = Calc_Solidficiaton_Secondary(t_sol, nodes, p);
			const vector<double> primaryParams = params[0];
			const vector<double> secondaryParams = params[0];
			Set_Solidficiaton_Secondary(primaryParams, secondaryParams, sim, p);
		}
	}

	double Grid::Calc_Solidification_time(vector<int>& start_seg, const double t, const Simdat& sim, const int p) {

		// Controls iteration
		bool runFlag = true;
		int runIter = 0;;
		int maxIter = sim.settings.max_iter;
		
		// Temperatures
		double T2 = T[p];
		double T1 = T_last[p];
		double T0 = sim.material.T_liq;

		// Times
		double t2 = t;
		double t1 = t - sim.param.dt;
		double t0 = 0.0;

		double T_temp = 0.0;
		double T_err = 0.0;
		double err_limit = sim.settings.dttest;

		double m = 0.0;

		Nodes nodes;

		
		while (runFlag) { //Continue until you get a temperature really close to the solidification temp
			// Slope between two temperatures
			m = (T2 - T1) / (t2 - t1);

			// Solve for when the line intersects T0
			t0 = t1 + ((T0 - T1) / m);

			// Get Quadrature information at t0
			Util::ClearNodes(nodes);
			Calc::Integrate_Serial(nodes, start_seg, sim, t0, 0);

			// Calculate Temperature at t0
			T_temp = Calc_T(t0, nodes, sim, 0, p);

			// How far is this temperature from T0?
			T_err = abs(1.0 - T_temp / T0);
			
			// Increment interation
			runIter++;

			// If error below error limit or iteration above iteration limit, stop iterating
			if ((T_err < err_limit) || (runIter >= maxIter)) {
				runFlag = false;
			}
			// Otherwise, use new point in root finding method 
			else {
				if (T_temp > T0) {
					t1 = t0;
					T1 = T_temp;
				}
				else {
					t2 = t0;
					T2 = T_temp;
				}
			}
		}

		return t0;
	}

	vector<vector<double>> Grid::Calc_Solidficiaton_Primary(const double t, const Nodes& nodes, const int p) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/
		
		const double xp = get_x(p);
		const double yp = get_y(p);
		const double zp = get_z(p);

		double dT = 0;
		double Gx_temp = 0;
		double Gy_temp = 0;
		double Gz_temp = 0;
		double Laplace = 0;
		double dT_t = 0;
		for (size_t iter = 0; iter < nodes.size; iter++) {

			const double dx = xp - nodes.xb[iter];
			const double dy = yp - nodes.yb[iter];
			const double dz = zp - nodes.zb[iter];

			const double phi = exp(-3.0 * ((dx * dx * nodes.phix[iter]) + (dy * dy * nodes.phiy[iter]) + (dz * dz * nodes.phiz[iter])) + nodes.expmod[iter]);

			const double dT_seg = nodes.dtau[iter] * phi;

			dT += dT_seg;
			
			const double ddpx = (-6.0 * nodes.phix[iter]);
			const double ddpy = (-6.0 * nodes.phiy[iter]);
			const double ddpz = (-6.0* nodes.phiz[iter]);

			const double dpx = (ddpx * dx);
			const double dpy = (ddpy * dy);
			const double dpz = (ddpz * dz);

			Gx_temp += dT_seg * dpx;	//x-gradient
			Gy_temp += dT_seg * dpy;	//y-gradient
			Gz_temp += dT_seg * dpz;	//z-gradient

			Laplace += dT_seg * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
			dT_t += (nodes.dtau[iter]== 0) ? phi : 0;  //Notice that the interval of time we are integrating over is of size, that is because we are looking the instantaneous change at that time
		}

		vector<double> primaryParams = { Gx_temp, Gy_temp, Gz_temp, Laplace, dT_t };

		return { primaryParams };
	}

	vector<vector<double>> Grid::Calc_Solidficiaton_Secondary(const double t, const Nodes& nodes, const int p) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/
		
		const double xp = get_x(p);
		const double yp = get_y(p);
		const double zp = get_z(p);

		double dT = 0;
		double Gx_temp = 0;
		double Gy_temp = 0;
		double Gz_temp = 0;
		double Laplace = 0;
		double dT_t = 0;

		double dGxdx = 0;
		double dGxdy = 0;
		double dGxdz = 0;
		double dGydy = 0;
		double dGydz = 0;
		double dGzdz = 0;

		for (size_t iter = 0; iter < nodes.size; iter++) {

			const double dx = xp - nodes.xb[iter];
			const double dy = yp - nodes.yb[iter];
			const double dz = zp - nodes.zb[iter];

			const double phi = exp(-3.0 * ((dx * dx * nodes.phix[iter]) + (dy * dy * nodes.phiy[iter]) + (dz * dz * nodes.phiz[iter])) + nodes.expmod[iter]);

			const double dT_seg = nodes.dtau[iter] * phi;

			dT += dT_seg;
			
			const double ddpx = (-6.0 * nodes.phix[iter]);
			const double ddpy = (-6.0 * nodes.phiy[iter]);
			const double ddpz = (-6.0* nodes.phiz[iter]);

			const double dpx = (ddpx * dx);
			const double dpy = (ddpy * dy);
			const double dpz = (ddpz * dz);

			Gx_temp += dT_seg * dpx;	//x-gradient
			Gy_temp += dT_seg * dpy;	//y-gradient
			Gz_temp += dT_seg * dpz;	//z-gradient

			Laplace += dT_seg * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
			dT_t += (nodes.dtau[iter]== 0) ? phi : 0;  //Notice that the interval of time we are integrating over is of size, that is because we are looking the instantaneous change at that time
			
			dGxdx += dT_seg * (dpx * dpx + ddpx);
			dGxdy += dT_seg * (dpx * dpy);
			dGxdz += dT_seg * (dpx * dpz);

			dGydy += dT_seg * (dpy * dpy + ddpy);
			dGydz += dT_seg * (dpy * dpz);

			dGzdz += dT_seg * (dpz * dpz + ddpz);
		}

		vector<double> primaryParams = { Gx_temp, Gy_temp, Gz_temp, Laplace, dT_t };
		vector<double> secondaryParams = { dGxdx,dGxdy,dGxdz,dGydy,dGydz,dGzdz };

		return { primaryParams,secondaryParams };
	}

	void Grid::Set_Solidficiaton_Primary(const vector<double>& primaryParams, const Simdat& sim, const int p) {
		const double Gx_temp = primaryParams[0];
		const double Gy_temp = primaryParams[1];
		const double Gz_temp = primaryParams[2];
		const double Laplace = primaryParams[3];
		const double dT_t = primaryParams[4];
		
		const double G_temp = sqrt(Gx_temp * Gx_temp + Gy_temp * Gy_temp + Gz_temp * Gz_temp);

		const double Gxu_temp = Gx_temp / G_temp;
		const double Gyu_temp = Gy_temp / G_temp;
		const double Gzu_temp = Gz_temp / G_temp;

		const double dTdt_temp = abs(sim.material.a * Laplace + dT_t);
		const double V_temp = dTdt_temp / G_temp;

		const double eqFrac_temp = (1.0 - exp((-4 * PI * sim.material.cet_N0 / 3) * pow(G_temp * (sim.material.cet_n + 1) / (pow(sim.material.cet_a * V_temp, 1 / sim.material.cet_n)), -3)));

		set_G(G_temp, p);
		set_Gx(Gxu_temp, p);
		set_Gy(Gyu_temp, p);
		set_Gz(Gzu_temp, p);
		set_dTdt(dTdt_temp, p);
		set_V(V_temp, p);
		set_eqFrac(eqFrac_temp, p);
		add_RDF_cr(dTdt_temp, p);
	}

	void Grid::Set_Solidficiaton_Secondary(const vector<double>& primaryParams, const vector<double>& secondaryParams, const Simdat& sim, const int p) {
		const double Gx_temp = primaryParams[0];
		const double Gy_temp = primaryParams[1];
		const double Gz_temp = primaryParams[2];
		const double Laplace = primaryParams[3];
		const double dT_t = primaryParams[4];

		const double dGxdx = secondaryParams[0];
		const double dGxdy = secondaryParams[1];
		const double dGxdz = secondaryParams[2];
		const double dGydy = secondaryParams[3];
		const double dGydz = secondaryParams[4];
		const double dGzdz = secondaryParams[5];
		
		const double G_temp = sqrt(Gx_temp * Gx_temp + Gy_temp * Gy_temp + Gz_temp * Gz_temp);

		const double Gxu_temp = Gx_temp / G_temp;
		const double Gyu_temp = Gy_temp / G_temp;
		const double Gzu_temp = Gz_temp / G_temp;

		const double dTdt_temp = abs(sim.material.a * Laplace + dT_t);
		const double V_temp = dTdt_temp / G_temp;

		const double eqFrac_temp = 1 - exp((-4 * PI * sim.material.cet_N0 / 3) * pow(G_temp * (sim.material.cet_n + 1) / (pow(sim.material.cet_a * V_temp, 1 / sim.material.cet_n)), -3));

		set_G(G_temp, p);
		set_Gx(Gxu_temp, p);
		set_Gy(Gyu_temp, p);
		set_Gz(Gzu_temp, p);
		set_dTdt(dTdt_temp, p);
		set_V(V_temp, p);
		set_eqFrac(eqFrac_temp, p);
		add_RDF_cr(dTdt_temp, p);

		const double temp_x1 = dGxdx * Gxu_temp + dGydy * Gyu_temp + dGxdz * Gzu_temp;
		const double temp_y1 = dGxdy * Gxu_temp + dGydy * Gyu_temp + dGydz * Gzu_temp;
		const double temp_z1 = dGxdz * Gxu_temp + dGydz * Gyu_temp + dGzdz * Gzu_temp;

		const double temp_xyz1 = (temp_x1 * Gxu_temp + temp_y1 * Gyu_temp + temp_z1 * Gzu_temp);

		const double temp_x2 = Gxu_temp * temp_xyz1;
		const double temp_y2 = Gxu_temp * temp_xyz1;
		const double temp_z2 = Gxu_temp * temp_xyz1;
		
		const double Hx_temp = temp_x1 - temp_x2;
		const double Hy_temp = temp_y1 - temp_y2;
		const double Hz_temp = temp_z1 - temp_z2;

		const double H_temp = sqrt(Hx_temp * Hx_temp + Hy_temp * Hy_temp + Hz_temp * Hz_temp);
		const double Hxu_temp = Hx_temp / H_temp;
		const double Hyu_temp = Hy_temp / H_temp;
		const double Hzu_temp = Hz_temp / H_temp;

		set_H(H_temp, p);
		set_Hx(Hxu_temp, p);
		set_Hy(Hyu_temp, p);
		set_Hz(Hzu_temp, p);
	}
}
