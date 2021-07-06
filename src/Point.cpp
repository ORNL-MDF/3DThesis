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

#include <vector>
#include <cmath>

#include "DataStructs.h"
#include "Point.h"
#include "Calc.h"

#define PI 3.14159265358979323846

//Constructor and destructor
Point::Point(void){}
Point::~Point(void){}

// Set initial parameters
void Point::Initialize(Simdat& sim){
	//Uses very small value instead of 0 for 
	//variables that will be displayed in log scale
	Gtemp = Gmag = G = V = 1e-100;
	Gx = Gy = Gz = 0;
	theta = 0;
	eq_frac = -1;
	t_last_liq = t_last_sol = 0.0;
	T_calc_flag = s_flag = check_flag = output_flag = 0;
	T = Tlast = sim.mat.Tinit;
}

void Point::Solidify(double t2, vector<path_seg>& segv, Simdat& sim) {
	if (!sim.setting.out_mode) { return; }

	if (sim.setting.infBeams) {
		vector<infBeam> infBeams = sim.infBeams;
		t_last_liq = Point::Calc_Time(t2, infBeams, sim, sim.mat.T_liq);
		Point::Calc_Sol(infBeams, sim);
	}
	else {
		t_last_liq = Point::Calc_Time(t2, segv, sim, sim.mat.T_liq);
		Point::Calc_Sol(segv, sim);
	}

	Point::Calc_Sol_Dirs(sim);
	return;
}

double Point::Calc_Time(double t2, vector<path_seg>& segv, Simdat& sim, double T_cross) {

	int flag = 1;
	int iter = 0;
	int use_pt = 0;
	int max_iter = sim.setting.max_iter;
	double T2 = T;
	double T1 = Tlast;
	double T0 = T_cross;
	double T_temp = 0.0;
	double T_err = 0.0;
	double err_limit = sim.setting.dttest;
	double time2 = t2;
	double time1 = t2 - sim.param.dt;
	double time0 = 0;
	double m = 0.0;
	double dT = 0.0;

	vector<int_seg> isegv;

	while (flag) { //Continue until you get a temperature really close to the solidification temp
		m = (T2 - T1) / (time2 - time1);
		time0 = time1 + ((T0 - T1) / m);
		//Calc::GaussIntegrate(isegv, segv, sim, time0, 0);
		Calc::Integrate_thread(isegv, segv, sim, time0, 1);
		T_temp = Temp_Calc_Pre_Path(time0, isegv, sim, 0, 0);
		T_err = abs(1.0 - T_temp / T0);
		iter++;
		if ((T_err < err_limit) || (iter >= max_iter)) {
			flag = 0;
			dT = 0;

			//Calc Quadratic Derivative fitting 3 points
			/*double alpha, tau, a, b;
			alpha = 1.0 / (time2 - time1);
			tau = (time0 - time1) / (time2 - time1);
			a = (T0 - (1 - tau) * T1 - tau * T2) / ((1 - tau) * tau);
			b = (T0 - (1 - tau*tau) * T1 - tau*tau * T2) / ((1 - tau) * tau);
			dTdt_sol = (2 * a * tau + b) * alpha;*/

		}
		else{
			if (T_temp > T0) {
				time1 = time0;
				T1 = T_temp;
			}
			else {
				time2 = time0;
				T2 = T_temp;
			}
		}		
		isegv.clear();
	}

	return time0;
}

void Point::Calc_Sol(vector<path_seg>& segv, Simdat& sim) {
	vector<int_seg> isegv;
	double t0 = t_last_liq;
	Calc::Integrate_thread(isegv, segv, sim, t0, 1);
	Temp_Calc_Pre_Path(t0, isegv, sim, 0, 1);
}

void Point::Calc_Sol_Dirs(Simdat& sim) {
	if (Gtemp > (1e-50)) {
		G = Gtemp;

		double dTdt = sim.mat.a * laplace + dT_cur;

		dTdt_sol = -abs(dTdt);
		V = -dTdt_sol / G;

		//if (V < 0) { 
		//	V = -V; 
		//	//std::cout << "Velocity Error\n"; 
		//}

		Gx = Gx_temp;
		Gy = Gy_temp;
		Gz = Gz_temp;

		Gxu = Gx / G;
		Gyu = Gy / G;
		Gzu = Gz / G;

		if (sim.setting.out_mode == 3) {
			double tempx = dGxdx * Gxu + dGxdy * Gyu + dGxdz * Gzu;
			double tempy = dGxdy * Gxu + dGydy * Gyu + dGydz * Gzu;
			double tempz = dGxdz * Gxu + dGydz * Gyu + dGzdz * Gzu;

			double temp2x = Gxu * (tempx * Gxu + tempy * Gyu + tempz * Gzu);
			double temp2y = Gyu * (tempx * Gxu + tempy * Gyu + tempz * Gzu);
			double temp2z = Gzu * (tempx * Gxu + tempy * Gyu + tempz * Gzu);

			Hx = tempx - temp2x;
			Hy = tempy - temp2y;
			Hz = tempz - temp2z;

			double Hnorm = sqrt(Hx * Hx + Hy * Hy + Hz * Hz);
			Hxu = Hx / Hnorm;
			Hyu = Hy / Hnorm;
			Hzu = Hz / Hnorm;
		}

		//Calculation of equiaxed grain fraction according to the model by Gaumann et al., Acta Materialia, 2001.
		if (sim.mat.calc_CET && V > 0) { eq_frac = 1 - exp((-4 * PI * sim.mat.cet_N0 / 3) * pow(G * (sim.mat.cet_n + 1) / (pow(sim.mat.cet_a * V, 1 / sim.mat.cet_n)), -3)); }
	}
}

#pragma region InfBeamVerions
double Point::Calc_Time(double t2, vector<infBeam>& infBeams, Simdat& sim, double T_cross) {

	int flag = 1;
	bool switchSearch = false;
	int iter = 0;
	int max_iter = sim.setting.max_iter;
	double T2 = T;
	double T1 = Tlast;
	double T0 = T_cross;
	double T_temp = 0.0;
	double T_err = 0.0;
	double err_limit = sim.setting.dttest;
	double time2 = t2;
	double time1 = t2 - sim.param.dt;
	double time0 = 0;
	double m = 0.0;

	//Continue until you get a temperature really close to the solidification temp
	//False Position Method
	while (flag) { 
		m = (T2 - T1) / (time2 - time1);
		time0 = time1 + ((T0 - T1) / m);
		if (time0>time2 || time0<time1) {switchSearch = true; break;}
		Calc::GaussIntegrateInfBeams(infBeams, sim, time0, 0);
		T_temp = Temp_Calc_Pre_Path(time0, infBeams, sim, 0, 0);
		if (T_temp > T0) {
			time1 = time0;
			T1 = T_temp;
		}
		else {
			time2 = time0;
			T2 = T_temp;
		}
		T_err = fabs(1 - T_temp / T0);
		iter++;
		if ((T_err < err_limit) || (iter >= max_iter)) {flag = 0;}
	}

	//Continue until you get a temperature really close to the solidification temp
	//Binary Search
	iter = 0; time2 = t2; time1 = t2 - sim.param.dt;
	
	if (switchSearch) {
		while (flag) { 
			time0 = (time1 + time2) / 2.0;
			Calc::GaussIntegrateInfBeams(infBeams, sim, time0, 0);
			T_temp = Temp_Calc_Pre_Path(time0, infBeams, sim, 0, 0);
			if (T_temp > T0) {time1 = time0; T1 = T_temp;}
			else {time2 = time0;T2 = T_temp;}
			T_err = fabs(1 - T_temp / T0);
			iter++;
			if ((T_err < err_limit) || (iter >= max_iter)) {flag = 0;}
		}
	}

	return time0;
}

void Point::Calc_Sol(vector<infBeam>& infBeams, Simdat& sim) {
	double t0 = t_last_liq;
	Calc::GaussIntegrateInfBeams(infBeams, sim, t0, 1);
	Temp_Calc_Pre_Path(t0, infBeams, sim, 0, 1);
}
#pragma endregion InfBeamVerions

double Point::Temp_Calc_Pre_Path(double t2, vector<int_seg>& isegv, Simdat& sim, int set_T, int is_sol) {
	double dT;
	if (sim.setting.infBeams) { dT = Point::Temp_Calc_Pre_Path(t2,sim.infBeams,sim,set_T,is_sol); return dT ; }
	else {
		if (!is_sol || !sim.setting.out_mode) {dT = Point::Calc_T(t2, isegv, sim);}
		else if (sim.setting.out_mode != 3) { dT = Point::Calc_Solidification(t2, isegv, sim);}
		else {dT = Point::Calc_Secondary_Solidification(t2, isegv, sim); }

		if (!sim.setting.out_mode) { Gtemp = 1e-100; }
		else { Gtemp = sqrt(Gx_temp * Gx_temp + Gy_temp * Gy_temp + Gz_temp * Gz_temp); }

		if (set_T == 1) { T = dT + sim.mat.Tinit; T_calc_flag = 1; }
		else if (set_T == -1) { G = 1e-100; V = 1e-100; eq_frac = -1; output_flag = 1; Tlast = dT + sim.mat.Tinit; }

		if (T >= sim.mat.T_liq) {
			G = 1e-100;
			V = 1e-100;
			eq_frac = -1;
			output_flag = 1;
		}
		return dT + sim.mat.Tinit;
	}
}

double Point::Calc_T(double t2, vector<int_seg>& isegv, Simdat& sim) {

	/**************************************************************************************************
	Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
	Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
	**************************************************************************************************/

	double tau, dx, dy, dz, ct, beta, phix, phiy, phiz, phi, dTp, dT = 0;
	beta = pow(3 / PI, 1.5) * sim.beam.q / (sim.mat.rho * sim.mat.cps);
	int size = isegv.size();
	for (int iter = 0; iter < size; iter++) {
		tau = t2 - isegv[iter].taui;
		dx = x - isegv[iter].xb;
		dy = y - isegv[iter].yb;
		dz = z - isegv[iter].zb;

		ct = 12.0 * sim.mat.a * tau;
		phix = (sim.beam.ax * sim.beam.ax + ct);
		phiy = (sim.beam.ay * sim.beam.ay + ct);
		phiz = (sim.beam.az * sim.beam.az + ct);
		phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

		dTp = isegv[iter].dtau * isegv[iter].qmod * beta * phi;
		dT += dTp;
	}

	return dT;
}

double Point::Calc_Solidification(double t2, vector<int_seg>& isegv, Simdat& sim) {
	double tau, dx, dy, dz, ct, phix, phiy, phiz, beta, phi, dTp, dT = 0;

	Gx_temp = 0; Gy_temp = 0; Gz_temp = 0; laplace = 0; dT_cur = 0;

	/**************************************************************************************************
	Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
	Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
	**************************************************************************************************/

	beta = pow(3 / PI, 1.5) * sim.beam.q / (sim.mat.rho * sim.mat.cps);	//Beam power pre-factor
	int size = isegv.size();
	for (int iter = 0; iter < size; iter++) {
		tau = t2 - isegv[iter].taui;
		dx = x - isegv[iter].xb;
		dy = y - isegv[iter].yb;
		dz = z - isegv[iter].zb;

		ct = 12.0 * sim.mat.a * tau;
		phix = (sim.beam.ax * sim.beam.ax + ct);
		phiy = (sim.beam.ay * sim.beam.ay + ct);
		phiz = (sim.beam.az * sim.beam.az + ct);
		phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

		dTp = isegv[iter].dtau * isegv[iter].qmod * beta * phi;
		dT += dTp;

		double dpx = (-6.0 * dx / phix);
		double dpy = (-6.0 * dy / phiy);
		double dpz = (-6.0 * dz / phiz);
		double ddpx = (-6.0 / phix);
		double ddpy = (-6.0 / phiy);
		double ddpz = (-6.0 / phiz);

		Gx_temp += dTp * dpx;	//x-gradient
		Gy_temp += dTp * dpy;	//y-gradient
		Gz_temp += dTp * dpz;	//z-gradient

		laplace += dTp * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
		if (isegv[iter].dtau == 0) { dT_cur += (isegv[iter].qmod * beta * phi); }  //Notice dtau is not used, that is because we are looking the instantaneous change at that time
	}

	return dT;
}

double Point::Calc_Secondary_Solidification(double t2, vector<int_seg>& isegv, Simdat& sim) {
	double tau, dx, dy, dz, ct, phix, phiy, phiz, beta, phi, dTp, dT = 0;

	Gx_temp = 0; Gy_temp = 0; Gz_temp = 0; laplace = 0; dT_cur = 0;

	dGxdx = 0; dGxdy = 0; dGxdz = 0; dGydy = 0; dGydz = 0; dGzdz = 0;

	/**************************************************************************************************
	Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
	Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
	**************************************************************************************************/

	beta = pow(3 / PI, 1.5) * sim.beam.q / (sim.mat.rho * sim.mat.cps);	//Beam power pre-factor
	int size = isegv.size();
	for (int iter = 0; iter < size; iter++) {
		tau = t2 - isegv[iter].taui;
		dx = x - isegv[iter].xb;
		dy = y - isegv[iter].yb;
		dz = z - isegv[iter].zb;

		ct = 12.0 * sim.mat.a * tau;
		phix = (sim.beam.ax * sim.beam.ax + ct);
		phiy = (sim.beam.ay * sim.beam.ay + ct);
		phiz = (sim.beam.az * sim.beam.az + ct);
		phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

		dTp = isegv[iter].dtau * isegv[iter].qmod * beta * phi;
		dT += dTp;

		double dpx = (-6.0 * dx / phix);
		double dpy = (-6.0 * dy / phiy);
		double dpz = (-6.0 * dz / phiz);
		double ddpx = (-6.0 / phix);
		double ddpy = (-6.0 / phiy);
		double ddpz = (-6.0 / phiz);

		Gx_temp += dTp * dpx;	//x-gradient
		Gy_temp += dTp * dpy;	//y-gradient
		Gz_temp += dTp * dpz;	//z-gradient

		laplace += dTp * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
		if (isegv[iter].dtau == 0) { dT_cur += (isegv[iter].qmod * beta * phi); }  //Notice dtau is not used, that is because we are looking the instantaneous change at that time

		dGxdx += dTp * (dpx * dpx + ddpx);
		dGxdy += dTp * (dpx * dpy);
		dGxdz += dTp * (dpx * dpz);

		dGydy += dTp * (dpy * dpy + ddpy);
		dGydz += dTp * (dpy * dpz);

		dGzdz += dTp * (dpz * dpz + ddpz);
	}

	return dT;
}

#pragma region InfBeamVerions
double Point::Temp_Calc_Pre_Path(double t2, vector<infBeam>& infBeams, Simdat& sim, int set_T, int is_sol) {
	double dT;
	if (!is_sol || !sim.setting.out_mode) {
		dT = Point::Calc_T(t2, infBeams, sim);
	}
	else if (sim.setting.out_mode != 3) {
		dT = Point::Calc_Solidification(t2, infBeams, sim);
	}
	else {
		dT = Point::Calc_Secondary_Solidification(t2, infBeams, sim);
	}

	if (!sim.setting.out_mode) { Gtemp = 1e-100; }
	else { Gtemp = sqrt(Gx_temp * Gx_temp + Gy_temp * Gy_temp + Gz_temp * Gz_temp); }

	if (set_T == 1) { T = dT + sim.mat.Tinit; T_calc_flag = 1; }
	else if (set_T == -1) { G = 1e-100; V = 1e-100; eq_frac = -1; output_flag = 1; Tlast = dT + sim.mat.Tinit; }

	if (T >= sim.mat.T_liq) {
		G = 1e-100;
		V = 1e-100;
		eq_frac = -1;
		output_flag = 1;
	}

	return dT + sim.mat.Tinit;
}

double Point::Calc_T(double t2, vector<infBeam>& infBeams,  Simdat& sim) {
	double tau, dx, dy, dz, ct, beta, phix, phiy, phiz, phi, dTp, dT = 0;
	for (infBeam beam : infBeams) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/
		beta = pow(3 / PI, 1.5) * beam.q / (sim.mat.rho * sim.mat.cps);
		int size = beam.issegv.size(); 
		for (int iter = 0; iter < size; iter++) {
			tau = t2 - beam.issegv[iter].taui;
			dx = x - beam.issegv[iter].xb;
			dy = y - beam.issegv[iter].yb;
			dz = z - beam.issegv[iter].zb;

			ct = 12.0 * sim.mat.a * tau;
			phix = (beam.issegv[iter].ax * beam.issegv[iter].ax + ct);
			phiy = (beam.issegv[iter].ay * beam.issegv[iter].ay + ct);
			phiz = (beam.issegv[iter].az * beam.issegv[iter].az + ct);
			phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

			dTp = beam.issegv[iter].dtau * beam.issegv[iter].qmod * beta * phi;
			dT += dTp;
		}
	}
	return dT;
}

double Point::Calc_Solidification(double t2, vector<infBeam>& infBeams, Simdat& sim) {
	double tau, dx, dy, dz, ct, phix, phiy, phiz, beta, phi, dTp, dT = 0;
	Gx_temp = 0; Gy_temp = 0; Gz_temp = 0; laplace = 0; dT_cur = 0;
	for (infBeam beam : infBeams) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/
		beta = pow(3 / PI, 1.5) * beam.q / (sim.mat.rho * sim.mat.cps);	//Beam power pre-factor
		int size = beam.issegv.size();
		for (int iter = 0; iter < size; iter++) {
			tau = t2 - beam.issegv[iter].taui;
			dx = x - beam.issegv[iter].xb;
			dy = y - beam.issegv[iter].yb;
			dz = z - beam.issegv[iter].zb;

			ct = 12.0 * sim.mat.a * tau;
			phix = (beam.issegv[iter].ax * beam.issegv[iter].ax + ct);
			phiy = (beam.issegv[iter].ay * beam.issegv[iter].ay + ct);
			phiz = (beam.issegv[iter].az * beam.issegv[iter].az + ct);
			phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

			dTp = beam.issegv[iter].dtau * beam.issegv[iter].qmod * beta * phi;
			dT += dTp;

			double dpx = (-6.0 * dx / phix);
			double dpy = (-6.0 * dy / phiy);
			double dpz = (-6.0 * dz / phiz);
			double ddpx = (-6.0 / phix);
			double ddpy = (-6.0 / phiy);
			double ddpz = (-6.0 / phiz);

			Gx_temp += dTp * dpx;	//x-gradient
			Gy_temp += dTp * dpy;	//y-gradient
			Gz_temp += dTp * dpz;	//z-gradient

			laplace += dTp * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
			if (beam.issegv[iter].dtau == 0) { dT_cur += (beam.issegv[iter].qmod * beta * phi); }  //Notice dtau is not used, that is because we are looking the instantaneous change at that time
		}
	}
	return dT;
}

double Point::Calc_Secondary_Solidification(double t2, vector<infBeam>& infBeams, Simdat& sim) {
	double tau, dx, dy, dz, ct, phix, phiy, phiz, beta, phi, dTp, dT = 0;

	Gx_temp = 0; Gy_temp = 0; Gz_temp = 0; laplace = 0; dT_cur = 0;

	dGxdx = 0; dGxdy = 0; dGxdz = 0; dGydy = 0; dGydz = 0; dGzdz = 0;
	for (infBeam beam : infBeams) {
		/**************************************************************************************************
		Heat transfer kernel for ellipsoidal Gaussian volumentric heat source
		Adapted from Nguyen et al., Welding Journal, 1999 (Eq. 7)
		**************************************************************************************************/
		beta = pow(3 / PI, 1.5) * beam.q / (sim.mat.rho * sim.mat.cps);	//Beam power pre-factor
		int size = beam.issegv.size();
		for (int iter = 0; iter < size; iter++) {
			tau = t2 - beam.issegv[iter].taui;
			dx = x - beam.issegv[iter].xb;
			dy = y - beam.issegv[iter].yb;
			dz = z - beam.issegv[iter].zb;

			ct = 12.0 * sim.mat.a * tau;
			phix = (beam.issegv[iter].ax * beam.issegv[iter].ax + ct);
			phiy = (beam.issegv[iter].ay * beam.issegv[iter].ay + ct);
			phiz = (beam.issegv[iter].az * beam.issegv[iter].az + ct);
			phi = exp(-3.0 * ((dx * dx / phix) + (dy * dy / phiy) + (dz * dz / phiz))) / (sqrt(phix * phiy * phiz));

			dTp = beam.issegv[iter].dtau * beam.issegv[iter].qmod * beta * phi;
			dT += dTp;

			double dpx = (-6.0 * dx / phix);
			double dpy = (-6.0 * dy / phiy);
			double dpz = (-6.0 * dz / phiz);
			double ddpx = (-6.0 / phix);
			double ddpy = (-6.0 / phiy);
			double ddpz = (-6.0 / phiz);

			Gx_temp += dTp * dpx;	//x-gradient
			Gy_temp += dTp * dpy;	//y-gradient
			Gz_temp += dTp * dpz;	//z-gradient

			laplace += dTp * (dpx * dpx + dpy * dpy + dpz * dpz + ddpx + ddpy + ddpz); //laplacian
			if (beam.issegv[iter].dtau == 0) { dT_cur += (beam.issegv[iter].qmod * beta * phi); }  //Notice dtau is not used, that is because we are looking the instantaneous change at that time

			dGxdx += dTp * (dpx * dpx + ddpx);
			dGxdy += dTp * (dpx * dpy);
			dGxdz += dTp * (dpx * dpz);

			dGydy += dTp * (dpy * dpy + ddpy);
			dGydz += dTp * (dpy * dpz);

			dGzdz += dTp * (dpz * dpz + ddpz);
		}
	}
	return dT;
}
#pragma endregion InfBeamVerions

void Point::find_t_last_liq(double t2, vector<int_seg>& isegv, Simdat& sim) {
	double T_temp = Temp_Calc_Pre_Path(t2, isegv, sim, 0, 0);
	if (T_temp > sim.mat.T_liq) {t_last_liq = t2;}
}

double Point::find_t_last_heat(double t, vector<path_seg>& segv, Simdat& sim) {

	vector<int_seg> isegv;
	Calc::GaussIntegrate(isegv, segv, sim, t, 0);
	double T2 = Temp_Calc_Pre_Path(t, isegv, sim, 0, 0);
	double T1;
	
	while (true) { //Continue until you are liquid
		isegv.clear();
		t -= sim.param.dt;
		Calc::GaussIntegrate(isegv, segv, sim, t, 0);
		T1 = Temp_Calc_Pre_Path(t, isegv, sim, 0, 0);
		if (T1 < T2) { break; }
	}
	return t;
}