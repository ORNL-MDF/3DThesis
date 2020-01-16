//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

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

#pragma once
#include <cmath>
class Point{
private:
	//Attributes
	double	x, y, z;							//Point location
	int		i, j, k;							//Point index
	double	T, Tlast, G, V, Gtemp, Gmag;		//Tempearture, previous temperature, thermal gradient, interface velocity, temporary G, global G
	double	Gx, Gy, Gz, Gx_temp, Gy_temp, Gz_temp, theta;
	double  laplace;
	double	Gxu, Gyu, Gzu;					//Components of unit vector in direction of solidification
	double  eq_frac;
	int		T_calc_flag;
	int		s_flag;
	int		check_flag;
	int		output_flag;
	double  dTdt_sol; //Rate of cooling
	double	dT_cur; //Instantaneous heat source
	double	t_last_liq, t_last_sol; //Last time the point was liquid

	double dGxdx, dGxdy, dGxdz;
	double dGydy, dGydz;
	double dGzdz;

	double Hx, Hy, Hz;
	double Hxu, Hyu, Hzu;

public:
	//Functions
	Point();
	~Point();

	//////////////////////////////////////////////////////////
	///////////// FUNCTIONS IN POINT.CPP FILE ////////////////
	//////////////////////////////////////////////////////////
	void	Initialize();																//Initialize values
	double	Temp_Calc_Pre_Path(double, std::vector<int_seg>&, Simdat&, int, int);		//Calculate the temperature of a point at some time
		double  Calc_T(double, std::vector<int_seg>&, Simdat&);						    //Sub_mode for Temp_Calc_Pre_Path
		double  Calc_Solidification(double, std::vector<int_seg>&, Simdat&);
		double  Calc_Secondary_Solidification(double, std::vector<int_seg>&, Simdat&);
	double 	Calc_Time(double, std::vector<path_seg>&, Simdat&, double);					//Find when a point crosses a specific temperature (must be getting tracked...so mode 2 must be >sim.t_liq
	void	CalcGo_2(double, std::vector<path_seg>&, Simdat&);							//Find when a point solidifies and calculate the solidification conditions of the point
	void	find_t_last_liq(double, std::vector<int_seg>&, Simdat&);					//Find the last time a point was liquid
	double	find_t_last_heat(double, std::vector<path_seg>&, Simdat&);					//Find the last time a point was heating up
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////


	void	set_xloc(double xp) { x = xp; }	
	void	set_yloc(double yp) { y = yp; }
	void	set_zloc(double zp) { z = zp; }
	void	set_i(int ip) { i = ip; }
	void	set_j(int jp) { j = jp; }
	void	set_k(int kp) { k = kp; }

	void	set_s_flag(int f) { s_flag = f; }
	void	set_T_calc_flag(){ T_calc_flag = 1; }
	void	set_check_flag(){ check_flag = 1; }

	void	set_T(double Ttemp) { T = Ttemp; }
	void	set_Tlast() { Tlast = T; }
	void	set_Tlast_val(double T) { Tlast = T; }
	void	set_G(double Gin){ G = Gin; }
	void	set_V(double Vin){ V = Vin; }
	void	set_output_flag(int out){ output_flag = out; }
	
	void	init_check_flag(){ check_flag = 0; }
	void	init_T_calc_flag(){ T_calc_flag = 0; }

	double	get_x(){ return x; }
	double	get_y(){ return y; }
	double	get_z(){ return z; }

	int		get_i(){ return i; }
	int		get_j(){ return j; }
	int		get_k(){ return k; }

	double	get_T(){ return T; }
	double	get_Tlast(){ return Tlast; }
	double	get_G(){ return G; }
	double	get_V(){ return V; }
	double	get_Gmag(){ return Gmag; }
	double	get_Gx() { return Gx; }
	double	get_Gy() { return Gy; }
	double	get_Gz() { return Gz; }
	double	get_theta() { return theta; }

	double	get_s_flag() { return s_flag; }
	double	get_Gxu() { return Gxu; }
	double	get_Gyu() { return Gyu; }
	double	get_Gzu() { return Gzu; }
	double  get_eq_frac() { return eq_frac; }

	void	set_eq_frac(Simdat& sim) {
		if (V > 0) { 
			eq_frac = 1 - exp((-4 * 3.1415 * sim.mat.cet_N0 / 3) * pow(G * (sim.mat.cet_n + 1) / (pow(sim.mat.cet_a * V, 1 / sim.mat.cet_n)), -3)); 
		}
	}

	int		get_T_calc_flag() { return T_calc_flag; }
	int		get_check_flag() { return check_flag; }
	int		get_output_flag() { return output_flag; }

	double  get_dTdt_sol() { return dTdt_sol; }
	void	set_dTdt_sol(double dTdttemp) { dTdt_sol = dTdttemp; }

	void	set_t_last_liq(double t) { t_last_liq = t; }
	void	set_t_last_sol(double t) { t_last_sol = t; }
	double  get_t_last_liq() { return t_last_liq; }
	double  get_t_last_sol() { return t_last_sol; }

	double get_H()	{return sqrt(Hx*Hx + Hy * Hy + Hz * Hz);}

	double get_Hx() { return Hx; }
	double get_Hy() { return Hy; }
	double get_Hz() { return Hz; }

	double get_Hxu() { return Hxu; }
	double get_Hyu() { return Hyu; }
	double get_Hzu() { return Hzu; }
};