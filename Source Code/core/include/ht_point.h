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

#ifndef POINT_H
#define POINT_H
#include <vector>
#include "data_structs.h"

namespace SAHTM
{
    class Point{
    public:
        Point();
        ~Point();
        void	Initialize(Simdat&);
        void	CalcGo_2(double, std::vector<PathSeg>&, Simdat&);
        double	Temp_Calc_Pre_Path(double t2, std::vector<int_seg>&, Simdat&, int, int);
        std::string Print();
        //Attributes
        double	x, y, z;							//Point location
        int		i, j, k;							//Point index
        //double	xb, yb, zb, qmod;				//Beam location and power multiplier
        double	T, Tlast, G, V, Gtemp, Gmag;		//Tempearture, previous temperature, thermal gradient, interface velocity, temporary G, global G
        double	Gx, Gy, Gz, Gx_temp, Gy_temp, Gz_temp, theta;
        double  laplace;
        double	Gxu, Gyu, Gzu;					//Components of unit vector in direction of solidification
        double	V_theta;							//dendrite tip velocity relative to underlying crystallographic direction
        double	TScale;							//Scale analysis variable
        double	kh, khmax;						//Varaibles for keyhole calculations
        bool	T_calc_flag;
        bool	s_flag;
        bool	check_flag;
        bool	output_flag;
        double  dTdt_sol; //Rate of cooling
        double	dT_cur; //Instantaneous heat source
    };
}
#endif // POINT_H
