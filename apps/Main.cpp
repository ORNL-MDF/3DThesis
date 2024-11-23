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

#include <Thesis_Core.hpp>

#ifdef Thesis_ENABLE_MPI
int main(int argc, char * argv[]) {	
	// Initialize MPI
    MPI_Init(&argc, &argv);
	Thesis::Run::Classic(argc, argv);
	// Finalize MPI
    MPI_Finalize();

	return 0;
}
#else
int main(int argc, char * argv[]) {	
	Thesis::Run::Classic(argc, argv);
	return 0;
}
#endif
