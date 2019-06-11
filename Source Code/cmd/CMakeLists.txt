cmake_minimum_required(VERSION 3.6.0)

project(SAHTM_CMD)
if(USE_OpenMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()
add_executable (SAHTM_CMD main.cpp)

target_link_libraries(SAHTM_CMD LINK_PUBLIC SAHTM_Lib)
