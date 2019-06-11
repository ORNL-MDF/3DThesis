cmake_minimum_required(VERSION 3.6.0)

project(SAHTM_Lib)

option(USE_OpenMP "Enable OpenMP for multithreaded speed ups" ON)
if(USE_OpenMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

file (GLOB_RECURSE SAHTM_SOURCES
      "src/*.cpp"
	  "include/*h"
	  )
	  
include_directories(include)

add_library (SAHTM_Lib ${SAHTM_SOURCES})

target_include_directories(SAHTM_Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
