configure_file(ThesisConfig.h.cmakein ThesisConfig.h)

set(THESIS_HEADERS
  Calc.h
  DataStructs.h
  Grid.h
  Init.h
  Melt.h
  Out.h
  Run.h
  Util.h
)
if(Thesis_ENABLE_MPI)
  list(APPEND THESIS_HEADERS MpiStructs.h)
endif()

file(GLOB THESIS_SOURCES GLOB *.cpp)

install(FILES ${THESIS_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(Thesis STATIC ${THESIS_SOURCES} ${THESIS_HEADERS})

target_include_directories(Thesis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS Thesis
  EXPORT Thesis_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(Thesis::Thesis ALIAS Thesis)
target_compile_features(Thesis PUBLIC cxx_std_11)
target_link_libraries(Thesis PUBLIC OpenMP::OpenMP_CXX)

if(Thesis_ENABLE_MPI)
  target_link_libraries(Thesis PUBLIC MPI::MPI_CXX)
endif()