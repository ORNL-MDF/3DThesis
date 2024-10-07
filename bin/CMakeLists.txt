add_executable(3DThesis Main.cpp)
target_link_libraries(3DThesis LINK_PUBLIC Thesis)

if(Thesis_ENABLE_MPI)
  target_link_libraries(3DThesis LINK_PUBLIC MPI::MPI_CXX)
endif()

install(TARGETS 3DThesis DESTINATION ${CMAKE_INSTALL_BINDIR})
