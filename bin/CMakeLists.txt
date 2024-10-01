add_executable(3DThesis Main.cpp)
target_link_libraries(3DThesis LINK_PUBLIC Thesis MPI::MPI_CXX)
install(TARGETS 3DThesis DESTINATION ${CMAKE_INSTALL_BINDIR})
