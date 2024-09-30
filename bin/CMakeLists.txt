add_executable(3DThesis Main.cpp)
target_link_libraries(3DThesis LINK_PUBLIC Thesis ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
install(TARGETS 3DThesis DESTINATION ${CMAKE_INSTALL_BINDIR})
