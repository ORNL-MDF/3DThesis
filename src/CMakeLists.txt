file(GLOB THESIS_HEADERS GLOB *.h)
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

add_executable(3DThesis Main.cpp)
target_link_libraries(3DThesis LINK_PUBLIC Thesis)
install(TARGETS 3DThesis DESTINATION ${CMAKE_INSTALL_BINDIR})
