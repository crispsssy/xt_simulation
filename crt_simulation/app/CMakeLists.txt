#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(GainSaturation)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)
find_package(ROOT REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(GainSaturation GainSaturation.cxx)
target_link_libraries(GainSaturation Garfield)
target_link_libraries(GainSaturation ${ROOT_LIBRARIES})

#---Copy all files locally to the build directory-------------------------------

