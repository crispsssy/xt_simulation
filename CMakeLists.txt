cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(xt_sim)
set(CMAKE_CXX_STANDARD 14)

find_package(Garfield REQUIRED)
find_package(ROOT REQUIRED)

add_library(XT ${CMAKE_CURRENT_SOURCE_DIR}/src/xt_simulation.cxx)
target_include_directories(XT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
         ${ROOT_INCLUDE_DIRS})
target_link_libraries(XT Garfield::Garfield)
target_link_libraries(XT ${ROOT_LIBRARIES})

add_executable(xt_sim ${CMAKE_CURRENT_SOURCE_DIR}/src/xt_simulation.cxx)
target_link_libraries(xt_sim PUBLIC XT)
#target_link_libraries(xt_sim Garfield)
#target_link_libraries(xt_sim ${ROOT_LIBRARIES})

install(TARGETS xt_sim
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
