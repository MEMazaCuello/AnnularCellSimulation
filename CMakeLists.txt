cmake_minimum_required(VERSION 3.15)

project(AnnularCell VERSION 2.2)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_DIR "src")
set(SOURCES
    src/main.cpp
    src/analysis.hpp
    src/analysis.cpp
    src/annularCell.hpp
    src/annularCell.cpp
    src/grid.hpp
    src/grid.cpp
    src/rod.hpp
    src/rod.cpp
    src/GlobalParameters.hpp
)
add_executable(AnnularCell ${SOURCES})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT AnnularCell)