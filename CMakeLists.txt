cmake_minimum_required(VERSION 3.22)

project(FMIndex-Project VERSION 0.2.0)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(FM-Index)
add_subdirectory(openbwt-v1.5)
add_subdirectory(UnitTests)
