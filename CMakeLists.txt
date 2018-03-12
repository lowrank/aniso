cmake_minimum_required(VERSION 3.5)
project(aniso2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DRUN_OMP -DBBFMM_CACHE -march=native -std=c++11 -fopenmp ")

set(SOURCE_FILES main.cpp
        Aniso.cpp Aniso.h
        Geometry.cpp Geometry.h
        KernelFactory.cpp KernelFactory.h
        Quadrature.cpp Quadrature.h
        bbfmm/bbfmm.h bbfmm/blas_wrapper.cpp
        bbfmm/linalg.h bbfmm/utils.h utility/Profiler.h
        utility/config.cpp utility/config.h
        gmres.h gmres.cpp
        )
add_executable(aniso2 ${SOURCE_FILES})


target_link_libraries(aniso2 m blas quadmath)