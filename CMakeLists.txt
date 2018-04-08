cmake_minimum_required(VERSION 2.8)
project(aniso2)


set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(CMakeScripts_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${CMakeScripts_MODULE_DIR})

find_package(MATLAB REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DMATLAB_MEX_FILE -DRUN_OMP -DBBFMM_CACHE \
    -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread -DMX_COMPAT_32 \
    -fomit-frame-pointer -malign-double -fstrict-aliasing -march=native -std=gnu++11 -fopenmp -ffast-math")

include_directories(${MATLAB_INCLUDE_DIR} mexplus bbfmm)
set(MATLAB_LIBS ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY} ${MATLAB_MAT_LIBRARY} )
#${MATLAB_ROOT}/bin/glnxa64/libmwblas.so ${MATLAB_ROOT}/bin/glnxa64/libmwlapack.so

set(SOURCE_FILES
        Aniso.cpp Aniso.h
        Geometry.cpp Geometry.h
        KernelFactory.cpp KernelFactory.h
        Quadrature.cpp Quadrature.h
        bbfmm/bbfmm.h bbfmm/blas_wrapper.cpp
        bbfmm/linalg.h bbfmm/utils.h utility/Profiler.h
        utility/config.cpp utility/config.h
        gmres.h gmres.cpp matlab_io.h AnisoWrapper.cpp AnisoWrapper.h bbfmm/common.h)

#add_executable(aniso2 ${SOURCE_FILES})
add_library(AnisoWrapper SHARED ${SOURCE_FILES})

target_link_libraries(AnisoWrapper m blas quadmath ${MATLAB_LIBS})
set_target_properties(AnisoWrapper PROPERTIES SUFFIX ".mexa64")
set_target_properties(AnisoWrapper PROPERTIES PREFIX "")

set_target_properties(AnisoWrapper
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_ROOT}/bin"
        )