cmake_minimum_required(VERSION 2.8)
project(dendro)

# version number 
set (Dendro_VERSION_MAJOR 4)
set (Dendro_VERSION_MINOR 0)

# so that it canfind PETSc
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(PETSc COMPONENTS CXX)

if(OPENMP_FOUND)
  set ($CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set ($CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# Options 
option(__USE_PVT_DA_IN_MG__ "Use private DA in Multigrid" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)

##------
# initially a copy of makefileCore ...

include_directories(${MPI_INCLUDE_PATH}
                    include
                    include/oct
                    include/fem
                    include/oda 
                    include/omg 
                    include/seq 
                    include/par 
                    include/omp_par
                    include/point 
                    include/test 
                    include/binOps 
                    include/pc 
                    include/sys 
                    examples/include 
                    include/hilbert 
                    include/visualization 
                    ${PETSC_INCLUDES}
                    )

# Sources
file(GLOB OMG_SOURCES src/omg/*.C)
file(GLOB ODA_SOURCES src/oda/*.C)
file(GLOB ODA_PRIVATE_SOURCES src/oda/private/*.C)
file(GLOB OCT_SOURCES src/oct/*.C)
file(GLOB HILBERT_SOURCES src/hilbert/*.cpp)


#### Libraries ...
add_library(sys src/sys/sys.C)
add_library(OMG ${OMG_SOURCES})
add_library(ODA ${ODA_SOURCES} ${ODA_PRIVATE_SOURCES})
add_library(Viz src/visualization/treenode2vtk.C)
add_library(Oct ${OCT_SOURCES})
add_library(Par src/par/parUtils.C)
add_library(BinOps src/binOps/binUtils.C)
add_library(Hilbert ${HILBERT_SOURCES})
add_library(Point src/point/Point.C)
add_library(Test src/test/testUtils.C)

##------

