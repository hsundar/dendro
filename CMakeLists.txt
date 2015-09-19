cmake_minimum_required(VERSION 2.8)
project(dendro)

# version number 
set (Dendro_VERSION_MAJOR 4)
set (Dendro_VERSION_MINOR 0)


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
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
                    $ENV{PETSC_DIR}/include 
                    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include 
                    )

link_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib)

# Sources
file(GLOB OMG_SOURCES src/omg/*.C)
file(GLOB ODA_SOURCES src/oda/*.C)
file(GLOB ODA_PRIVATE_SOURCES src/oda/private/*.C)
file(GLOB OCT_SOURCES src/oct/*.C)
file(GLOB HILBERT_SOURCES src/hilbert/*.cpp)


#### Libraries ...
set(DENDRO_LIBS OMG ODA Test Oct Par Point BinOps Sys PC Hilbert Viz)
add_library(OMG ${OMG_SOURCES})
add_library(ODA ${ODA_SOURCES} ${ODA_PRIVATE_SOURCES})
add_library(Viz src/visualization/treenode2vtk.C)
add_library(Oct ${OCT_SOURCES})
add_library(Par src/par/parUtils.C)
add_library(BinOps src/binOps/binUtils.C)
add_library(Hilbert ${HILBERT_SOURCES})
add_library(Point src/point/Point.C)
add_library(Test src/test/testUtils.C)
add_library(Sys src/sys/sys.C)
add_library(PC src/pc/blockDiag.C)


## executables 

##---------------------------------------------------------------------------------------
##  Helper APPS - No dendro
##---------------------------------------------------------------------------------------
add_executable(genGauss scripts/genGaussSeq.C)
add_executable(splitPoints scripts/splitPoints.C)
##---------------------------------------------------------------------------------------



##---------------------------------------------------------------------------------------
## 2:1-Balance Apps 
##---------------------------------------------------------------------------------------
add_executable(testRipple examples/src/drivers/tstRipple.C)
target_link_libraries(testRipple ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(runScal examples/src/drivers/runScal.C)
target_link_libraries(runScal ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(rippleBal examples/src/drivers/rippleBal.C)
target_link_libraries(rippleBal ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(seqRipple examples/src/drivers/seqRipple.C)
target_link_libraries(seqRipple ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(seqRipplePtr examples/src/drivers/seqRipplePtr.C)
target_link_libraries(seqRipplePtr ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(justBal examples/src/drivers/justBal.C)
target_link_libraries(justBal ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(testConAndBal examples/src/drivers/testConAndBal)
target_link_libraries(testConAndBal ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(testConBalAndCoarsen
  examples/src/drivers/testConBalAndCoarsen)
target_link_libraries(testConBalAndCoarsen ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)

add_executable(testPartition examples/src/drivers/testPartition.C)
target_link_libraries(testPartition ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)
##---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
##    OTHER APPS
##---------------------------------------------------------------------------------------
#add_executable(justBal examples/src/drivers/)
#target_link_libraries( ${DENDRO_LIBS} petsc ${MPI_LIBRARIES} m)
