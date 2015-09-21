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
option(USE_NCA_PROPERTY "Use the NCA based approach" ON)
option(BUILD_TESTS "Build test applications" ON)
option(BUILD_DA_EXAMPLES "Build examples using dendro::DA" ON)
option(BUILD_MG_EXAMPLES "Build test applications using dendro::MG" ON)

##------
include_directories(${MPI_INCLUDE_PATH}
                    include
                    include/fem
                    include/oda 
                    include/omg 
                    include/test 
                    examples/include 
                    $ENV{PETSC_DIR}/include 
                    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include 
                    )

link_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib)

# Sources
file(GLOB OMG_SOURCES src/omg/*.cpp)
file(GLOB ODA_SOURCES src/oda/*.cpp)
file(GLOB ODA_PRIVATE_SOURCES src/oda/private/*.cpp)

file(GLOB DENDRO_SOURCES src/*.cpp)
file(GLOB DENDRO_HEADERS include/*.h)
file(GLOB DENDRO_TEMPLATES include/*.tcc)

#### Libraries ...
set(DENDRO_LIBS dendroMG dendroDA dendroTest dendro)

add_library(dendroMG ${OMG_SOURCES})
add_library(dendroDA ${ODA_SOURCES} ${ODA_PRIVATE_SOURCES})
add_library(dendroTest src/test/testUtils.C)

add_library(dendro ${DENDRO_SOURCES} 
                   ${DENDRO_HEADERS}
                   ${DENDRO_TEMPLATES})

## executables 

##---------------------------------------------------------------------------------------
##  Helper APPS
##---------------------------------------------------------------------------------------
add_executable(genGauss scripts/genGaussSeq.C)
add_executable(splitPoints scripts/splitPoints.C)

add_executable(gen2spheres scripts/gen2spheres.C)
target_link_libraries(gen2spheres dendro ${MPI_LIBRARIES} m)

add_executable(genPrlFBMpts scripts/genPrlFBMpts.C)
target_link_libraries(genPrlFBMpts ${MPI_LIBRARIES} m)
##---------------------------------------------------------------------------------------


##---------------------------------------------------------------------------------------
## Octree Construct-Balance Apps
##---------------------------------------------------------------------------------------
add_executable(testRipple examples/src/drivers/tstRipple.C)
target_link_libraries(testRipple dendro petsc ${MPI_LIBRARIES} m)

add_executable(runScal examples/src/drivers/runScal.C)
target_link_libraries(runScal dendro petsc ${MPI_LIBRARIES} m)

add_executable(rippleBal examples/src/drivers/rippleBal.C)
target_link_libraries(rippleBal dendro petsc ${MPI_LIBRARIES} m)

add_executable(seqRipple examples/src/drivers/seqRipple.C)
target_link_libraries(seqRipple dendro petsc ${MPI_LIBRARIES} m)

add_executable(seqRipplePtr examples/src/drivers/seqRipplePtr.C)
target_link_libraries(seqRipplePtr dendro petsc ${MPI_LIBRARIES} m)

add_executable(justBal examples/src/drivers/justBal.C)
target_link_libraries(justBal dendro petsc ${MPI_LIBRARIES} m)

add_executable(testConAndBal examples/src/drivers/testConAndBal)
target_link_libraries(testConAndBal dendro petsc ${MPI_LIBRARIES} m)

add_executable(testConBalAndCoarsen
  examples/src/drivers/testConBalAndCoarsen)
target_link_libraries(testConBalAndCoarsen dendro petsc ${MPI_LIBRARIES} m)
##---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
##    Dendro DA APPS
##---------------------------------------------------------------------------------------
if(BUILD_DA_EXAMPLES)
    add_executable(testPartition examples/src/drivers/testPartition.C)
    target_link_libraries(testPartition dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(tstMatVec examples/src/drivers/tstMatVec.C examples/src/backend/odaJac.C examples/src/backend/handleType2Stencils.C)
    target_link_libraries(tstMatVec dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(buildDA examples/src/drivers/buildDA.C)
    target_link_libraries(buildDA dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(buildRgDA examples/src/drivers/buildRgDA.C)
    target_link_libraries(buildRgDA dendroDA dendro petsc ${MPI_LIBRARIES} m)

    #add_executable(octLaplacian examples/src/drivers/octLaplacian.C)
    #target_link_libraries(octLaplacian dendroDA dendro petsc ${MPI_LIBRARIES} m)
endif()
##---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
##    Dendro MG APPS
##---------------------------------------------------------------------------------------
if(BUILD_MG_EXAMPLES)
add_executable(checkError   examples/src/drivers/checkError.C
                            examples/src/backend/omgJac.C
                            examples/src/backend/matProp.C
                            examples/src/backend/odaJac.C
                            examples/src/backend/omgRhs.C
                            examples/src/backend/handleType2Stencils.C
                            examples/src/backend/handleType1Stencils.C)
target_link_libraries(checkError dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(newElasMesh  examples/src/drivers/newElasMesh.C
                            examples/src/backend/elasticityJac.C
                            examples/src/backend/elasticityRhs.C
                            examples/src/backend/omgJac.C
                            examples/src/backend/odaJac.C
                            examples/src/backend/omgRhs.C
                            examples/src/backend/vecMass.C
                            examples/src/backend/handleType2Stencils.C )
target_link_libraries(newElasMesh dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(newElasSolver    examples/src/drivers/newElasSolver.C
                                examples/src/backend/elasticityJac.C
                                examples/src/backend/elasticityRhs.C
                                examples/src/backend/omgJac.C
                                examples/src/backend/odaJac.C
                                examples/src/backend/omgRhs.C
                                examples/src/backend/vecMass.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(newElasSolver dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(elasticitySolver examples/src/drivers/elasticitySolver.C
                                examples/src/backend/elasticityJac.C
                                examples/src/backend/elasticityRhs.C
                                examples/src/backend/omgJac.C
                                examples/src/backend/odaJac.C
                                examples/src/backend/omgRhs.C
                                examples/src/backend/vecMass.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(elasticitySolver dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(tstMg    examples/src/drivers/tstMg.C
                        examples/src/backend/omgJac.C
                        examples/src/backend/matProp.C
                        examples/src/backend/odaJac.C
                        examples/src/backend/omgRhs.C
                        examples/src/backend/handleType2Stencils.C
                        examples/src/backend/handleType1Stencils.C )
target_link_libraries(tstMg dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(omgNeumann_ex1   examples/src/drivers/omgNeumann_ex1.C
                                examples/src/backend/omgNeumann.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(omgNeumann_ex1 dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(omgNeumann_ex2   examples/src/drivers/omgNeumann_ex2.C
                                examples/src/backend/omgNeumann.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(omgNeumann_ex2 dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(omgNeumann_ex3   examples/src/drivers/omgNeumann_ex3.C
                                examples/src/backend/omgNeumann.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(omgNeumann_ex3 dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(omgNeumann_2spheres  examples/src/drivers/omgNeumann_2spheres.C
                                    examples/src/backend/omgNeumann.C
                                    examples/src/backend/handleType2Stencils.C )
target_link_libraries(omgNeumann_2spheres dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(tstMgReg examples/src/drivers/tstRegularMg.C
                        examples/src/backend/omgJac.C
                        examples/src/backend/matProp.C
                        examples/src/backend/odaJac.C
                        examples/src/backend/omgRhs.C
                        examples/src/backend/handleType2Stencils.C
                        examples/src/backend/handleType1Stencils.C )
target_link_libraries(tstMgReg dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(testElasMatVec   examples/src/drivers/testElasMatVec.C
                                examples/src/backend/elasticityJac.C
                                examples/src/backend/elasticityRhs.C
                                examples/src/backend/omgJac.C
                                examples/src/backend/odaJac.C
                                examples/src/backend/omgRhs.C
                                examples/src/backend/vecMass.C
                                examples/src/backend/handleType2Stencils.C )
target_link_libraries(testElasMatVec dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(pts2Mesh examples/src/drivers/pts2Mesh.C
                        examples/src/backend/elasticityJac.C
                        examples/src/backend/elasticityRhs.C
                        examples/src/backend/omgJac.C
                        examples/src/backend/odaJac.C
                        examples/src/backend/omgRhs.C
                        examples/src/backend/vecMass.C
                        examples/src/backend/handleType2Stencils.C )
target_link_libraries(pts2Mesh dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

add_executable(createCoarser    examples/src/drivers/createCoarser.C)
target_link_libraries(createCoarser dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)
endif()
##---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
##    Dendro Test APPS
##---------------------------------------------------------------------------------------
if(BUILD_TESTS)
    add_executable(checkNl  src/test/compareNlFiles.C)
    target_link_libraries(checkNl dendro petsc ${MPI_LIBRARIES} m)

    add_executable(checkMl  src/test/compareMFiles.C)
    target_link_libraries(checkMl dendro petsc ${MPI_LIBRARIES} m)

    add_executable(dumpMesh examples/src/drivers/dumpMesh.C)
    target_link_libraries(dumpMesh dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(downSample examples/src/drivers/downSample.C)
    target_link_libraries(downSample dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(checkResults src/test/checkResults.C)
    target_link_libraries(checkResults dendroTest dendro petsc ${MPI_LIBRARIES} m)

    add_executable(checkScatterValues src/test/checkScatterValues.C)
    target_link_libraries(checkScatterValues dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(checkCoarsening src/test/checkCoarsening.C)
    target_link_libraries(checkCoarsening dendroTest dendro petsc ${MPI_LIBRARIES} m)

    add_executable(justCoarsen examples/src/drivers/justCoarsen.C)
    target_link_libraries(justCoarsen dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(testNodeAndVals src/test/testNodeAndVals.C)
    target_link_libraries(testNodeAndVals dendroTest dendro petsc ${MPI_LIBRARIES} m)

    add_executable(testPetscInt src/test/testPetscInt.C)
    target_link_libraries(testPetscInt dendroTest dendro petsc ${MPI_LIBRARIES} m)


    add_executable(tstFBM   examples/src/drivers/tstFBM.C
                            examples/src/backend/omgJac.C
                            examples/src/backend/matProp.C
                            examples/src/backend/odaJac.C
                            examples/src/backend/omgRhs.C
                            examples/src/backend/handleType2Stencils.C
                            examples/src/backend/handleType1Stencils.C )
    target_link_libraries(tstFBM dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(tstFBM2  examples/src/drivers/tstFBM2.C
                            examples/src/backend/omgJac.C
                            examples/src/backend/matProp.C
                            examples/src/backend/odaJac.C
                            examples/src/backend/omgRhs.C
                            examples/src/backend/handleType2Stencils.C
                            examples/src/backend/handleType1Stencils.C )
    target_link_libraries(tstFBM2 dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

    add_executable(tstDirichlet examples/src/drivers/tstDirichlet.C
                                examples/src/backend/omgJac.C
                                examples/src/backend/matProp.C
                                examples/src/backend/odaJac.C
                                examples/src/backend/omgRhs.C
                                examples/src/backend/handleType2Stencils.C
                                examples/src/backend/handleType1Stencils.C )
    target_link_libraries(tstDirichlet dendroMG dendroDA dendro petsc ${MPI_LIBRARIES} m)

endif()
