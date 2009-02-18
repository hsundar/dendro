
SHELL = /bin/sh
include ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf
include ${PETSC_DIR}/bmake/common/variables
CEXT = C
CFLAGS = -D__USE_PVT_DA_IN_MG__ -DPETSC_USE_LOG

#-D__SILENT_MODE__
#-D__USE_64_BIT_INT__
#-D__USE_MG_INIT_TYPE2__ 
#-D__USE_MG_INIT_TYPE3__
#-D__USE_A2A_FOR_MPI_ALLGATHER__ 
#-D__BLOCK_PART_EQUALS_MORTON_PART__
#-D__NO_GHOST_LOOP__
#-D__PROFILE_WITH_BARRIER__
#-D__MEASURE_FLAG_NODES__ 
#-D__MEASURE_BUILD_NLIST__
#-D__MEASURE_DA__
#-D__DEBUG__ 
#-D__DEBUG_DA_PUBLIC__
# -D__DEBUG_DA__
# -D__DEBUG_DAQ__
# -D__DEBUG_DA_NLIST__
#-D__DEBUG_DAQ_PUBLIC__
#-D__DEBUG_TN__ 
#-D__DEBUG_OCT__  
#-D__DEBUG_PAR__ 
#-D__DEBUG_MG__ 
#-D__USE_AGV_FOR_BAL__ 
#-D__MEASURE_BPART_COMM__
#-D__MEASURE_BAL_COMM__ 
#-Wall -Wold-style-cast -Woverloaded-virtual -Weffc++ -Wp64

GC = g++

include ./makefileCore

