
/**
  @file testUtils.txx
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  */

#include "dtypes.h"
#include "parUtils.h"
#include "seqUtils.h"
#include <iostream>

namespace seq { 
  namespace test {

    template <typename T>
      bool isSorted(const std::vector<T > & nodes) {
        for (unsigned int i = 1; i < nodes.size(); i++) {
          if ( nodes[i] < nodes[i-1] ) {
            std::cout<<"\n Local Sort Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl<<std::endl;
            return false;
          }
        }
        return true;
      }

    template <typename T>
      bool isSorted(T* nodes, unsigned int sz) {
        for (unsigned int i = 1; i < sz; i++) {
          if ( nodes[i] < nodes[i-1] ) {
            std::cout<<"\n Local Sort Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl<<std::endl;
            return false;
          }
        }
        return true;
      }

    template <typename T>
      bool isUniqueAndSorted(const std::vector<T > & nodes) {
        for (unsigned int i = 1; i < nodes.size(); i++) {
          if ( nodes[i] <= nodes[i-1] ) {
            std::cout<<"\n Local Sort+Unique Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl;
            return false;
          }
        }
        return true;
      }

    template <typename T>
      bool isUniqueAndSorted(T* nodes, unsigned int sz) {
        for (unsigned int i = 1; i < sz; i++) {
          if ( nodes[i] <= nodes[i-1] ) {
            std::cout<<"\n Local Sort+Unique Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl;
            return false;
          }
        }
        return true;
      }

  }//end namespace
}//end namespace

namespace par { 
  namespace test {

    template <typename T>
      bool isSorted(const std::vector<T > &nodes, MPI_Comm comm) {
        bool localPassed = seq::test::isSorted<T>(nodes);
        bool allLocalsPassed;

        par::Mpi_Allreduce<bool>(&localPassed, &allLocalsPassed, 1,
            par::Mpi_datatype<bool>::LAND(), comm);

        if(allLocalsPassed) {
          bool failedParCheck = false;

          MPI_Comm new_comm;
          par::splitComm2way(nodes.empty(), &new_comm, comm);

          if(!nodes.empty()) {
            int rank;
            int npes;
            MPI_Request request;
            MPI_Status status;
            MPI_Comm_rank(new_comm,&rank);
            MPI_Comm_size(new_comm,&npes);

            //Send last to the next proc.
            T end = nodes[nodes.size()-1];				
            if(rank < (npes -1)) {
              par::Mpi_Issend<T>(&end, 1, rank+1, 0, new_comm, &request );
            }

            T prev, me;
            me = nodes[0];

            //Recv prev from the prev proc.d
            if(rank) {
              par::Mpi_Recv<T>( &prev, 1, rank-1, 0, new_comm, &status );
              if(prev > me) {
                failedParCheck = true;			 	
              }
            }

            if(rank < (npes-1)) {
              MPI_Status statusWait;
              MPI_Wait(&request, &statusWait);
            }
          }//end if nodes not empty	

          bool anyProcFailed;
          par::Mpi_Allreduce<bool>(&failedParCheck, &anyProcFailed, 1,
              par::Mpi_datatype<bool>::LOR(), comm);

          return (!anyProcFailed);
        }

        return allLocalsPassed;     
      }
    template <typename T>
      bool isUniqueAndSorted(const std::vector<T > &nodes, MPI_Comm comm) {
        bool localPassed = seq::test::isUniqueAndSorted<T>(nodes);
        bool allLocalsPassed;

        par::Mpi_Allreduce<bool>(&localPassed, &allLocalsPassed, 1,
            par::Mpi_datatype<bool>::LAND(), comm);

        if(allLocalsPassed) {
          bool failedParCheck = false;

          MPI_Comm new_comm;
          par::splitComm2way(nodes.empty(), &new_comm, comm);

          if(!nodes.empty()) {
            int rank;
            int npes;
            MPI_Request request;
            MPI_Status status;
            MPI_Comm_rank(new_comm,&rank);
            MPI_Comm_size(new_comm,&npes);

            //Send last to the next proc.
            T end = nodes[nodes.size()-1];				
            if(rank < (npes-1)) {
              par::Mpi_Issend<T>(&end, 1, rank+1, 0, new_comm, &request );
            }

            T prev, me;
            me = nodes[0];

            //Recv prev from the prev proc.
            if(rank) {
              par::Mpi_Recv<T>( &prev, 1, rank-1, 0, new_comm, &status );
              if(prev >= me) {
                failedParCheck = true;			 	
              }
            }

            if(rank < (npes-1)) {
              MPI_Status statusWait;
              MPI_Wait(&request, &statusWait);
            }
          }//end if nodes not empty	

          bool anyProcFailed;
          par::Mpi_Allreduce<bool>(&failedParCheck, &anyProcFailed, 1,
              par::Mpi_datatype<bool>::LOR(), comm);

          return (!anyProcFailed);
        }

        return allLocalsPassed;     
      }

  }//end namespace
}//end namespace


