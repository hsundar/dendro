
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
      bool isSorted_all_pairs(const std::vector<T > & nodes) {
        for (unsigned int i = 0; i < nodes.size(); i++) {
          for(unsigned int j=i+1;j<nodes.size();j++)
          {
            if(nodes[i]>nodes[j]) {
              std::cout << "\n Local Sort_pair Check failed for: " << nodes[i] << " and " << nodes[j] << std::endl << std::endl;
              return false;
            }
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
            std::cout<<"nodes[i] < nodes[i-1]:"<<(nodes[i] < nodes[i-1])<<std::endl;
            std::cout<<"nodes[i] > nodes[i-1]:"<<(nodes[i] > nodes[i-1])<<std::endl;
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
            std::cout<<"nodes[i] <= nodes[i-1]:"<<(nodes[i] <= nodes[i-1])<<std::endl;
            std::cout<<"nodes[i] >= nodes[i-1]:"<<(nodes[i] >= nodes[i-1])<<std::endl;
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

//@author: Milinda Fernando
// School of Computing, University  of Utah
      template <typename T>
      bool  isComplete (const std::vector<T>& nodes, MPI_Comm comm) {

        int vol=0;

        int rank;
        MPI_Comm_rank(comm, &rank);
        int maxDepth=0;
        int maxDepth_g=0;

        if(nodes.size()!=0) {


          maxDepth = nodes[0].getMaxDepth();
          int dim = nodes[0].getDim();

          vol = 0;
          int len = 0;
          for (int i = 0; i < nodes.size(); i++) {
            len = 1u << (maxDepth - nodes[i].getLevel());
            if (i < (nodes.size() - 1) && (nodes[i].isAncestor(nodes[i + 1]))) {
              return false;
            }

            vol += len * len * len;
          }

        }else{
          vol=0;
        }


        int g_vol=0;
        par::Mpi_Reduce(&vol,&g_vol,1,MPI_SUM,0,comm);
        MPI_Allreduce(&maxDepth,&maxDepth_g,1,MPI_INT,MPI_MAX,comm);

        if(!rank)
        {
          assert(maxDepth_g!=0);
          std::cout<<"MaxDepth of the Octree:"<<maxDepth_g<<std::endl;
          int max_len=1u<<maxDepth_g;
          if(g_vol==(max_len*max_len*max_len))
          {
            return true;
          }else
          {
            std::cout<<"Volume of the Complete Octree:"<<(max_len*max_len*max_len)<<std::endl;
            std::cout<<"Computed octree volume:"<<g_vol<<std::endl;
            return false;
          }

        }

      }

  }//end namespace
}//end namespace


