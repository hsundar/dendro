/*
 * @author: Milinda Fernando
 * Schhol of Computing University of Utah
 *
 * This file align two octants and compare the differences in two octants.
 *
 *
 * */




#include "mpi.h"
#include "petsc.h"
#include "sys.h"

#include "parUtils.h"
#include "oda.h"


#include <execinfo.h>

#include <cxxabi.h>

#include "externVars.h"
#include "octreeStatistics.h"
#include "testUtils.h"



//Don't want time to be synchronized. Need to check load imbalance.
#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif


void split(std::vector<ot::TreeNode>& t, const ot::TreeNode& splitter,MPI_Comm  comm)
{

    int rank,size;
    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);

    std::vector<ot::TreeNode> sendNodes;
    std::vector<ot::TreeNode> recvNodes;

    long sendCnt =0;
    long recvCnt =0;


    MPI_Request send_req;
    //MPI_Request recv_req;
    MPI_Status  sts;

    if(rank!=size-1) {

        while ((sendCnt < t.size()) && t[sendCnt] < splitter) sendCnt++;
        if(sendCnt!=0) {
            for (int i = sendCnt; i < t.size(); i++) {
                sendNodes.push_back(t[i]);
                t.erase(t.begin()+i);
            }



            sendCnt = t.size() - sendCnt+1;
            sendNodes.resize(sendCnt);

            par::Mpi_Isend(sendNodes.data(), sendCnt, rank + 1, 1, comm, &send_req);
        }

    }

    long sendCnt_all[size];
    par::Mpi_Allgather(&sendCnt,sendCnt_all,1,comm);
//    for(int i=0;i<size;i++)
//        std::cout<<"SendCnt_all "<<i<<" :"<<sendCnt_all[i]<<std::endl;



    if(rank!=0) {
        recvCnt = sendCnt_all[rank - 1];
        if(recvCnt>0)
        {
            recvNodes.resize(recvCnt);
            par::Mpi_Irecv(recvNodes.data(),recvCnt,rank-1,1,comm,&send_req);
            MPI_Wait(&send_req,&sts);

//            for(int i=0;i<recvNodes.size();i++)
//                std::cout<<"Nodes Recieved for rank: "<<rank<<" :"<<i<<": "<<recvNodes[i]<<std::endl;

            std::vector<ot::TreeNode> t_aligned;
            t_aligned.reserve(recvNodes.size()+t.size());
            t_aligned.insert(t_aligned.end(),recvNodes.begin(),recvNodes.end());
            t_aligned.insert(t_aligned.end(),t.begin(),t.end());
            t=t_aligned;

        }

    }





}




void align_octrees(std::vector<ot::TreeNode>& t1, std::vector<ot::TreeNode>& t2, MPI_Comm comm)
{

    int rank,size;
    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);

    std::vector<ot::TreeNode> t1_cpy=t1;
    std::vector<ot::TreeNode> t2_cpy=t2;

  //  treeNodesTovtk(t1,rank,"t1_aligned");
  //  treeNodesTovtk(t2,rank,"t2_aligned");


    int l_sz_1=t1_cpy.size();
    int l_sz_2=t2_cpy.size();


    int g_sz_1;
    int g_sz_2;

    par::Mpi_Allreduce(&l_sz_1,&g_sz_1,1,MPI_SUM,comm);
    par::Mpi_Allreduce(&l_sz_2,&g_sz_2,1,MPI_SUM,comm);

    int g_index_1=0;
    int g_index_2=0;

    par::Mpi_Scan(&l_sz_1,&g_index_1,1,MPI_SUM,comm);
    par::Mpi_Scan(&l_sz_2,&g_index_2,1,MPI_SUM,comm);

    // sort the octrees using Hilbert or Morton.
    ot::TreeNode max_1;
    ot::TreeNode max_2;

    ot::TreeNode maxs[2];
    maxs[0]=t1.back();
    maxs[1]=t2.back();

    std::vector<ot::TreeNode> max_all(2*size);



    par::Mpi_Allgather(maxs, max_all.data(), 2, comm);
    // par::Mpi_Allgather(&max_2,max_all+size,1,comm);

        if(!rank)
        for(int i=0;i<2*size;i++)
            std::cout<<"max_all :"<<i<<": "<<max_all[i]<<std::endl;


    std::sort(max_all.begin(), max_all.end());

    long total_sz = g_sz_1 + g_sz_2;

    long *lc_rank = new long[2*size];
    long *gb_rank = new long[2*size];

    // !!================================
    long less_than_A=0, less_than_B=0;

    for (int i=0; i<max_all.size(); ++i) {
        while ((less_than_A < t1.size()) && t1[less_than_A] <= max_all[i]) less_than_A++;
        while ((less_than_B < t2.size()) && t2[less_than_B] <= max_all[i]) less_than_B++;

        lc_rank[i] = less_than_A + less_than_B;
    }
    par::Mpi_Allreduce(lc_rank, gb_rank,2*size, MPI_SUM, comm);


//    if(!rank)
//        for(int i=0;i<2*size;i++)
//            std::cout<<"gb_rank:"<<i<<":"<<gb_rank[i]<<std::endl;


    // !!================================

    std::vector<ot::TreeNode> splitter_all;

    int idx=0;
    for (int i=1; i<size; ++i) {
       long opt_split = ((long)i*total_sz)/size;
       while ( opt_split < gb_rank[idx] ) idx++;

       ot::TreeNode splitter = (opt_split - gb_rank[idx])<(gb_rank[idx+1]-opt_split)? max_all[idx]:max_all[idx+1];
       splitter_all.push_back(splitter);

    }

    if(!rank)
        for(int i=0;i<splitter_all.size();i++)
        {
            std::cout<<"Splitter of "<<i<<" is : "<<splitter_all[i]<<std::endl;
        }


    split(t1,splitter_all[rank],comm);
    split(t2,splitter_all[rank],comm);

    std::vector<ot::TreeNode> diff_A;
    std::vector<ot::TreeNode> diff_B;

    for(int i=0;i<t1.size();i++)
    {
        int lb=0;
        while((lb<t2.size()) && t2[lb]<t1[i]) lb++;

        if(t2[lb]!=t1[i]) {
            diff_B.push_back(t2[lb]);
            diff_A.push_back(t1[lb]);
        }
    }

    if(diff_A.size()>0)
        treeNodesTovtk(diff_A,rank,"diff_A");
    if(diff_B.size()>0)
        treeNodesTovtk(diff_B,rank,"diff_B");


    treeNodesTovtk(t1,rank,"t1_aligned");
    treeNodesTovtk(t2,rank,"t2_aligned");



//      treeNodesTovtk(diff,rank,"diff");


}




int main(int argc, char **argv)
{

    int size, rank;
    bool incCorner = 1;
    char ptsFileName[256];

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);



   // std::cout<<"Com Size:"<<size<<std::endl;

    std::vector<ot::TreeNode> oct_M;
    std::vector<ot::TreeNode> oct_H;

    sprintf(ptsFileName, "%s_H_%d_%d.oct", argv[1], rank, size);
    ot::readNodesFromFile_binary(ptsFileName,oct_H);

    sprintf(ptsFileName, "%s_M_%d_%d.oct", argv[1], rank, size);
    ot::readNodesFromFile_binary(ptsFileName,oct_M);

//    if(!rank)
//    {
//        std::cout<<"Input files read complete"<<std::endl;
//    }

    std::vector<ot::TreeNode> sorted_octs;
    par::sampleSort(oct_H,sorted_octs,MPI_COMM_WORLD);
    oct_H=sorted_octs;

//    if(!rank)
//    {
//        std::cout<<"Input Sorted"<<std::endl;
//    }

    align_octrees(oct_H,oct_M,MPI_COMM_WORLD);








    MPI_Finalize();

}

