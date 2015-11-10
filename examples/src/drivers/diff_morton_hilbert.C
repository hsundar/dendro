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

    max_1=t1[t1.size()-1];
    max_2=t2[t2.size()-1];

    ot::TreeNode max_2_all[size];
    par::Mpi_Allgather(&max_1,max_2_all,1,comm);

    std::vector<ot::TreeNode> max_2_all_v(max_2_all,max_2_all+size);

    int np=0;
    np= std::lower_bound(max_2_all_v.begin(),max_2_all_v.end(),max_1)-max_2_all_v.begin();

    int recieve_np[size];
    //MPI_Alltoall(&np,1,MPI_INT,recieve_np,size,MPI_INT,comm);
    par::Mpi_Allgather(&np,recieve_np,1,comm);

//    if(!rank)
//    for(int i=0;i<size;i++)
//        std::cout<<"recieve_np["<<i<<"]:"<<recieve_np[i]<<std::endl;



    int recieve_cnt=0;
    for(int i=0;i<size;i++)
    {
        if(recieve_np[i]==rank) {

            recieve_cnt++;
        }
    }

//    std::cout<<"Rank:"<<rank<<"\t Reciever Cnt:"<<recieve_cnt<<std::endl;

    ot::TreeNode search_key[recieve_cnt];

    MPI_Request  request;

    par::Mpi_Isend(&max_1,1,np,1,comm,&request);

    for(int i=0;i<size;i++)
    {
        if(recieve_np[i]==rank) {
            par::Mpi_Irecv((search_key+i),1,i,1,comm,&request);

        }
    }

    int rank_g_index_1[recieve_cnt];

    for(int i=0;i<recieve_cnt;i++)
    {
        rank_g_index_1[i]=std::lower_bound(t2.begin(),t2.end(),search_key[i])-t2.begin();
        rank_g_index_1 [i]=g_index_2-l_sz_2+rank_g_index_1[i];
    }

    int rank_g_index_1_all[size*recieve_cnt];

    par::Mpi_Allgather(rank_g_index_1,rank_g_index_1_all,recieve_cnt,comm);

    int ideal_ld=std::max(g_sz_1,g_index_2)/size;
    int splitter=0;
    int splitter_all[size];
    for(int i=0;i<recieve_cnt;i++)
    {
        if(rank_g_index_1[i]>splitter)
            splitter=rank_g_index_1[i];
    }

//    for(int i=0;i<recieve_cnt;i++)
//        std::cout<<"Rank of list 1: rank_g_index_1["<<i<<"]:"<<rank_g_index_1[i]<<std::endl;

    par::Mpi_Allgather(&splitter,splitter_all,1,comm);

//    std::cout<<"Splitter of Rank:"<<rank<<" is : "<<splitter_all[rank]<<std::endl;

    int min_sz=std::min(g_index_1,g_index_2);

    int diff_1=0;
    int diff_2=0;

    int diff_1_all[size];
    int diff_2_all[size];


    diff_1=g_index_1-splitter_all[rank];
    diff_2=g_index_2-splitter_all[rank];


    par::Mpi_Allgather(&diff_1,diff_1_all,1,comm);
    par::Mpi_Allgather(&diff_2,diff_2_all,1,comm);

//    if(!rank)
//    for(int i=0;i<size;i++) {
//        std::cout << "Rank:" << i << "  diff_1_all[" << i<< "]:" << diff_1_all[i] << std::endl;
//        std::cout << "Rank:" << i << "  diff_2_all[" << i << "]:" << diff_2_all[i] << std::endl;
//        std::cout<<std::endl;
//
//        std::cout<<"Spliter Rank:"<<i<<"\t splitter_all["<<i<<"]:"<<splitter_all[i]<<std::endl;
//
//    }

    MPI_Request  req;
    MPI_Status sts;
    if(rank!=(size-1) && diff_1>0)
    {


        MPI_Status status;
        ot::TreeNode send_1[diff_1];
        int count=0;
        for(int i=(splitter_all[rank]+(1))-g_index_1+l_sz_1;i<t1.size();i++)
        {
           send_1[count]=t1[i];
           //std::cout<<"Rank :"<<rank<<"  Sending Treenode Index:"<< i<<", which is :"<<send_1[count]<<std::endl;
            count++;
        }

        //for(int i=0;i<diff_1;i++)


        //if(rank!=size-1)
        par::Mpi_Issend(send_1,diff_1,rank+1,1,comm,&req);
        t1.erase(t1.begin()+((splitter_all[rank]+(1))-g_index_1+l_sz_1),t1.end());




    }


    if(rank!=0 && diff_1_all[rank-1]>0) {

        ot::TreeNode recv_1[diff_1_all[rank-1]];
        par::Mpi_Recv(recv_1,diff_1_all[rank-1],rank-1,1,comm,&sts);

//        if(rank==2)
//        for(int i=0;i<diff_1_all[rank-1];i++)
//            std::cout <<"Recived octants:"<<recv_1[i]<<"Dim:"<<recv_1[i].getDim()<<std::endl;


        std::vector<ot::TreeNode> t1_aligned;
        for (int i = 0; i < diff_1_all[rank-1]; i++)
            t1_aligned.push_back(recv_1[i]);

        treeNodesTovtk(t1,rank,"recieved");

        for (int i = 0; i < t1.size(); i++)
            t1_aligned.push_back(t1_cpy[i]);

//        std::cout<<"Rank:"<<rank<<"\tt1_aligned_size:"<<t1.size()<<std::endl;
//        if(rank==2)
//        for(int i=0;i<t1_aligned.size();i++)
//        {
//            std::cout<<"TreeNode in rank:"<<rank<<" is:"<<t1_aligned[i]<<std::endl;
//        }

        t1.clear();
        t1=t1_aligned;

    }

    if(rank!=(size-1) && diff_2>0)
    {


        ot::TreeNode send_2[diff_2];
        int count=0;
        for(int i=splitter_all[rank]+(1)-g_index_2+l_sz_2;i<t2.size();i++){
            //int lc_index=(splitter_all[rank]+(i+1))-g_index_1+l_sz_2;//g_index_2-(splitter_all[rank]+(i+1));
            send_2[count]=t2[i];
            count++;
        }

        par::Mpi_Issend(send_2,diff_2,rank+1,1,comm,&req);
        t2.erase(t2.begin()+((splitter_all[rank]+(1))-g_index_1+l_sz_2),t2.end());




    }

    if(rank!=0 && diff_2_all[rank-1]>0) {

        ot::TreeNode recv_2[diff_2_all[rank-1]];

        par::Mpi_Recv(recv_2,diff_2_all[rank-1],rank-1,1,comm,&sts);
//        for(int i=0;i<diff_1_all[rank-1];i++)
//            std::cout <<"Recived octants:"<<recv_2[i]<<std::endl;

        std::vector<ot::TreeNode> t2_aligned;


        for (int i = 0; i < diff_2_all[rank-1]; i++)
            t2_aligned.push_back(recv_2[i]);

        for (int i = 0; i < t2.size(); i++)
            t2_aligned.push_back(t2[i]);

        t2.clear();
        t2 = t2_aligned;
    }


//    std::cout<<"Rank:"<<rank<<" Original t1_size:"<<t1_cpy.size()<<"\t aligend size:"<<t1.size()<<std::endl;
//    std::cout<<"Rank:"<<rank<<" Original t2_size:"<<t2_cpy.size()<<"\t aligend size:"<<t2.size()<<std::endl;


      std::vector<ot::TreeNode> diff;

      for(int i=0;i<t1.size();i++)
      {
          int lb=std::lower_bound(t2.begin(),t2.end(),t1[i])-t2.begin();
          if(t2[lb]!=t1[i])
              diff.push_back(t1[i]);

      }






      treeNodesTovtk(t1,rank,"t1_aligned");
      treeNodesTovtk(t2,rank,"t2_aligned");
      treeNodesTovtk(diff,rank,"diff");
















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

