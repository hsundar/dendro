#include "dynamicPartition.h"

/*
 * Input:
 * @parameter partiton: Partition we need to compute the boundary (partition should be unique and sorted,)
 *
 * Output:
 * @parameter boundary: list of boundary octants (sorted)
 * @parameter bdy_surfaces: Number of boundary surfaces.
 *
 * */
long BoundaryCalculation(std::vector<ot::TreeNode>& partition,int begin,int end, char* bdy_data)
{

    assert(!partition.empty());
    unsigned int dim=partition[0].getDim();
    int maxDepth=partition[0].getMaxDepth();

    assert(begin<partition.size());
    assert(end<=partition.size());

    ot::TreeNode tmp;

    long bdy_surfaces=0;



    for(int i=begin;i<end;i++)
    {

        std::vector<ot::TreeNode> neighbourOct;
        bdy_data[i]=0;
        //neighbour octant returns root if the getBottom, getTop, getLeft, getRight, getFront and getBack functions are not in side the root octant,


        tmp=partition[i].getTop();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=partition[i].getBottom();
        if(!tmp.isRoot())
        neighbourOct.push_back(tmp);

        tmp=partition[i].getRight();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=partition[i].getLeft();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=partition[i].getFront();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=partition[i].getBack();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);


        int surf_per_oct=0;
        for(int j=0;j<neighbourOct.size();j++)
        {

            tmp=neighbourOct[j];
            if(tmp<partition[begin] || (tmp>partition[end-1] && !partition[end-1].isAncestor(tmp)) )
            {
               surf_per_oct++;
            }

        }
        if(surf_per_oct>0)
        {
            int l=BAL_CONST*surf_per_oct;
            bdy_data[i]=l;
            bdy_surfaces=bdy_surfaces+l;
        }


        neighbourOct.clear();

    }


    return  bdy_surfaces;



}

void UpdateBoundaryCalculation(std::vector<ot::TreeNode>& partition,int begin, int end,long& bdy_surfaces,char* bdy_data)
{

    assert(end<(partition.size()));
    assert(!partition.empty());

    ot::TreeNode new_octant=partition[end];

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //std::cout<<"Rank:"<<rank<<" Slack Node:"<<new_octant<<std::endl;

    std::vector<ot::TreeNode> neighbourOct;
    ot::TreeNode tmp;

    tmp=new_octant.getTop();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);

    tmp=new_octant.getBottom();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);

    tmp=new_octant.getLeft();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);


    tmp=new_octant.getRight();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);

    tmp=new_octant.getFront();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);

    tmp=new_octant.getBack();
    if(!tmp.isRoot())
        neighbourOct.push_back(tmp);


    int num_surf_oct=0;
    bdy_data[end]=0;

    for(int j=0;j<neighbourOct.size();j++)
    {

        tmp=neighbourOct[j];

        if(tmp<*(partition.begin()+begin) || tmp>*(partition.begin()+begin+end)) {
            num_surf_oct++;
            bdy_data[end]=bdy_data[end]+BAL_CONST;
            bdy_surfaces=bdy_surfaces+BAL_CONST;

        }else
        {

//            bool state= (tmp<*(partition.begin()+begin) || tmp>*(partition.begin()+begin+end));
//            if(state)
//                std::cout<<"Rank:"<<rank<<" Seraching for: "<<tmp<<std::endl;

            int found_oct=std::lower_bound((partition.begin()+begin),(partition.begin()+begin+end),tmp)-partition.begin();
            //std::cout<<"found_pt:"<<found_oct<<" \t "<<partition[found_oct]<<std::endl;
            int l=0;
            int pow=(2-2*abs(new_octant.getLevel()-partition[found_oct].getLevel()));

            l=(1u<<pow); // this gives us how much we need to reduce,
            assert(pow==0 || pow==2);
//            if(pow!=2 && pow!=0)
//            {
//                std::cout<<"Rank:"<<rank<<" POW:"<<pow<<" new oct_lev:"<<new_octant.getLevel()<<" partition lev:"<<partition[found_oct].getLevel()<<std::endl;
//            }

            if(partition[found_oct].getLevel()<=new_octant.getLevel())
            {
                bdy_data[found_oct]=bdy_data[found_oct]-l;

            }else
            {
                bdy_data[end]=bdy_data[end]-l;

            }
            bdy_surfaces=bdy_surfaces-l;
        }

      //  std::cout<<"Rank:"<<rank<<"\t Number of boundary faces:"<<bdy_surfaces<<std::endl;

    }
    neighbourOct.clear();

    if(bdy_data[end]<0)
        bdy_data[end]=0;



}


/*
 *Input: @parameter: partition: Original partitions, slack: should be [0,1]
 * Output: @parameter: New partition,
 *
 */

void DynamicPartitioning(std::vector<ot::TreeNode>& partition, double slack, MPI_Comm comm)
{

    int rank,size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);


    int local_sz=partition.size();
    int global_sz=0;

    par::Mpi_Allgather(&local_sz,&global_sz,1,comm);

    int slackCnt=global_sz*slack/size;
    if(!rank) std::cout << "slack size is " << slackCnt << " octants" << std::endl;

    assert(slackCnt < local_sz);

    //----------------------------------------------------------------------
    //   FLEX
    //----------------------------------------------------------------------

    // 1. each process sends slackCnt octs to next/prev.
    std::vector<ot::TreeNode> ghosted(local_sz+2*slackCnt);
    MPI_Status status;

    ghosted.insert(ghosted.begin()+slackCnt, partition.begin(), partition.end());

    int prev = (rank)?rank-1:MPI_PROC_NULL;
    int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

    ot::TreeNode* sendPrev = &(*(partition.begin()));
    ot::TreeNode* sendNext = &(partition[partition.size() - slackCnt]);

    ot::TreeNode* recvPrev = &(*(ghosted.begin()));
    ot::TreeNode* recvNext = &(*(ghosted.begin()+local_sz+slackCnt));

    // send to prev
    if (rank)
        par::Mpi_Sendrecv<ot::TreeNode>(sendPrev, slackCnt, prev, 0,
                                        recvPrev, slackCnt, prev, 0,
                                        MPI_COMM_WORLD, &status);
    // send to next
    if (rank != (size-1))
        par::Mpi_Sendrecv<ot::TreeNode>(sendNext, slackCnt, next, 0,
                                        recvNext, slackCnt, next, 0,
                                        MPI_COMM_WORLD, &status);


    std::vector<ot::TreeNode> boundary;
    char* bdy_data=new char[ghosted.size()];
    long local_nm_faces=0;



    local_nm_faces = BoundaryCalculation(ghosted,slackCnt,slackCnt+partition.size(),bdy_data);
    long local_nm_faces_prev=local_nm_faces/4;
    long global_opt=0;

//    for(int i=0;i<bdy_data.size();i++)
//        std::cout<<"Octant:"<<i<<"\t boundary octant data:"<<(int)bdy_data[i]<<std::endl;

    //std::cout<<"Rank:"<<rank<<"num_nm_faces:"<<local_nm_faces<<std::endl;



    long global_nm_faces_min=0;
    long global_nm_faces_max=0;
    long global_nm_faces_sum=0;
    double global_nm_faces_mean=0.0;

    par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_min,1,MPI_MIN,comm);
    par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_max,1,MPI_MAX,comm);
    par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_sum,1,MPI_SUM,comm);

    global_nm_faces_mean=global_nm_faces_sum/(double)size;

    global_opt=global_nm_faces_min;


    if (!rank) {
        std::cout << YLW << "========Dynamic Partitionign START========" << NRM << std::endl;
        std::cout << RED " Boundary Surfaces (min):"<<global_nm_faces_min<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (max):"<<global_nm_faces_max<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (mean):"<<global_nm_faces_mean<< NRM << std::endl;
        std::cout << YLW << "===============================================" << NRM << std::endl;
    }



    int slackCnt_opt=0;
    long global_opt_stat[3];
    for(int i=0;i<slackCnt;i++)
    {
        if(rank<(size-1)) {
            UpdateBoundaryCalculation(ghosted, slackCnt, slackCnt + partition.size() + i, local_nm_faces, bdy_data);
            local_nm_faces_prev=local_nm_faces/4;
        }
        //std::cout<<"Rank:"<<rank<<"\t Updated Boundary faces:"<<local_nm_faces<<std::endl;

        //local_nm_faces=local_nm_faces/BAL_CONST;

        par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_min,1,MPI_MIN,comm);
        par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_max,1,MPI_MAX,comm);
        par::Mpi_Allreduce(&local_nm_faces_prev,&global_nm_faces_sum,1,MPI_SUM,comm);

        global_nm_faces_mean=global_nm_faces_sum/size;

//        if(!rank && global_opt>global_nm_faces_min)
//        {
//            slackCnt_opt=i;
//            global_opt_stat[0]=global_nm_faces_min;
//            global_opt_stat[1]=global_nm_faces_max;
//            global_opt_stat[2]=global_nm_faces_mean;
//        }


        if (!rank) {
            std::cout << YLW << "======= SLACK NODE:"<<ghosted[slackCnt+partition.size()+i]<<"========" << NRM << std::endl;
            std::cout << RED " Boundary Surfaces (min):"<<global_nm_faces_min<< NRM << std::endl;
            std::cout << RED " Boundary Surfaces (max):"<<global_nm_faces_max<< NRM << std::endl;
            std::cout << RED " Boundary Surfaces (mean):"<<global_nm_faces_mean<< NRM << std::endl;
            std::cout << YLW << "===============================================\n" << NRM << std::endl;
        }


    }

//    MPI_Barrier(MPI_COMM_WORLD);
//
//    std::cout<<"Rank:"<<rank<<" Optiaml Slack index:"<<slackCnt_opt<<std::endl;
//
//
//    if (!rank) {
//            std::cout << YLW << "=======OPTIMAL NUMBER OF BOUNDARY SURFACES========" << NRM << std::endl;
//            std::cout << RED " Boundary Surfaces (min):"<<global_opt_stat[0]<< NRM << std::endl;
//            std::cout << RED " Boundary Surfaces (max):"<<global_opt_stat[1]<< NRM << std::endl;
//            std::cout << RED " Boundary Surfaces (mean):"<<global_opt_stat[2]<< NRM << std::endl;
//            std::cout << YLW << "===============================================\n" << NRM << std::endl;
//    }


}
