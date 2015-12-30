//@author: Milinda Fernando.
// School of Computing
// University of Utah

#include "octreeStatistics.h"

//@author: Hari Sundar
void flexiblePartitionCalculation(std::vector<ot::TreeNode>& balOct,double slack,int q,MPI_Comm comm)
{

    int rank,size;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);


    unsigned long localSz=balOct.size(), globalSz;
    MPI_Allreduce(&localSz, &globalSz, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    int slackCnt = slack*globalSz/size/q;
    int part_size = localSz/q;

    if(!rank) std::cout << "slack size is " << slackCnt << " octants" << std::endl;

    assert(slackCnt < part_size);


    //----------------------------------------------------------------------
    //   FLEX
    //----------------------------------------------------------------------

    // 1. each process sends slackCnt octs to next/prev.
    std::vector<ot::TreeNode> ghosted(localSz+2*slackCnt);
    MPI_Status status;

    ghosted.insert(ghosted.begin()+slackCnt, balOct.begin(), balOct.end());


    int prev = (rank)?rank-1:MPI_PROC_NULL;
    int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

    ot::TreeNode* sendPrev = &(*(balOct.begin()));
    ot::TreeNode* sendNext = &(balOct[balOct.size() - slackCnt]);

    ot::TreeNode* recvPrev = &(*(ghosted.begin()));
    ot::TreeNode* recvNext = &(*(ghosted.begin()+localSz+slackCnt));

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

    // Have the extra octants ...
    // 2. compute partitions.
    int* partitions   = new int[q+1];
    int* num_faces    = new int[q+1];
    int* num_elems    = new int[q+1];
    int faces;

    // partition 0
    partitions[0] = slackCnt;
    if (!rank) {
        // find first
        num_faces[0] = 6*localSz; // calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
        for (int j = 0; j < 2*slackCnt; ++j) {
            faces = calculateBoundaryFaces(ghosted.begin()+j, ghosted.begin()+slackCnt+part_size);
            if (faces < num_faces[0]) {
                num_faces[0] = faces; // we will update this later ...
                partitions[0] = j;
            }
        }
    } else {
        num_faces[0] = 6*localSz; // calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
        for (int j = 0; j < slackCnt; ++j) {
            faces = calculateBoundaryFaces(ghosted.begin()+slackCnt+j, ghosted.begin()+slackCnt+part_size);
            if (faces < num_faces[0]) {
                num_faces[0] = faces; // we will update this later ...
                partitions[0] = slackCnt+j;
            }
        }
    }
    std::vector<ot::TreeNode>::const_iterator first, last;

    // local partitions
    for (int i=1; i<q; ++i) {
        first = ghosted.begin() + partitions[i-1];
        num_faces[i] = 6*localSz;
        last = ghosted.begin() + i*part_size;
        for (int j = 0; j < 2*slackCnt; ++j) {
            faces = calculateBoundaryFaces(first, last + j);
            if (faces < num_faces[i]) {
                num_faces[i] = faces;
                partitions[i] = i*part_size + j;
            }
        }
    }
    // partition q
    if (rank == size-1) {
        first = ghosted.begin() + partitions[q-1];
        num_faces[q] = 6*localSz;
        last = ghosted.begin() + q*part_size;
        for (int j = 0; j < slackCnt; ++j) {
            faces = calculateBoundaryFaces(first, last + j);
            if (faces < num_faces[q]) {
                num_faces[q] = faces;
                partitions[q] = q*part_size + j;
            }
        }
    } else {
        first = ghosted.begin() + partitions[q-1];
        num_faces[q] = 6*localSz;
        last = ghosted.begin() + q*part_size;
        for (int j = 0; j < 2*slackCnt; ++j) {
            faces = calculateBoundaryFaces(first, last + j);
            if (faces < num_faces[q]) {
                num_faces[q] = faces;
                partitions[q] = q*part_size + j;
            }
        }
    }

    // @milinda now update num_faces
    char vtkFilename[256];
    sprintf(vtkFilename,"%s_%f_%d","flex_part",slack,rank);
    std::vector<ot::TreeNode>flexParts[q+1];
    for (int i=0; i<q; ++i) {
        num_faces[i] = calculateBoundaryFaces(ghosted.begin()+partitions[i], ghosted.begin()+partitions[i+1]);
        num_elems[i] = partitions[i+1] - partitions[i];
        //flexParts[i].insert(ghosted.begin()+partitions[i], ghosted.begin()+partitions[i+1]);
        for(int k=partitions[i];k<=partitions[i+1];k++)
            flexParts[i].push_back(*(ghosted.begin()+k));

        treeNodesTovtk(flexParts[i],i,vtkFilename);

    }





   int  faces_min=INT8_MAX;
   int  faces_max=0;
   int  faces_sum=0;

   double max_min_r=INT8_MAX;

   int  element_min=0;
   int  element_max=0;
   int elems_sum=0;

    for(int i=0;i<(q);i++)
    {
        //std::cout<<"Face Status: "<<face_status[i]<<std::endl;

        if( faces_min>num_faces[i])
            faces_min=num_faces[i];

        if(faces_max<num_faces[i])
            faces_max=num_faces[i];


        if(element_min>num_elems[i])
            element_min=num_elems[i];

        if(element_max<num_elems[i]);
            element_max=num_elems[i];

        faces_sum=+num_faces[i];
        elems_sum=+num_elems[i];
    }

    max_min_r=(double)faces_max/(double)faces_min;

    if(!rank)
        std::cout<<"Flexible part calculation completed"<<std::endl;
    //std::cout<<"Rank:"<<rank<<" Faces:"<<faces<<std::endl;

    int faces_min_g;
    int faces_max_g;
    int faces_sum_g;

    int elems_min_g;
    int elems_max_g;
    int elems_sum_g;

    double max_min_r_g;

    double mean_num_faces;
    double mean_num_elems;

    MPI_Allreduce(&faces_min,&faces_min_g,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(&faces_max,&faces_max_g,1,MPI_INT,MPI_MAX,comm);
    MPI_Allreduce(&faces_sum,&faces_sum_g,1,MPI_INT,MPI_SUM,comm);
    MPI_Allreduce(&max_min_r,&max_min_r_g,1,MPI_DOUBLE,MPI_MIN,comm);

    MPI_Allreduce(&element_min,&elems_min_g,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(&element_max,&elems_max_g,1,MPI_INT,MPI_MAX,comm);
    MPI_Allreduce(&elems_sum,&elems_sum_g,1,MPI_INT,MPI_SUM,comm);

    mean_num_faces=(double)faces_sum_g/(size*q);
    mean_num_elems=(double)elems_sum_g/(size*q);
   // MPI_Allgather(&faces,1,MPI_INT,faces_all,1,MPI_INT,comm);

//  if(!rank)
//    for(int i=0;i<size;i++)
//      std::cout<<"faces_all:"<<i<<" :"<<faces_all[i]<<std::endl;


    if (!rank) {
        std::cout << BLU << "===============================================" << NRM << std::endl;
        std::cout << RED " Boundary Surfaces (min):"<<faces_min_g<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (max):"<<faces_max_g<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (mean):"<<mean_num_faces<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (max/min parition wise):"<<max_min_r_g<< NRM << std::endl;

        std::cout << YLW " Partition elements (min):"<<elems_min_g<< NRM << std::endl;
        std::cout << YLW " Partition elements (max):"<<elems_max_g<< NRM << std::endl;
        std::cout << YLW " Partition elements (mean):"<<mean_num_elems<< NRM << std::endl;

        std::cout << BLU << "===============================================" << NRM << std::endl;
    }




    if (!rank) {
        std::cout << GRN << "Finalizing ..." << NRM << std::endl;
    }

    delete [] partitions;
    delete [] num_faces;


}

// Contains the functions to calculate the mesh statistics.
// This function is to calculate the boundary faces
// Assume that the given octree vector is sorted.
void calculateBoundaryFaces(const std::vector<ot::TreeNode> &partition, int q, double* stat) {


    int com_size=0;
    MPI_Comm_size(MPI_COMM_WORLD,&com_size);
    int local_sz= partition.size();

//    stat[3]=0;  //min SVR : Surface to Volume Ratio
//    stat[4]=0;  //max SVR
//    stat[5]=0;  // mean SVR

    stat[0]=0;  //min faces
    stat[1]=0;  //max faces
    stat[2]=0;  // mean faces

    assert(local_sz>q);

    int boundary_faces[q];
    int found_pt;
    int num_boundary_faces=0;

    int temp=0;
    int begin=0;
    int end=0;
    unsigned long total_boundary_faces=0;

    unsigned long min_faces=1u<<25;
    unsigned long max_faces=0;

    std::vector<ot::TreeNode> neighbourOct;
    ot::TreeNode tmp;
    int surf_per_oct=0;
    for (int j=0;j<q;j++){

        boundary_faces[j]=0;
        begin=j*local_sz/q;
        end=(j+1)*local_sz/q;
        num_boundary_faces=0;

        for(int temp=begin;temp<end;temp++)
        {

            tmp=partition[temp].getTop();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);

            tmp=partition[temp].getBottom();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);

            tmp=partition[temp].getRight();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);

            tmp=partition[temp].getLeft();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);

            tmp=partition[temp].getFront();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);

            tmp=partition[temp].getBack();
            if(!tmp.isRoot())
                neighbourOct.push_back(tmp);


            surf_per_oct=0;
            for(int k=0;k<neighbourOct.size();k++)
            {
                tmp=neighbourOct[k];

//                int found_pt=std::lower_bound(partition.begin()+begin,partition.begin()+begin+end-1,tmp)-partition.begin()-begin;
//                if((found_pt==0 && partition[found_pt]!=tmp )|| (found_pt==(end-1) && !partition[found_pt].isAncestor(tmp) && partition[found_pt]<tmp))
//                {
//
//                    surf_per_oct++;
//                    if(!(tmp<partition[begin] || (tmp>partition[end-1] && !partition[end-1].isAncestor(tmp)))) {
//                        std::cout<<"Condition Mismatch: "<<found_pt<<"\t tmp:"<<tmp<<"\t begin:"<<partition[begin]<<"end: "<<partition[end-1]<<std::endl;
//                    }
//
//
//                }

                if(tmp<partition[begin] || (tmp>partition[end-1]) )
                {
                    surf_per_oct++;
                }

            }
            neighbourOct.clear();

            num_boundary_faces=num_boundary_faces+surf_per_oct;


        }

        boundary_faces[j]=num_boundary_faces;
        total_boundary_faces += num_boundary_faces;
        if (min_faces > num_boundary_faces) min_faces = num_boundary_faces;
        if (max_faces < num_boundary_faces) max_faces = num_boundary_faces;



    }

    unsigned long global_sum, global_max, global_min;

    MPI_Reduce(&total_boundary_faces, &global_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&max_faces, &global_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&min_faces, &global_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

    stat[0]=global_min;
    stat[1]=global_max;
    stat[2]=global_sum/(double)(q*com_size);


}

int calculateBoundaryFaces(const std::vector<ot::TreeNode>::const_iterator &beg, const std::vector<ot::TreeNode>::const_iterator &end)
{


    std::vector<ot::TreeNode>::const_iterator found_pt;
    int num_boundary_faces=0;
    ot::TreeNode tmp;
    std::vector<ot::TreeNode> neighbourOct;
    //std::cout<<"CalCulate Boundary Loop started"<<std::endl;
    for(std::vector<ot::TreeNode>::const_iterator it =beg ; *it!=*end ;it++)
    {

        neighbourOct.clear();

        tmp=(*it).getTop();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=(*it).getBottom();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=(*it).getRight();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=(*it).getLeft();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=(*it).getFront();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);

        tmp=(*it).getBack();
        if(!tmp.isRoot())
            neighbourOct.push_back(tmp);



        for(int k=0;k<neighbourOct.size();k++)
        {
            tmp=neighbourOct[k];
            if(tmp<*beg || (tmp>*end))
            {
                num_boundary_faces++;
            }
        }




    }



return num_boundary_faces;


}

