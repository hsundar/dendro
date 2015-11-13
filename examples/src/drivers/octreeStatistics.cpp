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
    for (int i=0; i<q; ++i) {
        num_faces[i] = calculateBoundaryFaces(ghosted.begin()+partitions[i], ghosted.begin()+partitions[i+1]);
        num_elems[i] = partitions[i+1] - partitions[i];
    }




   int  faces_min=INT_MAX;
   int  faces_max=0;
   int  faces_sum=0;

   double max_min_r=INT_MAX;

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




//@author: Milinda Fernando.
// This function is to calculate the boundary faces
// Assume that the given octree vector is sorted.
int calculateBoundaryFaces(const std::vector<ot::TreeNode> & mesh, int q) {

    ot::TreeNode first=mesh[0];
    ot::TreeNode last=mesh[mesh.size()-1];
    ot::TreeNode R (first.getDim(), first.getMaxDepth());


    int com_size=q;
    int mesh_nodes=mesh.size();
    assert(mesh_nodes>q);
    int local_mesh_size=mesh_nodes/q;
    int boundary_faces[q];



    int found_pt;
    int num_boundary_faces=0;

    ot::TreeNode top;
    ot::TreeNode bottom;
    ot::TreeNode left;
    ot::TreeNode right;
    ot::TreeNode front;
    ot::TreeNode back;
    int temp=0;
    int begin=0;
    int end=0;
    unsigned long total_boundary_faces=0;

    unsigned long min_faces=1<<30;
    unsigned long max_faces=0;

    for (int j=0;j<q;j++){

        boundary_faces[j]=0;
        begin=j*local_mesh_size;

        temp=begin;
        end=begin + local_mesh_size;

        if(end+local_mesh_size>mesh_nodes)
            end=mesh_nodes;


        num_boundary_faces=0;

        while(temp<end)//for(int i=0;i<mesh.size();i++)
        {

            if(mesh[temp].getX()!=0 && mesh[temp].getY()!=0 && mesh[temp].getZ()!=0 ) {
                top = mesh[temp].getTop();
                bottom = mesh[temp].getBottom();
                left = mesh[temp].getLeft();
                right = mesh[temp].getRight();
                front = mesh[temp].getFront();
                back = mesh[temp].getBack();

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], top, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                // std::cout<<"top:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(top)))) {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], bottom, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                //std::cout<<"bottom:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(bottom)))) {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;


                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], left, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                //std::cout<<"left:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(left)))) {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], right, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                //std::cout<<"right:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(right)))) {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

                    num_boundary_faces++;
                }
                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], front, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                //std::cout<<"front:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(front)))) {
//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], back, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                //std::cout<<"back:"<<found_pt<<std::endl;
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(back)))) {

//        if((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top)))
//          std::cout<<"end condition"<<std::endl;

                    num_boundary_faces++;
                }

            }
            temp++;
            if(temp%local_mesh_size==0)
            {
                break;
            }

        }
        boundary_faces[j]=num_boundary_faces;
        total_boundary_faces += num_boundary_faces;
        if (min_faces > num_boundary_faces) min_faces = num_boundary_faces;
        if (max_faces < num_boundary_faces) max_faces = num_boundary_faces;

        //std::cout<<"q:"<<j<<"\t "<<"number of boundary faces:"<<boundary_faces[j]<<std::endl;

    }

    unsigned long global_sum, global_max, global_min;

    MPI_Reduce(&total_boundary_faces, &global_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&max_faces, &global_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&min_faces, &global_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

    return total_boundary_faces;
}




// Contains the functions to calculate the mesh statistics.
// This function is to calculate the boundary faces
// Assume that the given octree vector is sorted.
void calculateBoundaryFaces(const std::vector<ot::TreeNode> & mesh, int q,double* stat) {


    ot::TreeNode first=mesh[0];
    ot::TreeNode last=mesh[mesh.size()-1];
    ot::TreeNode R (first.getDim(), first.getMaxDepth());
    int com_size=0;
    MPI_Comm_size(MPI_COMM_WORLD,&com_size);
    int mesh_nodes=mesh.size();
    stat[0]=0;  //min faces
    stat[1]=0;  //max faces
    stat[2]=0;  // mean faces

//    stat[3]=0;  //min SVR : Surface to Volume Ratio
//    stat[4]=0;  //max SVR
//    stat[5]=0;  // mean SVR

//    std::cout<<"NUmber of Fake Proc:"<<q<<std::endl;
//    std::cout<<"NUmber of mesh nodes:"<<mesh_nodes<<std::endl;
    assert(mesh_nodes>q);

    int local_mesh_size=mesh_nodes/q;
    int boundary_faces[q];

    int found_pt;
    int num_boundary_faces=0;

    ot::TreeNode top;
    ot::TreeNode bottom;
    ot::TreeNode left;
    ot::TreeNode right;
    ot::TreeNode front;
    ot::TreeNode back;
    int temp=0;
    int begin=0;
    int end=0;
    unsigned long total_boundary_faces=0;

    unsigned long min_faces=1<<28;
    unsigned long max_faces=0;



    for (int j=0;j<q;j++){

        boundary_faces[j]=0;
        begin=j*local_mesh_size;

        temp=begin;
        end=begin + local_mesh_size;

        if(end+local_mesh_size>mesh_nodes)
            end=mesh_nodes;


        num_boundary_faces=0;

        while(temp<end)//for(int i=0;i<mesh.size();i++)
        {

            if(mesh[temp].getX()!=0 && mesh[temp].getY()!=0 && mesh[temp].getZ()!=0 ) {

                top = mesh[temp].getTop();
                bottom = mesh[temp].getBottom();
                left = mesh[temp].getLeft();
                right = mesh[temp].getRight();
                front = mesh[temp].getFront();
                back = mesh[temp].getBack();

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], top, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(top)))) {
                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], bottom, std::less<ot::TreeNode>()) -
                            &mesh[begin]);

                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(bottom)))) {
                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], left, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(left)))) {
                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], right, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(right)))) {
                    num_boundary_faces++;
                }
                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], front, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(front)))) {
                    num_boundary_faces++;
                }

                found_pt = (std::lower_bound(&mesh[begin], &mesh[end - 1], back, std::less<ot::TreeNode>()) -
                            &mesh[begin]);
                if (found_pt == 0 || ((found_pt == (end - 1 - begin)) && (!mesh[end - 1].isAncestor(back)))) {


                    num_boundary_faces++;
                }

            }
            temp++;
            if(temp%local_mesh_size==0)
            {
                break;
            }

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

    ot::TreeNode top;
    ot::TreeNode bottom;
    ot::TreeNode left;
    ot::TreeNode right;
    ot::TreeNode front;
    ot::TreeNode back;

    //std::cout<<"CalCulate Boundary Loop started"<<std::endl;
    for(std::vector<ot::TreeNode>::const_iterator it =beg ; *it!=*end ;it++)
    {

        top=(*(it)).getTop();
        bottom=(*(it)).getBottom();
        left=(*(it)).getLeft();
        right=(*(it)).getRight();
        front=(*(it)).getFront();
        back=(*(it)).getBack();

        if(it->getX()==0 || it->getY()==0 || it->getZ()==0)
            continue;

       // if(top.getX() && top.getY() && top.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && top < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(top)))) {
                num_boundary_faces++;
            }
        //}

        //if(bottom.getX() && bottom.getY() && bottom.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && bottom < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(bottom)))) {
                num_boundary_faces++;
            }
        //}

        //if(left.getX() && left.getY() && left.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && left < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(left)))) {
                num_boundary_faces++;
            }
        //}

        //if(right.getX() && right.getY() && right.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && right < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(right)))) {
                num_boundary_faces++;
            }
        //}

        //if(front.getX() && front.getY() && front.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && front < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(front)))) {
                num_boundary_faces++;
            }
        //}

        //if(back.getX() && back.getY() && back.getZ()>0) {
            found_pt = beg;
            while (found_pt < end && back < *found_pt) found_pt++;
            if (found_pt == beg || ((found_pt == end) && (!(*it).isAncestor(back)))) {
                num_boundary_faces++;
            }
        //}


    }



return num_boundary_faces;


}

