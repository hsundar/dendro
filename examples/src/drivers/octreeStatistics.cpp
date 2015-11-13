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
    MPI_Allreduce(&localSz, &globalSz, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);

    int slackCnt = slack*globalSz/size/q;

    if(!rank) std::cout << "slack size is " << slackCnt << " octants" << std::endl;

    //----------------------------------------------------------------------
    //   FLEX
    //----------------------------------------------------------------------

    // 1. each process sends slackCnt octs to next/prev.
    std::vector<ot::TreeNode> slack_next(slackCnt), slack_prev(slackCnt);
    MPI_Status status;

    int prev = (rank)?rank-1:MPI_PROC_NULL;
    int next = (rank+1==size)?MPI_PROC_NULL:rank+1;

    ot::TreeNode* sendPrev = &(*(balOct.begin()));
    ot::TreeNode* sendNext = &(balOct[balOct.size() - slackCnt]);

    ot::TreeNode* recvPrev = &(*(slack_prev.begin()));
    ot::TreeNode* recvNext = &(*(slack_next.begin()));

//    if(!rank)
//        std::cout<<"TreeNode Communication started"<<std::endl;

    //send to prev

    if (rank)
        par::Mpi_Sendrecv<ot::TreeNode>(sendPrev, slackCnt, prev, 0,
                                        recvPrev, slackCnt, prev, 0,
                                        comm, &status);
    // send to next
    if (rank != (size-1))
        par::Mpi_Sendrecv<ot::TreeNode>(sendNext, slackCnt, next, 0,
                                        recvNext, slackCnt, next, 0,
                                        comm, &status);

//  for(int i=0;i<slack_prev.size();i++)
//  std::cout<<"Recive Previous:"<<i<<" ::"<<slack_prev[i]<<std::endl;



    // Have the extra octants ...
    // 2. compute partitions.
    int part_size = localSz/q;
    int* partitions   = new int[q+1];
    int* num_faces    = new int[q+1];
    bool * face_status=new bool[q+1];
    int faces;
    int faces_min=0;
    int faces_max=0;
    int faces_sum=0;

    int faces_all[size];
//    if(!rank)
//        std::cout<<"TreeNode Communication Completed"<<std::endl;

    // partition 0

    for(int i=0;i<(q+1);i++)
        face_status[i]= false;


    partitions[0] = slackCnt;

    if (!rank) {
        num_faces[0] = calculateBoundaryFaces(balOct.begin(), balOct.begin()+part_size);
        face_status[0]=true;
        for(int i=0;i<part_size;i++)
            slack_prev.push_back(balOct[i]);
            //slack_prev.insert(slack_prev.end(), balOct.begin(), balOct.begin()+part_size);


        for (int j = 0; j < slackCnt; ++j) {
            faces = calculateBoundaryFaces((slack_prev.begin()+j), slack_prev.end());
            //std::cout<<"j:"<<j<<std::endl;
            //std::cout<<"Faces from rank:"<<rank<<": is :"<<faces<<std::endl;
            //std::cout<<"Rank:"<<rank<<"Slack_Cnt:"<<slackCnt<<std::endl;
            if (faces < num_faces[0]) {
                num_faces[0] = faces;
                partitions[0] = j;
                face_status[0]=true;
            }
        }
        //std::cout<<"Rank 0 finished"<<std::endl;
    }

    std::vector<ot::TreeNode>::const_iterator first, last;
    // local partitions
    for (int i=1; i<q; ++i) {
        first = balOct.begin()+i*part_size;
        last = (i==(q-1))?balOct.end():balOct.begin()+(i+1)*part_size;
        //assert(first<last);
        partitions[i] = i*part_size + slackCnt;
        num_faces[i] = calculateBoundaryFaces(first,last);
        face_status[i]=true;
        first -= slackCnt;
        for (int j = 0; j < slackCnt; ++j) {
            //assert(j<balOct.size());
            faces = calculateBoundaryFaces(first+j, last);
           // std::cout<<"Faces from rank:"<<rank<<": is :"<<faces<<std::endl;
           // std::cout<<"Rank:"<<rank<<"Slack_Cnt:"<<slackCnt<<std::endl;
            if (faces < num_faces[i]) {
                num_faces[i] = faces;
                partitions[i] = j;
                face_status[i]=true;
            }

        }

       // std::cout<<"Rank "<<rank<<" Loop 2 ended"<<std::endl;
    }



    faces_min=INT_MAX;
    faces_max=0;
    faces_sum=0;

//    std::cout << BLU << "===============================================" << NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (min):"<<faces_min<< NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (max):"<<faces_max<< NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (mean):"<<faces_sum<< NRM << std::endl;
//    std::cout << BLU << "===============================================" << NRM << std::endl;




    for(int i=0;i<(q+1);i++)
    {
        //std::cout<<"Face Status: "<<face_status[i]<<std::endl;

        if( face_status[i] && faces_min>num_faces[i])
            faces_min=num_faces[i];

        if(face_status[i] && faces_max<num_faces[i])
            faces_max=num_faces[i];

        if(face_status[i])
            faces_sum=+num_faces[i];
    }

//    std::cout << BLU << "===============================================" << NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (min):"<<faces_min<< NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (max):"<<faces_max<< NRM << std::endl;
//    std::cout << RED " Boundary Surfaces (mean):"<<faces_sum<< NRM << std::endl;
//    std::cout << BLU << "===============================================" << NRM << std::endl;


    if(!rank)
        std::cout<<"Flexible part calculation completed"<<std::endl;
    //std::cout<<"Rank:"<<rank<<" Faces:"<<faces<<std::endl;

    int faces_min_g;
    int faces_max_g;
    int faces_sum_g;



    double mean_num_faces;

    MPI_Allreduce(&faces_min,&faces_min_g,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(&faces_max,&faces_max_g,1,MPI_INT,MPI_MAX,comm);
    MPI_Allreduce(&faces_sum,&faces_sum_g,1,MPI_INT,MPI_SUM,comm);

    mean_num_faces=(double)faces_sum_g/(size*q);

   // MPI_Allgather(&faces,1,MPI_INT,faces_all,1,MPI_INT,comm);

//  if(!rank)
//    for(int i=0;i<size;i++)
//      std::cout<<"faces_all:"<<i<<" :"<<faces_all[i]<<std::endl;


    if (!rank) {
        std::cout << BLU << "===============================================" << NRM << std::endl;
        std::cout << RED " Boundary Surfaces (min):"<<faces_min_g<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (max):"<<faces_max_g<< NRM << std::endl;
        std::cout << RED " Boundary Surfaces (mean):"<<mean_num_faces<< NRM << std::endl;
        std::cout << BLU << "===============================================" << NRM << std::endl;
    }



    // partition q
    partitions[q] = balOct.size();

    //----------------------------------------------------------------------
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

