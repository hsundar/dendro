//@author: Milinda Fernando.
// School of Computing
// University of Utah

#include "octreeStatistics.h"

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

    unsigned long min_faces=1<<50;
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

            top=mesh[temp].getTop();
            bottom=mesh[temp].getBottom();
            left=mesh[temp].getLeft();
            right=mesh[temp].getRight();
            front=mesh[temp].getFront();
            back=mesh[temp].getBack();

            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], top, std::less<ot::TreeNode>()) - &mesh[begin]);
            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(top))))
            {
                num_boundary_faces++;
            }

            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], bottom, std::less<ot::TreeNode>()) - &mesh[begin]);

            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(bottom))))
            {
                num_boundary_faces++;
            }

            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], left, std::less<ot::TreeNode>()) - &mesh[begin]);
            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(left))))
            {
                num_boundary_faces++;
            }

            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], right, std::less<ot::TreeNode>()) - &mesh[begin]);
            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(right))))
            {
                num_boundary_faces++;
            }
            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], front, std::less<ot::TreeNode>()) - &mesh[begin]);
            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(front))))
            {
                num_boundary_faces++;
            }

            found_pt=(std::lower_bound(&mesh[begin], &mesh[end-1], back, std::less<ot::TreeNode>()) - &mesh[begin]);
            if(found_pt==0 || ((found_pt==(end-1-begin)) && (!mesh[end-1].isAncestor(back)))) {


                num_boundary_faces++;
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

//@author: Milinda Fernando.
// This function is to calculate the boundary faces with some flexibility allowed.
// Assume that the given octree vector is sorted.


double calculateOptimalBoundaries(const std::string& fileprefix, int q,double slack,int opt)
{

    int com_size;
    MPI_Comm_size(MPI_COMM_WORLD,&com_size);
    std::vector<ot::TreeNode> globalNodes;
    char ptsFileName[256];
//    for(int i=0;i<com_size;i++)
//    {
//        ot::readNodesFromFile()
//        ot::writePtsToFile()
//    }






    return 0.0;


}