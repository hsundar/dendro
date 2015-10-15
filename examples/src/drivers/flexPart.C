//@author: Milinda Fernando.
// School of Computing
// University of Utah

// Partitioning with flexibility to find out the optimal partition blocks w.r.t surface area.

#include "TreeNode.h"
#include "octreeStatistics.h"
#include "hcurvedata.h"

#define OPT_LOW_BLOCK 1
#define OPT_SUM_BLOCK 2


void optimalPartition(const std::vector<ot::TreeNode> global_mesh,double slack,int q);
void calculateBoundaryFaces_Serial(const std::vector<ot::TreeNode> & mesh, int q,double* stat);

void calculateBoundaryFaces_Serial(const std::vector<ot::TreeNode> & mesh, int q,double* stat) {


    ot::TreeNode first=mesh[0];
    ot::TreeNode last=mesh[mesh.size()-1];
    ot::TreeNode R (first.getDim(), first.getMaxDepth());
    int com_size=1;
    int mesh_nodes=mesh.size();
    stat[0]=0;  //min faces
    stat[1]=0;  //max faces
    stat[2]=0;  // mean faces


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

    unsigned long min_faces=1<<20;
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
            //std::cout<<"Boundary Calculation For Node :"<<temp<<std::endl;
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

//    unsigned long global_sum, global_max, global_min;
//
//    MPI_Reduce(&total_boundary_faces, &global_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&max_faces, &global_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&min_faces, &global_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

    stat[0]=min_faces;
    stat[1]=max_faces;
    stat[2]=total_boundary_faces;///(double)(q*com_size);

}




void optimalPartition(const std::vector<ot::TreeNode> global_mesh,double slack,int q)
{
    unsigned long int mesh_nodes=global_mesh.size();
    unsigned long int slack_cnt=(unsigned long int) mesh_nodes*slack;
    unsigned long int fixed_par_mesh_size=mesh_nodes/q;

    assert(fixed_par_mesh_size > slack_cnt);

    std::vector<ot::TreeNode>* fixed_par=new std::vector<ot::TreeNode>[q] ;
    int begin=0;
    int end=0;
    for(int i=0;i<q;i++)
    {
        fixed_par[i]=std::vector<ot::TreeNode>(fixed_par_mesh_size);

        begin=i*q;
        end=begin+fixed_par_mesh_size;

        if(end+fixed_par_mesh_size>mesh_nodes)
            end=mesh_nodes;

        fixed_par[i].insert(fixed_par[i].begin(),global_mesh.begin()+begin,global_mesh.begin()+end);

        std::cout<<"block:"<<i<<"Size:"<<fixed_par[i].size()<<std::endl;

    }

    std::vector<ot::TreeNode>* flex_par=new std::vector<ot::TreeNode>[q] ;

    for(int j=0;j<q;j++)
    {
        flex_par[j]=std::vector<ot::TreeNode>(fixed_par_mesh_size);
        flex_par[j]=fixed_par[j];
    }

    std::cout<<"fixed blocks assigned to flex blocks"<<std::endl;

    double block1[3];
    double block2[3];
    int global_min=1<<20;
    int opt_slack=0;
    std::cout<<"Main Slack Loop Stat"<<std::endl;
    for(int i=1;i<slack_cnt;i++)
    {
      for(int j=0;j<q;j++)
      {
          flex_par[j]=fixed_par[j];
          //std::cout<<"flex par:"<<j<<flex_par[j].size();
      }

        std::cout<<"flex_blocks Restored:"<<std::endl;
      for(int j=0;j<q;j++)
      {
          // case 1: leftmost block
          if(j==0) {
              for(int k=0;k<i;k++)
                   flex_par[j].push_back(flex_par[j+1][k]);
              flex_par[j + 1].erase(flex_par[j + 1].begin(), flex_par[j].begin() + i);

              calculateBoundaryFaces_Serial(flex_par[j],1,block1);
              calculateBoundaryFaces_Serial(flex_par[j+1],1,block2);

          }else {
              //case 2: right most block

              //flex_par[j-1].insert(flex_par[j-1].end(),flex_par[j].begin(),flex_par[j].begin()+i);
              for(int k=0;k<i;k++)
                  flex_par[j-1].push_back(flex_par[j][k]);
              flex_par[j].erase(flex_par[j].begin(), flex_par[j].begin() + i);

              calculateBoundaryFaces_Serial(flex_par[j - 1], 1, block1);
              calculateBoundaryFaces_Serial(flex_par[j], 1, block2);

//              if (global_min > std::min(block1[2], block2[2])) {
//                  global_min = std::min(block1[2], block2[2]);
//                  opt_slack = i;
//              }
          }
          std::cout << RED " Boundary Surfaces (min):"<<block1[0]<< NRM << std::endl;
          std::cout << RED " Boundary Surfaces (max):"<<block1[1]<< NRM << std::endl;
          std::cout << RED " Boundary Surfaces (mean):"<<block1[2]<< NRM << std::endl;

          std::cout << RED " Boundary Surfaces (min):"<<block2[0]<< NRM << std::endl;
          std::cout << RED " Boundary Surfaces (max):"<<block2[1]<< NRM << std::endl;
          std::cout << RED " Boundary Surfaces (mean):"<<block2[2]<< NRM << std::endl;

          if(global_min>std::min(block1[2],block2[2])) {
              global_min = std::min(block1[2], block2[2]);
              opt_slack=i;
          }

      }
        std::cout<<"Opt_Slac_count:"<<opt_slack<<"\t after "<<i<<" slack"<<std::endl;

    }

    std::cout<<"Opt_ Slack Count:"<<opt_slack<<std::endl;
    std::cout<<"Min Max Surface Count:"<<global_min<<std::endl;

}




int main(int argc, char **argv) {

    int com_size=atoi(argv[1]);
    int dim=3;
    std::vector<ot::TreeNode> globalNodes;
    _InitializeHcurve();
    std::vector<ot::TreeNode> bal_oct;
    for (int i = 0; i < com_size; i++) {
        char ptsFileName[256];
        sprintf(ptsFileName, "%s_%d", argv[2], i);
        ot::readNodesFromFile(ptsFileName,bal_oct);
        for(int i=0;i<bal_oct.size();i++) {
            globalNodes.push_back(bal_oct[i]);
            //std::cout<<"GlobalNodes:"<<globalNodes[i]<<std::endl;
        }
        bal_oct.clear();
    }
    //std::cout << "NUmber of pts read:" << bal_pts.size()<<std::endl;
    std::cout<<"Total Number of Nodes:"<<globalNodes.size()<<std::endl;
    std::sort(globalNodes.begin(),globalNodes.end());
    optimalPartition(globalNodes,0.001,8);





}