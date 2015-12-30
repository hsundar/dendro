/*
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * milinda@cs.utah.edu

 * Functions needed for Flexible partition calculation, Where we trade off load balance for achieve lower communication costs.
 *
 * */

#include "TreeNode.h"
#include "parUtils.h"
#include <vector>


#define BAL_CONST 4

/*
 * Input:
 * @parameter partiton: Partition we need to compute the boundary (partition should be unique and sorted and balanced)
 *
 * Output:
 * @parameter boundary: list of boundary octants (sorted)
 * @parameter bdy_surfaces: Number of boundary surfaces.
 * @parameter bdy_surface_data: How many boundary surfaces for each octant.
 *
 * */
long BoundaryCalculation(std::vector<ot::TreeNode>& partition,int begin,int end, char* bdy_data);


/*
 * Input:
 * @parameter new_octant: List of new octants
 *
 * Output:
 * @parameter boundary: Updated boundary octants based on the new_octants
 * @parameter bdy_surfaces: Updated Number of boundary surfaces.
 *
 * */


void UpdateBoundaryCalculation(std::vector<ot::TreeNode>& partition,int begin, int end,long& bdy_surfaces,char* bdy_data);


/*
 *Input: @parameter: partition: Original partitions, slack: should be [0,1]
 * Output: @parameter: New partition,
 *
 * */
void DynamicPartitioning(std::vector<ot::TreeNode>& partition, double slack, MPI_Comm comm);











