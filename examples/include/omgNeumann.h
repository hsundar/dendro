#ifndef OMG_NEUMANN
#define OMG_NEUMANN

void solve_neumann(
    /* input parameters: */
     std::vector<double>& pts,
     void (*CalcCoef)(const std::vector<double> & pts, std::vector<double> & values),
     void (*CalcRHS)(const std::vector<double> & pts, std::vector<double> & values),
     int numMultigridLevels,
     /* output parameters */
     Vec& sol,
     ot::DAMG * & damg
    );

void getPrivateMatricesForKSP_Shell_Jac2(Mat mat, Mat *AmatPrivate, Mat *PmatPrivate, MatStructure* pFlag);

void solve_neumann_oct(
    /* input parameters: */
     std::vector<ot::TreeNode>& octs,
     void (*CalcCoef)(const std::vector<double> & pts, std::vector<double> & values),
     void (*CalcRHS)(const std::vector<double> & pts, std::vector<double> & values),
     int numMultigridLevels,
     /* output parameters */
     Vec& sol,
     ot::DAMG * & damg
    );

#endif

