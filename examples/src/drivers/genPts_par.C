/*
 *
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @date: 11/10/2015
 *
 *
 * Contains the Normal(Gaussian) and Logarithmic normal random number (octants) generator based on new c+11 rand engines.
 *
 *
 * */

#include "genPts_par.h"


void genGauss(const double& sd, const int numPts, int dim, char * filePrefix,MPI_Comm comm)
{


    int rank,size;

    char ptsFileName[256];

    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);


    std::random_device rd;
    std::mt19937 gen(rd());

    sprintf(ptsFileName, "%s%d_%d.pts", filePrefix, rank, size);

    std::normal_distribution<> d(0,sd);
    double* xyz = new double[dim*numPts];
    double temp=0;
    for(long i=0;i<(long)dim*numPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }

    FILE* outfile;
    outfile = fopen(ptsFileName,"wb");
    fwrite(&numPts, sizeof(unsigned int), 1, outfile);
    fwrite(xyz, sizeof(double), (dim*numPts), outfile);
    fclose(outfile);

    delete [] xyz;


}


void genLogarithmicGauss(const double& sd, const int numPts, int dim, char * filePrefix,MPI_Comm comm)
{


    int rank,size;

    char ptsFileName[256];

    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);


    std::random_device rd;
    std::mt19937 gen(rd());

    sprintf(ptsFileName, "%s%d_%d.pts", filePrefix, rank, size);

    std::lognormal_distribution<>d(0,sd);

    
    double* xyz = new double[dim*numPts];
    double temp=0;
    for(long i=0;i<(long)dim*numPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }

    FILE* outfile;
    outfile = fopen(ptsFileName,"wb");
    fwrite(&numPts, sizeof(unsigned int), 1, outfile);
    fwrite(xyz, sizeof(double), (dim*numPts), outfile);
    fclose(outfile);

    delete [] xyz;



}
