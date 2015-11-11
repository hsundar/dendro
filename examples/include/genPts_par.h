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


#ifndef GEN_GAUSS_H
#define GEN_GAUSS_H


#include <iostream>
#include "mpi.h"
#include <iostream>
#include <random>
#include <chrono>


void genGauss(const double& sd, const int numPts, int dim, char * filePrefix,MPI_Comm comm);
void genLogarithmicGauss(const double& sd, const int numPts, int dim, char * filePrefix,MPI_Comm comm);


#endif