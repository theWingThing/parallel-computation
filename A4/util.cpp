// Utilities for cardiac electrophysiology simulator

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iomanip>
#include "apf.h"
using namespace std;

// Allocate a 2D array

double **alloc2D(int m,int n){

   double **E;
   int nx=n, ny=m;
   E = (double**)malloc(sizeof(double*)*ny + sizeof(double)*nx*ny);
   assert(E);
   for(int j=0;j<ny;j++) E[j] = (double*)(E+ny) + j*nx;
   return(E);
}

void init (double **E,double **E_prev,double **R, int m0, int n0, int m,int n, int m_global, int n_global){
    int i,j;
    // Initialization
    for (j=1; j<=m + 1; j++)
        for (i=1; i<= n+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }

    for (j=1; j<=m + 1; j++)
    {
        for (i=n/2+2; i<= n+1 ; i++){
            int i1 = i + n0;
            if (i1 >= n_global/2+1)
            {
                E_prev[j][i] = 1.0;
            }

        }

    }

    for (j=1; j<=m+1; j++){
        int j1 = j+m0 ;
        if (j1 >=m_global/2+1)
        for (i=1; i<=n+1; i++)
            R[j][i] = 1.0;
    }
}

 //
 // Compute the timestep dt
double ComputeDt(int n, double& alpha){
// We should always use double precision values for the folowing variables:
//    rp, dte, dtr, ddt
//
// This ensures that the computation of dte and especially dt
// will not lose precision (i.e. if computed as single precision values)
// 
 const double dx = 1.0/(double)n;
 const double rp= kk*(b+1)*(b+1)/4;
 const double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
 const double dtr=1/(epsilon+((M1/M2)*rp));
 const double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;

 double dt = (double) ddt;
 alpha = d*dt/(dx*dx);
 return dt;
 }

