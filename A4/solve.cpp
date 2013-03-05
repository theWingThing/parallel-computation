/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "Plotting.h"
using namespace std;

void repNorms(ofstream& logfile, double l2norm, double mx, double dt, int m,int n, int niter, int stats_freq);


// Reports statistics about the computation: the L2 Norm and the Infinity NOrm
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
// The Linf norm is simply the maximum (absolute) value over
// all points in the array

 double stats(double **E, int m, int n, double *_mx)
 {
     double mx = -1;
     double l2norm = 0;

     for (int j=1; j<=m+1; j++)
     {
         for (int i=1; i<=n+1; i++) 
         {
             l2norm += E[j][i]*E[j][i];
             double fe = fabs(E[j][i]);
             if (fe > mx)
             { 
                 mx = fe;
             }
         }
     }

     // In the parallel version, you must sum all the local contributoins
     // before dividing by (m+1)*(n+1)
     l2norm /= (double) ((m+1)*(n+1));
     l2norm = sqrt(l2norm);

     *_mx = mx;
     return l2norm;
 }

int solve(ofstream& logfile, double ***_E, double ***_E_prev, double **R, int m, int n, 
    int niters, double alpha, double dt, int plot_freq, Plotter *plotter, int stats_freq)
{

     // Simulated time is different from the integer timestep number
     double t = 0.0;

     double **E = *_E, **E_prev = *_E_prev;
     int niter;

     // We continue to sweep over the mesh until the simulation has reached
     // the desired simulation Time
     // This is different from the number of iterations
     for (niter = 0; niter < niters; niter++)
     {
      
#ifdef DEBUG
         double mx;
         double l2norm = stats(E_prev,m,n,&mx);
         repNorms(logfile,l2norm,mx,dt,m,n,niter, stats_freq);
         if (plot_freq)
         { 
             plotter->updatePlot(E,  niter, m+1, n+1, WAIT);
         }
    //    splot(E_prev,niter,m+1,n+1,WAIT);
#endif
       
         /* 
          * Copy data from boundary of the computational box to the
          * padding region, set up for differencing computational box's boundary
          *
          * These are physical boundary conditions, and are not to be confused
          * with ghost cells that we would use in an MPI implementation
          *
          * The reason why we copy boundary conditions is to avoid
          * computing single sided differences at the boundaries
          * which increase the running time of solve()
          *
          */
          
          for (int j=1; j<=m+1; j++) 
          {
             E_prev[j][0] = E_prev[j][2];
          }

          for (int j=1; j<=m+1; j++) 
          { 
              E_prev[j][n+2] = E_prev[j][n];
          }
       
          for (int i=1; i<=n+1; i++) 
          { 
              E_prev[0][i] = E_prev[2][i];
          }

          for (int i=1; i<=n+1; i++) 
          {
              E_prev[m+2][i] = E_prev[m][i];
          }
         
         // Solve for the excitation, a PDE
          for (int j=1; j<=m+1; j++)
          {
              for (int i=1; i<=n+1; i++) 
              { 
                  E[j][i] = E_prev[j][i] + alpha * 
                      ( E_prev[j][i+1]+E_prev[j][i-1] - 
                        4*E_prev[j][i]+E_prev[j+1][i] +
                        E_prev[j-1][i] );
              }
          }

         /* 
          * Solve the ODE, advancing excitation and recovery variables
          *     to the next timtestep
          */
         for (int j=1; j<=m+1; j++)
         {
             double *RR = &R[j][1];
             double *EE = &E[j][1];

             for (int i=1; i<=n+1; i++, EE++,RR++)
             { 
                 EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
             }
         }

         for (int j=1; j<=m+1; j++)
         {
             double *RR = &R[j][1];
             double *EE = &E[j][1];

             for (int i=1; i<=n+1; i++, EE++, RR++)
             { 
                 RR[0] += dt * ( epsilon + 
                                 M1 * RR[0] / (EE[0]+M2)) * (-RR[0] - kk * EE[0] * (EE[0] - b - 1)
                                );
             }
         }

         if (stats_freq)
         {
             double mx;
             double l2norm = stats(E_prev,m,n,&mx);
             repNorms(logfile,l2norm,mx,dt,m,n,niter, stats_freq);
         }

         if (plot_freq)
         {
             if (!(niter % plot_freq))
             {
               //splot(E,niter,m+1,n+1,WAIT); 
               plotter->updatePlot(E,  niter, m+1, n+1, WAIT);
             }
          }

         // Swap current and previous
         double **tmp = E; E = E_prev; E_prev = tmp;
     }

      // Store them into the pointers passed in
      *_E = E;
      *_E_prev = E_prev;

      return niter;
}
