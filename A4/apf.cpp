/* 
 * Driver for cardiac electrophysiology simulation using Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <string>
#include <math.h>
#include "apf.h"
#include "Plotting.h"
#ifdef _MPI_
#include <mpi.h>
#endif
using namespace std;


// External functions
void cmdLine(int argc, char *argv[], int& n, int& stats_freq, int& plot_freq, int& px, int& py, bool &noComm, int &niters);
int solve(ofstream& logfile, double ***_E, double ***_E_prev, double **R, int m, int n, int niters, double alpha, double dt, int plot_freq, Plotter *plotter, int stats_freq);
void printTOD(ofstream& logfile, string mesg);
double stats(double **E, int m, int n, double *_mx);
void ReportStart(ofstream& logfile, double dt, int niters, int m, int n, int px, int py, bool noComm);
void ReportEnd(ofstream& logfile, int niters, double l2norm, double mx, int m,int n, double t0, int px, int py);
double **alloc2D(int m,int n);
double ComputeDt(int n, double& alpha);
void init(double **E,double **E_prev,double **R, int m0, int n0, int m,int n, int global_m, int global_n);
double getTime();

// Main program
int main(int argc, char** argv)
{
     /*
      *  Solution arrays
      *   E is the "Excitation" variable, a voltage
      *   R is the "Recovery" variable
      *   E_prev is the Excitation variable for the previous timestep,
      *      and is used in time integration
      */
     double **E, **R, **E_prev;

     // Default values for the command line arguments
     int m=100,n=100;
     int stats_freq = 0;
     int plot_freq = 0;
     int px = 1, py = 1;
     int niters=100;
     bool noComm = false;

#ifdef _MPI_
     MPI_Init(&argc,&argv);

#endif

    // Parse command line arguments
     cmdLine( argc, argv, n, stats_freq,  plot_freq, px, py, noComm, niters);

     if (n < 26)
     {
        cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
        exit(-1);
     }

     m = n;

#ifdef _MPI_
     int nprocs=1, myrank=0;
     int local_n = n, local_m = m; 
         
     MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank > 0)
     {
         // fix for x-axis for now
         local_m /= nprocs;

         if(myrank == nprocs -1) 
         {
             local_m += (m % nprocs); 
         }
     }

     ofstream logfile("Log.txt",ios::out);
     printTOD(logfile, "Simulation begins");

     // for each process we create array for computation size
     // n over p by n in 1D except process 0 which has a whole array
     E = alloc2D(local_m+3,local_n+3);
     E_prev = alloc2D(local_m+3,local_n+3);
     R = alloc2D(local_m+3,local_n+3);

	 //cout << "alloc array to " << myrank << endl;

     // now set m and n to be the size of sub array
     // for each task assigned to a process
     // m = x axis, n = y axis
     
#else
     // The log file
     // Do not change the file name or remove this call
     ofstream logfile("Log.txt",ios::out);
     printTOD(logfile, "Simulation begins");

     // Allocate contiguous memory for solution arrays
     // The computational box is defined on [1:m+1,1:n+1]
     // We pad the arrays in order to facilitate differencing on the 
     // boundaries of the computation box
     E = alloc2D(m+3,n+3);
     E_prev = alloc2D(m+3,n+3);
     R = alloc2D(m+3,n+3);
#endif

     // To initialize the meshes, a process needs to know the coordinates of
     // the local  partition with respect to the global coordinate system
     // We need this information because initial values depend on  absolute
     // (global) coordinates
     // When we are running with 1 process, the mapping is trivial:
     // the local and global coordinates are the same
     // But with more than 1 process, we need global coordinates (m0,n0)
     // corresponding to the local coordinates (0,0)
     // We also need to know the global bounds of the array (m_global, n_global)
     //
     // For example, if N=999 (m_global=n_global=999) and px=2, py=2
     // (m0,n0) should be (0,0) for the process of rank 0
     // (500,0) for the process of rank 1.
     // For the process of rank 2, (m0,n0)= (0,500), for rank 3, (m0,y0)=(500,500)
     //
     // We provide (m0,n0) for 1 process. You must write code to compute
     // the global coordinate
     
     int m0 = 0, n0 = 0;

     int m_global=m, n_global=n;
#ifdef _MPI_

     // in 1D we set m to n in n by n global array 
     // m0 = 0 for all process
     // m is x-axis
     m0 = (n/nprocs) * myrank;
#endif
    if(myrank == 0)
    {
        init(E, E_prev,R, m0, n0, m, n, m_global, n_global);
#if 0
         cout << "rank--->"<<myrank<<endl;
         for(int i = 0; i < local_m + 3; i++)
         {
             cout << i << ", ";
         }
         cout << endl;

         for(int i = 0; i < m + 3; i++)
         {
             cout << "row: " << i << endl;
             for(int j = 0; j < n +3; j++)
             {
                 cout << R[i][j] <<", ";
             }
             cout <<endl;
         }
#endif
    }


    if(myrank > 0)
    {
         init(E, E_prev, R, m0, n0, local_m, local_n, m_global,n_global);
#if 0
         cout << "rank--->"<<myrank<<endl;
         for(int i = 0; i < local_m + 3; i++)
         {
             cout << i << ", ";
         }
         cout << endl;

         for(int i = 0; i < local_m + 3; i++)
         {
             cout << "row: " << i << endl;
             for(int j = 0; j < local_n +3; j++)
             {
                 cout << R[i][j] <<", ";
             }
             cout <<endl;
         }
#endif
    }
     //
     // Initialize two simulation parameters: timestep and alpha
     // Do not remove this call or the code will not run correctly
     //
     double alpha;
     double dt = ComputeDt(n,alpha);
     // Report various information
     // Do not remove this call, it is needed for grading
     ReportStart(logfile, dt, niters, m, n, px, py, noComm);

     Plotter *plotter = NULL;

     if(myrank==0)
     {
         if (plot_freq)
         {
             plotter = new Plotter();
             assert(plotter);
         }
     }
     // Start the timer
#ifdef _MPI_
     //cout << "barrier #1 start" << endl;
     MPI_Barrier(MPI_COMM_WORLD);
     //cout << "barrier #1 end" << endl;
     double local_t1 = -MPI_Wtime();
     double t0;
#else
     double t0 = -getTime();
     int niter = solve(logfile, &E, &E_prev, R, m, n, niters, alpha, dt, plot_freq, plotter, stats_freq);
#endif

#ifdef _MPI_
     //cout << "solve #1 start" << endl;
     //cout << "myrank: "<< myrank << endl;
     //cout << "m: "<< m << " n: " << n << endl;
    int niter = solve(logfile, &E, &E_prev, R, m, n, niters, alpha, dt, plot_freq, plotter, stats_freq);
     
     // find the max from each process
     local_t1 += MPI_Wtime();

     MPI_Reduce(&local_t1, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else
     t0 += getTime();
#endif

     if (niter != niters)
        cout << "*** niters should be equal to niters" << endl;
      // Report various information
      // Do not remove this call, it is needed for grading
      double mx;
      double l2norm = 0; 

      if(myrank == 0)
      {
          l2norm = stats(E_prev,m,n,&mx);
          ReportEnd(logfile,niters,l2norm,mx, m,n,t0,px, py);
      }
      else 
      { 
          ReportEnd(logfile,niters,l2norm,mx,local_m,local_n,t0,px, py);
      }

#ifdef _MPI_
     if(myrank == 0)
     {
         if (plot_freq)
         {
            cout << "\n\nEnter any input to close the program and the plot...";
            int resp;
            cin >> resp;
          }
    }
#endif
     logfile.close();

     free (E);
     free (E_prev);
     free (R);

#ifdef _MPI_
     if(myrank == 0)
     {
#endif
         if (plot_freq)
             delete plotter;
#ifdef _MPI_
     }

     MPI_Finalize();
#endif
}
