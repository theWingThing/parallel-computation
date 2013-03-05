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
void init (double **E,double **E_prev,double **R,int m,int n);
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
 if (n < 26){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    exit(-1);
 }
 m = n;
 int nprocs=1, myrank=0;
#ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif

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

 init(E,E_prev,R,m,n);

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
 if (plot_freq){
     plotter = new Plotter();
     assert(plotter);
 }

 // Start the timer
#ifdef _MPI_
 double t0 = -MPI_Wtime();
#else
 double t0 = -getTime();
#endif
 int niter = solve(logfile, &E, &E_prev, R, m, n, niters, alpha, dt, plot_freq, plotter, stats_freq);

#ifdef _MPI_
 t0 += MPI_Wtime();
#else
 t0 += getTime();
#endif

if (niter != niters)
   cout << "*** niters should be equal to niters" << endl;
 // Report various information
 // Do not remove this call, it is needed for grading
 double mx;
 double l2norm = stats(E_prev,m,n,&mx);
 ReportEnd(logfile,niters,l2norm,mx,m,n,t0,px, py);

 if (plot_freq){
    cout << "\n\nEnter any input to close the program and the plot...";
    int resp;
    cin >> resp;
  }

 logfile.close();
 free (E);
 free (E_prev);
 free (R);
 if (plot_freq)
     delete plotter;
#ifdef _MPI_
 MPI_Finalize();
#endif
}
