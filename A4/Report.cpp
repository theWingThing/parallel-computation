// 
// Performs various reporting functions
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#ifdef _MPI_
#include "mpi.h"
#endif
using namespace std;

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem

void ABEND()
{
   cout.flush();
   cerr.flush();
#ifdef _MPI_
   MPI_Abort(MPI_COMM_WORLD,-1);
#else
   exit(-1);
#endif
}
 
void Stop(){
   cout.flush();
   cerr.flush();
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(-1);
}

// Report statistics periodically
void repNorms(ofstream& logfile, double l2norm, double mx,  double dt, int m,int n, int niter, int stats_freq){

     int myrank;
#ifdef _MPI_
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
     myrank = 0;
#endif
     if (!myrank){
        if ( !(niter % stats_freq)) {
            cout <<      setw(6);
            cout.setf(ios::fixed);
            cout << "iteration= " << niter << ", t= ";
            cout.unsetf(ios::fixed);
            cout.setf(ios::scientific);
            cout.precision(6);
            cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
            logfile <<          setw(6);
            logfile.setf(ios::fixed);
            logfile << "iteration= " << niter << ", t= ";
            logfile.unsetf(ios::fixed);
            logfile.setf(ios::scientific);
            logfile.precision(6);
            logfile << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
        }
    }
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
void printTOD(ofstream& logfile,string mesg)
{
     time_t tim = time(NULL);
     string s = ctime(&tim);
     int myrank;
#ifdef _MPI_
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
     myrank = 0;
#endif
     if (!myrank){
        cout << endl;
        logfile << endl;
        if (mesg.length() ==  0) {
            cout << "Time of day: " << s.substr(0,s.length()-1) << endl;
            logfile << "Time of day: " << s.substr(0,s.length()-1) << endl;
        }
        else {
            cout << "[" << mesg << "] " ;
            cout << s.substr(0,s.length()-1) << endl;
            logfile << "[" << mesg << "] " ;
            logfile << s.substr(0,s.length()-1) << endl;
        }
        cout << endl;
        logfile << endl;
    }
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


// Computes the gigaflops rate

double gflops(int n, int niter, double time){
    int n2 = n*n;
    int64_t updates = (int64_t) n2 * (int64_t) niter;
    int64_t flops = 28 * updates;
    double flop_rate = (double) flops / time;
    return ( flop_rate/1.0e9);
}


void ReportEnd(ofstream& logfile, int niters, double l2norm, double mx, int m,int n, double t0, int px, int py){
    printTOD(logfile,"Simulation completes");    
    int myrank;

#ifdef _MPI_
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
    myrank = 0;
#endif
    if (!myrank){
        double gf = gflops(n+1, niters, t0);
	cout << "End at";
        cout <<          setw(6);
        cout.setf(ios::fixed);
        cout << " iteration " << niters << endl;
        cout.unsetf(ios::fixed);
        cout.setf(ios::scientific);
        cout.precision(6);
        cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
        cout.unsetf(ios::scientific);
        cout.unsetf(ios::fixed);
        cout.precision(6);
        cout << "Running Time: " << t0 << " sec." << endl;
        cout.precision(3);
        cout << "GFlop rate: " << gf << endl << endl;

        cout << "  M x N, [px x py]       [#iter] {T_p, Gflops} Linf, L2" << endl;
        cout << "># " << m  << " " << n << " ";
        cout.precision(3);
        cout << "[" << px << " " << py << "] ";
        cout.precision(6);
        cout <<  "\t[ " << niters << "] ";
        cout.precision(4);
//        cout << t0 << " "  << gf << " ";
        cout << "{" << t0 << " "  << gf << "} ";

        cout.unsetf(ios::fixed);
        cout.setf(ios::scientific);
        cout.precision(5);
        cout << mx << " " << l2norm << endl;

	logfile << "End at";
        logfile <<          setw(6);
        logfile.setf(ios::fixed);
        logfile << " iteration " << niters << endl;
        logfile.unsetf(ios::fixed);
        logfile.setf(ios::scientific);
        logfile.precision(5);
        logfile << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
        logfile.unsetf(ios::scientific);
        logfile.unsetf(ios::fixed);
        logfile.precision(6);
        logfile << "Running Time: " << t0 << " sec." << endl;
        logfile << "{" << t0 << " "  << gf << "} ";
        logfile.precision(3);
        logfile << "GFlop rate: " << gf << endl << endl;
        logfile.precision(5);

        logfile << "  M x N, [px x py, NT] [#iter] {T_p, Gflops} Linf, L2" << endl;
        logfile << "># " << m  << " " << n << " ";
        logfile.precision(3);
        logfile << "[" << px << " " << py << "] ";
        logfile.precision(6);
        logfile <<  "\t[ " << niters << "] ";
        logfile.precision(4);
        logfile << "{t0" << " "  << gf << "} ";

        logfile.unsetf(ios::fixed);
        logfile.setf(ios::scientific);
        logfile.precision(5);
        logfile << mx << " " << l2norm << endl;
    }

#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void ReportStart(ofstream& logfile,double dt, int niters, int m, int n, int px, int py, bool noComm){
  
    int myrank;
#ifdef _MPI_
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
     myrank = 0;
#endif
     if (!myrank){
        cout << "dt= " << dt << ", ";
	cout << "# iters = " << niters << endl;
        cout << "m x n = " << m << " x " << n << endl;
        cout << "processor geometry: " << px << " x " << py << endl;
        cout << endl;

        logfile << "dt= " << dt << ", ";
	logfile << "# iters = " << niters << endl;
        logfile << "m x n = " << m << " x " << n << endl;
        logfile << "processor geometry: " << px << " x " << py << endl;
        logfile << endl;


#ifdef _MPI_
        cout << "Compiled with MPI ENabled\n";
        logfile << "Compiled with MPI ENabled\n";
        if (noComm){
            cout << "Communication shut off" << endl;
            logfile << "Communication shut off" << endl;
        }
#else
        cout << "Compiled with MPI DISabled\n";
        logfile << "Compiled with MPI DISabled\n";
#endif

     }
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
