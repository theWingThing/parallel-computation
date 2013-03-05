// Process command line arguments
// 
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#ifdef _MPI_
#include <mpi.h>
#endif
using namespace std;

void Stop();

void cmdLine(int argc, char *argv[], int& n, int& stats_freq, int& plot_freq, int& px, int& py, bool &noComm, int &niter){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"stats-freq", required_argument, 0, 's'},
        {"plot", required_argument, 0, 'p'},
	{"px", required_argument, 0, 'x'},
	{"py", required_argument, 0, 'y'},
	{"niter", required_argument, 0, 'i'},
	{"nocomm", no_argument, 0, 'k'},

 };
    int nprocs=1, myrank=0;
#ifdef _MPI_
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif

    // Process command line arguments
 for(int ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:x:y:i:j:s:kp:",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                n = atoi(optarg);
                break;
	   //
           // X processor geometry
            case 'x':
                px = atoi(optarg);
                break;

            // X processor geometry
            case 'y':
                py = atoi(optarg);
                break;

            // # of iterations
	    // Use this option control the number of mesh sweeps
            case 'i':
                niter = atoi(optarg);
                break;


	    // Print statistics for assessing correctness
            case 's':
                stats_freq = atoi(optarg);
                break;


	    // Plot the excitation variable
            case 'p':
                plot_freq = atoi(optarg);
                break;

            // Shut off communication
            case 'k':
                noComm = true;
                break;

	    // Error
            default:
                cout << "Usage: apf [-n <domain size>] [-i <# iterations>]";
                cout << "\n\t    ";
                cout << " [-s <stats frequency>[-p <plot frequency>]\n\t";
		cout << "     [-x <x processor geometry>]";
		cout << "     [-y <x processor geometry>]";
		cout << "     [-k <no communication>]" << endl;
                cout << endl;
                exit(-1);
            }
    }
 }
 if ((px * py) != nprocs){
    if (!myrank)
        cout << "\n *** The number of processes in the geometry (" << px*py << ")  is not the same as the number requested (" << nprocs << ")" << endl << endl;
    Stop();
 }
 if (!myrank)
    if ((px * py) > 1){
        cout << "\n *** The number of processes in the geometry > 1, but you have not enabled MPI\n";
        Stop();
    }
}
