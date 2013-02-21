#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include "particle.h"
#include "common.h"
#include "Plotting.h"
#include "Bins.h"
#ifdef  _OPENMP
#include "omp.h"
#endif

using namespace std;

extern double dt;
extern double size;

//
//  Tuned constants
//
double density = 0.0005;
double  mass =    0.01;
double cutoff =  0.01;
double min_r  =  (cutoff/100);

int cmdLine( int argc, char **argv, int&n, int& nt, int& chunk, int& nsteps, int& ntlot, int& sd, bool& imbal, char** savename, int& x, int& y);
void SimulateParticles(int nsteps, particle_t *particles, int n, int nt, int chunk, int nplot, bool imbal, double &uMax, double &vMax, double &uL2, double &vL2, Plotter *plotter, FILE *fsave, int nx, int ny );
void RepNorms(double uMax,double vMax,double uL2,double vL2);

int main( int argc, char **argv )
{    
	int nt;     // number of threads
    int n;      // # of particles
    int nsteps; // # of timesteps
    int nplot;  // Plotting frequency
    int sd;     // Random seed
    char *savename;     // Name of file to save output
    bool imbal; // Use an imbalanced initial particle distribition
    int chunk;  // Chunk size for dynamic scheduling
                // When chunk=0, we use static BLOCK scheduling
	int nx, ny; // number of columns and rows of boxes

	int OK =  cmdLine( argc, argv, n, nt, chunk, nsteps, nplot, sd, imbal, &savename, nx, ny);
	
    // If there was a parsing error, exit
    if (!OK){
        exit(-1);
    }

    set_size( n );

    //sanity check
    if( size / (double) nx < cutoff )
    {
	    cout << "nx is too large" << endl;
	    return 0;
    }

    cout << endl << "dt: " << dt << endl;
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    particle_t *particles = new particle_t[ n ];
    assert(particles);
    
    init_particles( n, particles,sd, imbal);
    double uMax, vMax, uL2, vL2;
    Plotter *plotter = NULL;
    if (nplot){
        plotter = new Plotter();
        assert(plotter);
        VelNorms(particles,n,uMax, vMax, uL2, vL2);
        plotter->updatePlot(particles,n,0,uMax,vMax,uL2,vL2);
    }
    cout << "# particles : " << n << endl;
    if (chunk == -1)
       cout << "Using static BLOCK decomposition" << endl;
    else
       cout << "Using DYNAMIC scheduling with a chunk size of " << chunk << endl;
    if (imbal)
       cout << " Irregular particle distribution" << endl;
    cout << "Nsteps: " << nsteps << endl;
    cout << "Partition: " << nx << " x "<< ny << endl;
#ifdef  _OPENMP
   printf("Compiled with openmp enabled\n");
   omp_set_num_threads(nt);
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if (tid == 0) {
        int nthreads = omp_get_num_threads();
        cout << "Number of openMP threads:  " << nthreads << endl;
    }
  }
#ifdef DYN
    cout << "Compiled with dynamic scheduling.";
#ifdef CHUNK
    cout << "  Chunksize = " << CHUNK << endl;
#else
    cout << endl;
#endif
#endif
#else
    cout << "Number of C++11 threads:  " << nt << endl;
#endif

    // Box the particles into nx by ny regions
    Bins bins(size, nx, ny, nt, n, particles);
  
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    bins.SimulateParticles(nsteps,particles,n,nt,chunk, nplot, imbal, uMax, vMax, uL2, vL2, plotter, fsave, nx, ny, dt );
    simulation_time = read_timer( ) - simulation_time;
   
    cout << endl;
    cout <<  "n = " << n << ", nsteps = " << nsteps << endl;
    VelNorms(particles,n,uMax,vMax,uL2,vL2);

    RepNorms(uMax,vMax,uL2,vL2);
    cout <<  "Running time = " << simulation_time << " sec.\n";

    if( fsave )
        fclose( fsave );
    
    if (nplot)
        delete plotter;
    delete [ ] particles;
    
    return 0;
}
