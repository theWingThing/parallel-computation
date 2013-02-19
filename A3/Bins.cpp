//
// A partition is a 2D array of boxes, which store the particles
//
#include <iostream>
#include "Bins.h"
#include "Box.h"

Bins::Bins(double size, int nx, int ny, int np, int n, particle_t* particles) : _size(size), _nx(nx), _ny(ny), _np(np), _n(n)
{
	boxesCount = nx * ny;
	boxes = new Box[boxesCount];
	boxWidth = size / (double) nx;
	boxHeight = size / (double) ny;
	for(int i=0; i < nx; ++i) {
		for(int j=0; j < ny; ++j) {			
			int binID = j*nx + i;	      
			boxes[binID].I = i;
			boxes[binID].J = j;
			boxes[binID].bins = this;
		}
	}
	SortParticles(n, particles);
}

//
// Sort only once at the start
//
void Bins::SortParticles(int n, particle_t* particles)
{
	for(int p=0; p < n; ++p) {
		int i = (int) (particles[p].x / boxWidth);
		int j  = (int) (particles[p].y / boxHeight);
		int boxID = j* _nx + i;

		boxes[boxID].AddParticle(&particles[p]);
	}
}

void call_apply_forces(int i, int nt, int boxesCount, Box* boxes){
  for(int j = i * (boxesCount/nt); j < (i+1) * (boxesCount/nt) - 1; j++)
    boxes[j].apply_forces();
} 

//
// Apply forces to all the bins
//
void Bins::apply_forces(int nt)
{
/*#ifdef  _OPENMP
#ifdef DYN
#if CHUNK
#pragma omp parallel for  schedule(dynamic,CHUNK)
#else
#pragma omp parallel for  schedule(dynamic)
#endif
#else
#pragma omp parallel for
#endif
#endif
  for(int i=0; i < boxesCount; ++i)
    boxes[i].apply_forces();
*/ 
  thread *thrds = new thread[nt];
  for(int i = 0; i < nt; ++i)
      thrds[i] = thread(call_apply_forces, i, nt, boxesCount, boxes);
  for(int i = 0; i < nt; ++i)
    thrds[i].join();

}

void call_move_particles(int i, int nt, int boxesCount, Box* boxes, double dt){
  for(int j = i * (boxesCount/nt); j < (i+1) * (boxesCount/nt) - 1; j++)
    boxes[j].move_particles(dt);
} 


//
// Move the particles in all the bins
//
void Bins::move_particles(double dt, int nt)
{
/*#ifdef  _OPENMP
#ifdef DYN
#if CHUNK
#pragma omp parallel for shared(dt) schedule(dynamic,CHUNK)
#else
#pragma omp parallel for shared(dt) schedule(dynamic)
#endif
#else
#pragma omp parallel for shared(dt)
#endif
#endif
  for(int i=0; i < boxesCount; ++i)
    boxes[i].move_particles(dt);
*/ 
  thread *thrds = new thread[nt];
  for(int i = 0; i < nt; ++i)
      thrds[i] = thread(call_move_particles, i, nt, boxesCount, boxes, dt);
  for(int i = 0; i < nt; ++i)
    thrds[i].join();


//
// After moving the particles, we check each particle
// to see if it moved outside its curren box
// If, so we append to the inbound partcle list for the new bin
// As written, this code is not threadsafe
//
  for(int i=0; i < boxesCount; ++i)
    boxes[i].UpdateParticlesBox();


//
// After we've updated all the inbound lists
// We then make a new pass, incorporating into this bin,
// any particles on the inbound list
// This work parallelizes
//

#ifdef  _OPENMP
#ifdef DYN
#if CHUNK
#pragma omp parallel for schedule(dynamic,CHUNK)
#else
#pragma omp parallel for schedule(dynamic)
#endif
#else
#pragma omp parallel for 
#endif
#endif
  for(int i=0; i < boxesCount; ++i)
    boxes[i].UpdateInboundParticles();
  

#ifdef _DEBUG
//
// May come in handy during debugging
//
  int particleCount = 0;
  for(int i=0; i < boxesCount; ++i) {
      particleCount += (int) boxes[i].boxParticles.size();
      for( unsigned int x=0; x <  boxes[i].boxParticles.size(); ++x) {
	particle_t* p = boxes[i].boxParticles[x];
	for( int j=0; j < boxesCount; ++j) {
	  if( i == j)
	    continue;
	  for(  unsigned int y=0; y < boxes[j].boxParticles.size(); ++y) {
	    if(p == boxes[j].boxParticles[y])
	      cout << "same particle detected in different box\n" ;
	  }
	}
      }
    }
  if(particleCount != _n)
    cout << "particle count = " << particleCount << endl;
#endif
}

void Bins::SimulateParticles(int nsteps, particle_t* particles, int n, int nt, int chunk, int nplot, bool imbal, double &uMax, double &vMax, double &uL2, double &vL2, Plotter *plotter, FILE *fsave, int nx, int ny, double dt ){
    for( int step = 0; step < nsteps; step++ ) {
    //
    //  compute forces
    //
	apply_forces(nt);

//     Debugging output
//      list_particles(particles,n);
    
    //
    //  move particles
    //
	move_particles(dt, nt);


	if (nplot && ((step % nplot ) == 0)){

	// Computes the absolute maximum velocity 
	    VelNorms(particles,n,uMax,vMax,uL2,vL2);
	    plotter->updatePlot(particles,n,step,uMax,vMax,uL2,vL2);
	}

//
// Might come in handy when debugging
// prints out summary statistics every time step
//
	//VelNorms(particles,n,uMax,vMax,uL2,vL2);
	
    //
    //  if we asked, save to a file
    //
	if( fsave && (step%SAVEFREQ) == 0 )
	    save( fsave, n, particles );
    }
}


Bins::~Bins()
{
	delete [] boxes;
}
