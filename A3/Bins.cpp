//
// A partition is a 2D array of boxes, which store the particles
//
#include <iostream>
#include "Bins.h"
#include "Box.h"
#include <mutex>
#include <condition_variable>
mutex m;
condition_variable cond_var;

int count = 0;
int num_threads = 1;
int chunk_size = 0;

struct ThreadArgs{
    int tid;
    int queue_size;
    int* job_queue;
    int special;
} *args;

class SelfScheduler
{
    private:
        int problem_size;
        int chunk_size;
        int num_processors;
        int counter;
        mutex critical_section;

        SelfScheduler(SelfScheduler& s){};

    public:
        SelfScheduler():problem_size(0), chunk_size(0),
            num_processors(0), counter(0){};

        SelfScheduler(int n,int P, int chunk):problem_size(n),
            chunk_size(chunk), num_processors(P), counter(0){};

        ~SelfScheduler(){};
        bool isDone()
        {
            return counter > problem_size;
        }
        bool getChunk(int& min, int& max)
        {
            critical_section.lock();
            int k;
            k = counter;
            counter += chunk_size;
            critical_section.unlock();

            if(counter > (problem_size))
            {
                return false;
            }

            min = k;
            max = k + chunk_size;

            return true;
        }
        void reset()
        {
            counter = 0;
        }
};

SelfScheduler *app_forces, *upd_part, *move_part;

void barrier(SelfScheduler* sche)
{
    std::unique_lock<std::mutex> lock(m);
    count++;
    if(count == num_threads || sche->isDone())
    {
        count = 0;
        cond_var.notify_all();
    }
    else
    {
        cond_var.wait(lock);
    }
    lock.unlock();

}

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

void dynamic_apply_forces(int id, Box* boxes){
  for(int k = 0; k < args[id].queue_size; k++){
    boxes[args[id].job_queue[k]].apply_forces();
  }
}
//
// Apply forces to all the bins
//
void Bins::apply_forces(int nt)
{/*
#ifdef  _OPENMP
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
  if(!chunk_size){
      for(int i = 0; i < nt; ++i){
          thrds[i] = thread(call_apply_forces, i, nt, boxesCount, boxes);
      }
  }
  else{
      for(int i = 0; i < nt; ++i){
          thrds[i] = thread(dynamic_apply_forces, i, boxes);
      }
  }
  for(int i = 0; i < nt; ++i)
    thrds[i].join();
  delete [] thrds;
}

void call_move_particles(int i, int nt, int boxesCount, Box* boxes, double dt){
  for(int j = i * (boxesCount/nt); j < (i+1) * (boxesCount/nt) - 1; j++)
    boxes[j].move_particles(dt);
} 

void dynamic_move_particles(int id, Box* boxes, double dt){
  for(int k = 0; k < args[id].queue_size; k++){
    boxes[args[id].job_queue[k]].move_particles(dt);
  }
}

void call_update_particles(int i, int nt, int boxesCount, Box* boxes){
  for(int j = i * (boxesCount/nt); j < (i+1) * (boxesCount/nt) - 1; j++)
    boxes[j].UpdateInboundParticles();
}
void dynamic_update_particles(int id, Box* boxes){
  for(int k = 0; k < args[id].queue_size; k++){
    boxes[args[id].job_queue[k]].UpdateInboundParticles();
  }
}

//
// Move the particles in all the bins
//
void Bins::move_particles(double dt, int nt)
{/*
#ifdef  _OPENMP
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
   if(!chunk_size){
  for(int i = 0; i < nt; ++i)
      thrds[i] = thread(call_move_particles, i, nt, boxesCount, boxes, dt);
  }
  else{
      for(int i = 0; i < nt; ++i){
          thrds[i] = thread(dynamic_move_particles, i, boxes, dt);
      }
  }
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
/*
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
*/
  thrds = new thread[nt];
   if(!chunk_size){
  for(int i = 0; i < nt; ++i)
      thrds[i] = thread(call_update_particles, i, nt, boxesCount, boxes);
  }
  else{
      for(int i = 0; i < nt; ++i){
          thrds[i] = thread(dynamic_update_particles, i, boxes);
      }
  }
  for(int i = 0; i < nt; ++i)
    thrds[i].join();

  delete [] thrds;
 

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
    num_threads = nt;
    args = new ThreadArgs[nt];
    chunk_size = chunk;
    if(chunk){
	for(int k = 0; k < nt; k++){
	    int extra = (boxesCount/chunk)%nt;
	    int queue_size;
	    if(k < extra || (boxesCount%chunk != 0 && k == extra))
		queue_size = (boxesCount / (chunk * nt)) + 1;
	    else
		queue_size = boxesCount / (chunk * nt);
	    args[k].job_queue = (int*) malloc(sizeof(int)*queue_size);
	    //assert(args[k].job_queue);

	    args[k].queue_size = queue_size;
	    args[k].tid = k;

	    for(int j = 0; j < args[k].queue_size; j++){
		args[k].job_queue[j] = k * chunk_size + (chunk_size * nt) * j;
		if(j == queue_size - 1 && boxesCount % chunk_size != 0)
		    args[k].special = boxesCount % chunk_size;
		else
		    args[k].special = 0;
	    }
	}
    }
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
