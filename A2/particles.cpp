// All the parallelism  arises in code within this file
//
// The following two functions update the state of the particles,
// and will perform the computation in parallel across multiple threads
// The openMP code specifies how you'll parallelize the code
// including options (set on the make command line) to
// enable dynamic scheduling and chunking
// Note that dynamic scheduling and chunking are specified at compile time
// in OpenMP, but at run time in your multithreaded code
// Keep this in mind when building the code and don't use the compile time
// flags dyn= chunk= unless you are using OpenMP
// You should not enable OpenMP and C++11 threads in the same code
// as the behavior is unpredictable
//
//    apply_forces( )
//    move_particles( )
//
// Your performance optimizations may change the order in
// which arithmetic gets done, affecting the results
// This is acceptable, as discussed in class, so long as the differences
// are to within roundoff errors 
// To assess correctness, examine the values of vMax and vL2
// reported at the simulation's end
///
//
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include "particle.h"
#include "common.h"
#include "Plotting.h"
#include <iostream>
#ifdef  _OPENMP
#include "omp.h"
#endif
using namespace std;

extern double size;
extern double dt;

mutex done;
bool notDone = true;

int count = 0;
int num_threads = 0;

//void barrier(void);

class SelfScheduler
{ 
    private:
        int problem_size;
        int chunk_size;
        int num_processors;
        int counter;
        mutex critical_section, m, d;
        condition_variable cond_var, cond_done;

        SelfScheduler(SelfScheduler& s){};

    public:
        SelfScheduler():problem_size(0), chunk_size(0),
            num_processors(0), counter(0){};

        SelfScheduler(int n,int P, int chunk):problem_size(n), 
            chunk_size(chunk), num_processors(P), counter(0){};

        ~SelfScheduler(){};

        void unlock()
        {
            std::unique_lock<std::mutex> lock(m);
            lock.unlock();
        }

        void begin_work()
        {
            cond_var.notify_all();
        }

        void cond_wait()
        {
            std::unique_lock<std::mutex> lock(m);
            cond_var.wait(lock);
        }

        bool isDone()
        {
            return counter > problem_size;
        }

        void wait_for_done()
        {        
            std::unique_lock<std::mutex> done_lock(d);
            cond_done.wait(done_lock);
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

SelfScheduler *scheduler, *scheduler1;

void imbal_particles(particle_t *particles, int n);

void barrier(SelfScheduler* sche)
{ 
    count++;
    
    if(count == num_threads || sche->isDone())
    {
        count = 0;
        sche->begin_work();
    }
    else 
    {
        sche->cond_wait();
    }
    sche->unlock();
}

void callback(particle_t* particles, SelfScheduler* self_sche, int tid, int NT, int n, void (*func)(particle_t* particles, int min, int max, int n))
{        
    int interval = n/NT;
    int next_tid = tid + 1;

    int min = tid*interval;
    int max = next_tid*interval;
    bool notDone = true;

    while(notDone)
    {
        if(self_sche)
        { 
            bool doMore = true;

            while(doMore)
            {      
                doMore = self_sche->getChunk(min,max); 
                barrier(self_sche);
                if(doMore)
                { 
                    func(particles, min, max, n);
                }
            }
        }
        else
        {
            func(particles, min, max, n);
        }

        //wait here untill the signal
        self_sche->cond_wait();
    }
}     

thread* spawn_threads (particle_t* pParticles, SelfScheduler* self_sche, int n, void (* func) (particle_t* particles, int min, int max, int n))
{
    thread *thrds = new thread[num_threads];

    for(int t = 0; t < num_threads; t++)
    {
        thrds[t] = thread(callback, pParticles, self_sche, t, num_threads, n, func);
    }

    return thrds;
}

// spawn and join threads
// with call back thread function
void thread_controller (particle_t* pParticles, SelfScheduler* self_sche, int n, void (* func) (particle_t* particles, int min, int max, int n))
{
    // wake up threads
    if(scheduler)
    {
        scheduler->begin_work();
    }

    self_sche->wait_for_done();
    self_sche->reset();
}

void compute_forces (particle_t* particles, int min, int max, int n)
{
    for (int i = min; i < max; i++) 
    {
        particles[i].ax = particles[i].ay = 0;

        if ((particles[i].vx != 0) || (particles[i].vy != 0))
        {    
            for (int j = 0; j < n; j++ )
            {
                if (i == j)
                {
                    continue;
                }

                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double r2 = dx * dx + dy * dy;

                if( r2 > cutoff*cutoff )
                {
                    continue;
                }

                r2 = fmax( r2, min_r*min_r );
                double r = sqrt( r2 );

                //  very simple short-range repulsive force
                double coef = ( 1 - cutoff / r ) / r2 / mass;
                particles[i].ax += coef * dx;
                particles[i].ay += coef * dy;
            }
        }
    }
}
//  compute force between two particles
//  You should not modify the numerical computations
//  other than to optimize them
//  Be careful in optimizing, as some optimizations
//  can subtly affect the computed answers
void apply_forces( particle_t* particles, int n)
{
/*
#ifdef DYN
#if CHUNK
#pragma omp parallel for shared(particles,n) schedule(dynamic,CHUNK)
#else
#pragma omp parallel for shared(particles,n) schedule(dynamic)
#endif
#else
#pragma omp parallel for shared(particles,n)
#endif
*/

    thread_controller (particles, scheduler, n, compute_forces);

}

void compute_particles_loc(particle_t* particles, int min, int max, int n)
{
    (void) n;

    for( int i = min; i < max; i++ ) 
    {
        //  slightly simplified Velocity Verlet integration
        //  conserves energy better than explicit Euler method
        particles[i].vx += particles[i].ax * dt;
        particles[i].vy += particles[i].ay * dt;
        particles[i].x  += particles[i].vx * dt;
        particles[i].y  += particles[i].vy * dt;

        //  bounce off the walls
        while( particles[i].x < 0 || particles[i].x > size ) 
        {
            particles[i].x  = particles[i].x < 0 ? -particles[i].x : 2*size-particles[i].x;
            particles[i].vx = -particles[i].vx;
        }

        while( particles[i].y < 0 || particles[i].y > size ) 
        {
            particles[i].y  = particles[i].y < 0 ? -particles[i].y : 2*size-particles[i].y;
            particles[i].vy = -particles[i].vy;
        }
    }
}
//
//  integrate the ODE, advancing the positions of the particles
//
void move_particles( particle_t* particles, int n)
{
    thread_controller (particles, scheduler1, n, compute_particles_loc);

}

// This is the main driver routine that runs the simulation
void SimulateParticles (int nsteps, particle_t *particles, int n, int nt, int chunk, int nplot, bool imbal, double &uMax, double &vMax, double &uL2, double &vL2, Plotter *plotter, FILE *fsave )
{
    thread* thrds, *thrds1;
    // set global variable of number of threads
    num_threads = nt;

    done.lock();
    notDone = true;
    done.unlock();

    if(0 < chunk)
    {
        scheduler = new SelfScheduler(n, nt, chunk);
        scheduler1 = new SelfScheduler(n, nt, chunk);
    }
    else
    {
        scheduler = NULL;
        scheduler1 = NULL;
    }

    thrds = spawn_threads (particles, scheduler, n, compute_forces);
    thrds1 = spawn_threads(particles, scheduler1, n, compute_particles_loc);

    for(int step = 0; step < nsteps; step++ ) 
    {

        //  compute forces
        apply_forces(particles,n);

        if(scheduler)
        {
            scheduler->reset();
        }

        // If we asked for an imbalanced distribution
        if (imbal)
            imbal_particles(particles,n);
        //  Debugging output
        //  list_particles(particles,n);
        

        //  move particles
        move_particles(particles, n);

        if (nplot && ((step % nplot ) == 0))
        {
            // Computes the absolute maximum velocity 
            VelNorms(particles,n,uMax,vMax,uL2,vL2);
            plotter->updatePlot(particles,n,step,uMax,vMax,uL2,vL2);
        }

        VelNorms(particles,n,uMax,vMax,uL2,vL2);

        //  save if necessary
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    
    done.lock();
    notDone = false;
    done.unlock();

    //join_threads
    scheduler->begin_work();

    for(int t = 0; t < num_threads; t++)
    { 
        thrds[t].join();
    }

    scheduler1->begin_work();

    for(int t = 0; t < num_threads; t++)
    { 
        thrds1[t].join();
    }

    delete [] thrds;
    delete [] thrds1;

    if(scheduler)
    {
        delete scheduler;
    }
    if(scheduler1)
    {
        delete scheduler1;
    }
}

