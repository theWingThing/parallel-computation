#ifndef _ST_Box_H
#define _ST_Box_H

#include "particle.h"
#include <vector>
#include <list>

using namespace std;

// We use ParticleLists to represent the list of particles
// owned by a bin
//
// We use a particle_t* to store the global list of particles
//
// Both refer to the same set of data in memory
//
typedef vector<particle_t*> ParticleList;
typedef ParticleList::iterator ParticleIterator;

struct Rect
{
	double x_start, x_end, y_start, y_end;
};

class Bins;

class Box
{
public:  
  ParticleList boxParticles; // the particles inside this box
  ParticleList inboundParticles; // the particles coming into this box 
  				 // We use this list to keep track of
				 // all incoming migrating particles
  int I, J;  	// Index for this box
  Bins* bins; // the bins object used to access other boxes
    
  void AddParticle(particle_t* p);
  void AddInboundParticle(particle_t* p);
  void apply_forces();
  void apply_forces(ParticleList& reactors, ParticleList& actors);
  void move_particles(double dt);
  void UpdateParticlesBox();
  void UpdateInboundParticles();
};


#endif
