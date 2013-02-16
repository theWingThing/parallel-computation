#include "Box.h"
#include "Bins.h"
#include <list>
#include <math.h>
#include <iostream>
#ifdef _DEBUG
#include <assert.h>
#endif
#if defined(_WIN32) || defined(_WIN64)
#define fmax max
#define fmin min
#pragma warning (disable:4996)
#define snprintf sprintf_s
#endif

using namespace std;

extern double cutoff;
extern double cutoff2;
extern double mass;
extern double size;
extern double min_r;

// May come in handy when debugging
inline bool is_inside(double x, double y, const Rect& r)
{
	return ( x >= r.x_start && x <= r.x_end && y >= r.y_start && y <= r.y_end);
}

void Box::AddParticle(particle_t* p)
{		
	boxParticles.push_back(p);
}

void Box::AddInboundParticle(particle_t* p)
{	
	inboundParticles.push_back(p);
}

// apply force from actors to reactors
void Box::apply_forces(ParticleList& reactors, ParticleList& actors) {
	for( ParticleIterator reactor = reactors.begin(); reactor != reactors.end(); reactor++ ) {
	if ((*reactor)->tag < 0)
	    continue;

        for (ParticleIterator actor = actors.begin(); actor != actors.end(); actor++ ){
            if ( *reactor == *actor)
                continue;
            double dx = (*actor)->x - (*reactor)->x;
            double dy = (*actor)->y - (*reactor)->y;
            double r2 = dx * dx + dy * dy;
            if( r2 > cutoff * cutoff )
                continue;
            r2 = fmax( r2, min_r*min_r );
            double r = sqrt( r2 );

            //  very simple short-range repulsive force
            double coef = ( 1 - cutoff / r ) / r2 / mass;
            (*reactor)->ax += coef * dx;
            (*reactor)->ay += coef * dy;
        }
    }
}

// Apply forces to the box
// The work divides into two parts as shown below
void Box::apply_forces()
{
	// reset force
	for( ParticleIterator p = boxParticles.begin(); p != boxParticles.end(); p++ ) 
        (*p)->ax = (*p)->ay = 0.0;
	// 1. Apply forces from particles inside this box
	//    to particles inside this Box
	apply_forces(boxParticles, boxParticles);

	// 2. Apply forces from particles from neighboring boxes
	//    to particles inside this Box
	for(int dx=-1; dx <= 1; ++dx)
		for(int dy=-1; dy <= 1; ++dy) {
			if( dx == 0 && dy == 0)
				continue;
			int x = I + dx;
			int y = J + dy;
			if( !(x < 0 || y < 0 || x >= bins->_nx || y >= bins->_ny) )
				apply_forces(boxParticles, bins->boxes[y*bins->_nx + x].boxParticles );
		}
}

//
//  integrate the ODE, advancing the positions of the particles
//
void Box::move_particles(double dt)
{
	// update particles' locations and velocities
	for( ParticleIterator particle = boxParticles.begin(); particle != boxParticles.end(); particle++ ) {
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
	if ((*particle)->tag < 0)
	    continue;
        (*particle)->vx += (*particle)->ax * dt;
        (*particle)->vy += (*particle)->ay * dt;
        (*particle)->x  += (*particle)->vx * dt;
        (*particle)->y  += (*particle)->vy * dt;

    //
    //  bounce from walls
    //
        while( (*particle)->x < 0 || (*particle)->x > size ) {
            (*particle)->x  = (*particle)->x < 0 ? -(*particle)->x : 2*size-(*particle)->x;
            (*particle)->vx = -(*particle)->vx;
        }
        while( (*particle)->y < 0 || (*particle)->y > size ) {
            (*particle)->y  = (*particle)->y < 0 ? -(*particle)->y : 2*size-(*particle)->y;
            (*particle)->vy = -(*particle)->vy;
        }		
    }	
	
}

//
// We check each particle to see if it moved outside the present box
// If, so we append to the inbound particle list for the new bin
//
void Box::UpdateParticlesBox()
{
	// Move to another Box if needed
	for( ParticleIterator particle = boxParticles.begin(); particle != boxParticles.end(); ) {
		int i = (int) ((*particle)->x / bins->boxWidth);
		int j  = (int) ((*particle)->y / bins->boxHeight);
		int newBoxID = j*bins->_nx + i;
		
		if(i != I || j != J) {
			bins->boxes[newBoxID].AddInboundParticle((*particle));
			particle = boxParticles.erase(particle);
		}
		else
			particle++;
	}
}

//
// After we've updated all the inbound lists
// We then make a new pass, incorporating into this bin,
// any particles on the inbound list
//

void Box::UpdateInboundParticles()
{
	for( ParticleIterator particle = inboundParticles.begin(); particle != inboundParticles.end(); particle++) {
		boxParticles.push_back(*particle);
	}
	inboundParticles.clear();
}
