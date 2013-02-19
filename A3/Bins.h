#ifndef _BINS_H
#define _BINS_H

#include "particle.h"
#include "Plotting.h"
#include <thread>

using namespace std;

class Box;

void call_apply_forces(int i, int nt, int boxesCount, Box* boxes);

class Bins
{
public:
	Bins(double size, int nx, int ny, int np, int n, particle_t* particles);
	~Bins();
	void SortParticles(int n, particle_t* particles);
	void apply_forces(int nt);
	void move_particles(double dt);
	void SimulateParticles(int nsteps, particle_t* particles, int n, int nt, int chunk, int nplot, bool imbal, double &uMax, double &vMax, double &uL2, double &vL2, Plotter *plotter, FILE *fsave, int nx, int ny, double dt );

	Box* boxes; // boxes inside the bins
	int boxesCount; // number of boxes
	double _size; 
	int _nx, _ny, _np, _n;
	double boxWidth, boxHeight; //dimensions of each box
	//thread *thrds;
};

#endif
