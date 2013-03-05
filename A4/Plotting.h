#ifndef __PLOTTING_H
#define __PLOTTING_H
#include <stdlib.h>
#include <cstdio>


class Plotter {
public:
    Plotter();
    ~Plotter();
    void updatePlot(double **U,  int niter, int m, int n, bool WAIT);


private:
    FILE *gnu_pipe;
};

#endif
