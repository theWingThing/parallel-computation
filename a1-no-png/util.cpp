#include "util.h"
#include <iostream>
#include <cstring>
#include <assert.h>

using namespace std;

// Allocate a 2D array
int **alloc2D(int m,int n){

   int **E;
   int nx=n+1, ny=m+1;
   E = (int**)malloc(sizeof(int*)*ny + sizeof(int)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++) E[j] = (int*)(E+ny) + j*nx;
   return(E);
}

