#include <pthread.h>

// Implemented for you in smdb.h
int ComputeMandelbrotPoint(int x, int y, int dimX, int dimY);

// Implement the thread function which will be called by the spawned threads
//
//
// This function is called by the main thread
// It will spawn and join threads, calling on your thread function
// to perform the former task
//
void Mandelbrot_pthreads(int** pts, int dimX, int dimY, int numThreads, int chunkSize) 
{
  //////////// add your code here //////////////////////////////////
}


