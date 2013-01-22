#include <pthread.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

// holds thread arguments
struct ThreadArgs
{
    int tid;          // thread id
    int NT;           // number of thread
    int** plot;       // mandelbrot points storage
    int start_x;      // start coordinate of x
    int start_y;      // start coordinate of y
    int64_t width;    // horizontal witdh
    int64_t height;   // horizontal height
};

// Implemented for you in smdb.h
// this function calcullate mandel brot point given x and y coordinates
int ComputeMandelbrotPoint(int x, int y, int dimX, int dimY);

// Implement the thread function which will be called by the spawned threads
//
//
// This function is called by the main thread
// It will spawn and join threads, calling on your thread function
// to perform the former task
//
// given a partition size, 
void* fill_array_with_mandelbrotpoint(void* arg)
{
    ThreadArgs *args = (ThreadArgs*) arg;

    for(int i = 0; i < args->height; i++)
    {
        for(int j = 0; j < args->width; j++)
        {
            int x = args->start_x + j;
            int y = args->start_y + i;

            args->plot[x][y] =  
                ComputeMandelbrotPoint (x, y, 0, 0);
        }
    }

    return NULL;
}

void Mandelbrot_pthreads(int** pts, int dimX, int dimY, int numThreads, int chunkSize) 
{
    //////////// add your code here //////////////////////////////////
    //// creates #numThreads of threads
    //// for test just evenly divde the works among the threads
    
    if(1 > numThreads)
    {
        fprintf(stderr, "invalid argument of numThreads = %d\n", numThreads);
        return;
    }

    int nThreads = numThreads;

    pthread_t* threads = new pthread_t[nThreads]; 
    ThreadArgs* args = new ThreadArgs[nThreads]; 

    // thread creation
    for(int t=0; t<nThreads; t++)
    {
        
        args[t].tid = t; 
        args[t].NT = nThreads; 

        // we need to figure out how to partition the plot spaces 
        // to each threads 
        args[t].start_x = 0;
        args[t].start_y = 0;
        args[t].width   = 0; 
        args[t].height  = 0; 

        assert(!pthread_create(&threads[t], NULL, fill_array_with_mandelbrotpoint, &args[t]));
    }

    for(int t=0; t < nThreads; t++)
    {
        void *value;
        pthread_join(threads[t], &value);
    }
}


