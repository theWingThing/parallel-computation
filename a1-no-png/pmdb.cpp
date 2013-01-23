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
    if(-1 != args->start_x)
    {
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
    // for block partition, there will be alway four blocks
    // the four blocks are defined as a retangular with
    // following dimensions: dimY/4 X dimX
    // Each thread will be given a number of block 
    // untill all threads have work
    // if there is one thread all work is givne the thread
    // if there is two threads each thread gets the half
    // if there is three threads each thread get one
    // and the first one get one more
    // if there are four threads
    // each gets one
    // threads more than 4 do not receive anywork
    
    for(int k = 0; k < nThreads; k++){
        args[k].start_x = 0;
        args[k].start_y = k/nThreads;
        args[k].width = dimX;
        args[k].height = (k+1)/nThreads;
    }
/*
    switch(nThreads)
    {
        case 1:
            args[0].start_x = 0;
            args[0].start_y = 0;
            args[0].width   = dimX; 
            args[0].height  = dimY; 
            break;
        case 2:
            args[0].start_x = 0;
            args[0].start_y = 0;
            args[0].width   = dimX; 
            args[0].height  = dimY/ 2; 

            args[1].start_x = 0;
            args[1].start_y = dimY;
            args[1].width   = dimX; 
            args[1].height  = dimY/ 2; 
            break;    
        case 3:
            args[0].start_x = 0;
            args[0].start_y = 0;
            args[0].width   = dimX; 
            args[0].height  = dimY / 4; 

            args[1].start_x = 0;
            args[1].start_y = dimY / 4;
            args[1].width   = dimX; 
            args[1].height  = dimY / 4; 

            args[2].start_x = 0;
            args[2].start_y = dimY / 2;
            args[2].width   = dimX; 
            args[2].height  = dimY/2; 
            break;
        default:
            args[0].start_x = 0;
            args[0].start_y = 0;
            args[0].width   = dimX; 
            args[0].height  = dimY / 4; 

            args[1].start_x = 0;
            args[1].start_y = dimY / 4;
            args[1].width   = dimX; 
            args[1].height  = dimY / 4; 

            args[2].start_x = 0;
            args[2].start_y = dimY / 2;
            args[2].width   = dimX; 
            args[2].height  = dimY / 4; 

            args[3].start_x = 0;
            args[3].start_y = dimY / 3;
            args[3].width   = dimX; 
            args[3].height  = dimY / 4; 
            break;    
    }
*/
    for(int t = 0; t < nThreads; t++)
    {
        args[t].tid = t; 
        args[t].NT = nThreads; 
        args[t].plot = pts;
        if(t > 3)
        {
            args[t].start_x = -1;
            args[t].start_y = -1;
            args[t].width   = -1;
            args[t].height  = -1;
        }

        assert(!pthread_create(&threads[t], NULL, fill_array_with_mandelbrotpoint, &args[t]));
    }

    for(int t=0; t < nThreads; t++)
    {
        void *value;
        pthread_join(threads[t], &value);
    }
}


