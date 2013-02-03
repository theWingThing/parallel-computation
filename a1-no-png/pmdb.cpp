#include <pthread.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

// holds thread arguments
struct ThreadArgs
{
    int tid;          // thread id
    int NT;           // number of thread
    int** plot;       // mandelbrot points storage
    int64_t width;    // horizontal witdh
    int64_t height;   // horizontal height
    int* job_queue;   // contains index to the first row of each chunk queued
    int queue_size;   // num of elements in the array
    int chunk_size;   // size of chunks in number of rows
    int special;      // if this is the last chunk of an uneven input
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

    for(int k = 0; k < args->queue_size; k++)
    {
        int start_y = args->job_queue[k];
        int height  = args->chunk_size;
	if(args->special && k == args->queue_size-1)
	    height  = args->special;
        int start_x = 0;
        int width   = args->width;

        for(int i = 0; i < height; i++)
        {
            for(int j = 0; j < width; j++)
            {
                int x = start_x + j;
                int y = start_y + i;

                args->plot[y][x] =  
                    ComputeMandelbrotPoint (x, y, 0, 0);
            }
        }
    }
    return NULL;
}

void Mandelbrot_pthreads(int** pts, int dimX, int dimY, int numThreads, int chunkSize) 
{
    if(1 > numThreads)
    {
        fprintf(stderr, "invalid argument of numThreads = %d\n", numThreads);
        return;
    }

    int nThreads = numThreads;

    pthread_t* threads = new pthread_t[nThreads]; 
    ThreadArgs* args = new ThreadArgs[nThreads]; 
    
    if(chunkSize)
    {
        for(int k = 0; k < nThreads; k++)
        {
	    int extra = (dimY/chunkSize)%numThreads;
            int queue_size;
	    if(k < extra || (dimY % chunkSize != 0 && k == extra))
		queue_size = (dimY / (chunkSize * numThreads))+1;
	    else
		queue_size = dimY / (chunkSize * numThreads);
            /*
            if(dimY % (chunkSize * numThreads) && (nThreads - 1 == k))
            { 
                queue_size++;
            }
            */

            args[k].job_queue = (int*) malloc(sizeof(int)*(queue_size)); 
            assert(args[k].job_queue);

            args[k].queue_size = queue_size;

            for(int i = 0; i < queue_size; i++)
            {
                args[k].job_queue[i] = k * (chunkSize) + (chunkSize*nThreads) * i;
		if(i == queue_size - 1 && dimY % chunkSize != 0)
		    args[k].special = dimY % chunkSize;
		else
		    args[k].special = 0;
            }
            args[k].width = dimX;
            args[k].chunk_size = chunkSize;
            args[k].tid = k; 
            args[k].NT = nThreads; 
            args[k].plot = pts;
            assert(!pthread_create(&threads[k], NULL, fill_array_with_mandelbrotpoint, &args[k]));
        }
    }
    else
    {
        for(int k = 0; k < nThreads; k++)
        {
            args[k].job_queue = (int*) malloc(sizeof(int)); 
            assert(args[k].job_queue);
            args[k].queue_size = 1;
            args[k].job_queue[0] = (dimY/nThreads)*k;
            args[k].width = dimX;
            args[k].chunk_size = (dimY/nThreads);
            args[k].tid = k; 
            args[k].NT = nThreads; 
            args[k].plot = pts;
            assert(!pthread_create(&threads[k], NULL, fill_array_with_mandelbrotpoint, &args[k]));
        }
    }


    for(int t=0; t < nThreads; t++)
    {
        void *value;
        pthread_join(threads[t], &value);
        free(args[t].job_queue);
    }
}


