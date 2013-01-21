//////////// DO NOT MODIFY THIS FILE ////////////////////////

#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <string.h>
#include <fstream>

using namespace std;

double getTime();
void Mandelbrot_pthreads(int** points, int dimX, int dimY, int numThreads, int chunkSize );
void ComputeMandelbrot(int **points, int dimX, int dimY);
void ReportMdb(ofstream& logfile, bool runSerial, int threads, int chunkSize, bool verify, int** pts, int dimX, int dimY, double ptime);
void splot(int **image, int x, int y);
int **alloc2D(int m,int n);


int cmdLine(int argc, char** argv, int &numThreads, bool &runSerial, int &imgX, int &imgY, float& re_min, float& im_min, float& re_max, float& im_max, int &doPlot, int &chunkSize, char *filename, bool &verify, int &maxIterations);

//global variables for the range of Mandelbrot calculation
int maxIterations = 200;

float scaleRe, scaleIm;
float re_min, im_min, re_max, im_max;

int main(int argc, char* argv[]) 
{
	// number of worker threads for ComputeMandelbrot_pthreads
	int numThreads = 1;	
	bool runSerial = false;
	bool verify = false;
	int imgX = 1024;  // size of x/real dimension
	int imgY = 1024;  // size of y/imaginary dimension
	int doPlot = false;
	int chunkSize = 0; // 0 means BLOCK partitioning, otherwise it is the size of each chunk with cylic paritioning

	double sTime = 0.0;
	double pTime = 0.0;
	double t1, t2;
	char filename[100] = "\0";

	re_min=-2.5;
	im_min = -1;
	re_max=1;
	im_max = 1;

	// processing command line arguments
	if(!cmdLine(argc, argv, numThreads, runSerial, imgX, imgY, re_min, im_min, re_max, im_max, doPlot, chunkSize, filename, verify, maxIterations) )
	  return 0;

	scaleRe = (re_max - re_min)/imgX;
	scaleIm = (im_max - im_min)/imgY;
      
	// print out parameters
	cout << "max iterations: " << maxIterations << endl;
	cout << "image size: " << imgX << " x " << imgY << endl;
        cout << "Window: [" << re_min << " : " << re_max << "] x [";
        cout << im_min << " : " << im_max << "]" << endl;

	int **image = alloc2D(imgX,imgY);

	// parallel Mandelbrot
	t1 = getTime();

	if(runSerial)
	  ComputeMandelbrot(image, imgX, imgY);
	else
	//////////////////////////////////////////////////////////////////////
	//
	// Implement this routine, contained in pmdb.cpp
	// Do not change the function signature, which is defined in this file
	//
	  Mandelbrot_pthreads(image, imgX, imgY, numThreads, chunkSize);
	//////////////////////////////////////////////////////////////////////
	t2 = getTime();
	pTime += t2 - t1;
	
	ofstream logfile("log");
	ReportMdb(logfile, runSerial, numThreads, chunkSize, verify, image, imgX, imgY, pTime);
	
	// write the visualization of Mandelbrot set to a png image file
	if (doPlot){
	  splot(image, imgX, imgY);
	}

	delete [] image;

	return 0;
}
