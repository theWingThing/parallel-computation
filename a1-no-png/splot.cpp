/* **********************************************************
 *  Author : Urvashi R.V. [04/06/2004]
 * 	Modified by Scott Baden [10/8/06]
 * 	Modified by Pietro Cicoctti [10/8/08]
 * 
 *************************************************************/

#include <stdio.h>
#include <math.h>

#include <algorithm>
#include <fstream>
#include "util.h"

/* Function to plot the 2D array
 * 'gnuplot' is instantiated via a pipe and 
 * the values to be plotted are passed through, along 
 * with gnuplot commands */

using namespace std;

extern float re_min, re_max, im_min, im_max;

void ConvertMdbToRGB(int** image, int max, int x, int y,ColorB** rgb);

FILE* gnu = NULL;

void splot(int **image, int x, int y)
{
  if(!gnu) gnu = popen("gnuplot -persist","w");
  int mx= -1, mn = 32768;
  for (int j=0; j<y; j++)
     for (int i=0; i<x; i++){
        if (image[j][i] > mx)
            mx = image[j][i];
        if (image[j][i] < mn)
            mn = image[j][i];
  }
  /*  for (int j=0; j<y; j++){
       for (int i=0; i<x; i++)
           printf("%d %d %d\n", i, j, image[j][i]);
       fprintf(gnu,"\n");
       }*/
  fprintf(gnu,"\n");
  fprintf(gnu,"\n");
  fprintf(gnu,"set size square\n");
  fprintf(gnu,"set key off\n");
  fprintf(gnu,"set pm3d map\n");
  fprintf(gnu,"set palette defined (-3 \"blue\", 0 \"white\", 1 \"red\")\n");
  // fprintf(gnu, "set xlabel 'Real'\n");
  //fprintf(gnu, "set ylabel 'Imaginary'\n");
    /* Various color schemes
     * fprintf(gnu,"set palette rgbformulae 22, 13, 31\n");
     * fprintf(gnu,"set palette rgbformulae 30, 31, 32\n");
    */

  fprintf(gnu,"splot [0:%d] [0:%d][%d:%d] \"-\"\n",x-1,y-1,mn,mx);
  for (int j=0; j<y; j++){
       for (int i=0; i<x; i++)
           fprintf(gnu,"%d %d %d\n", i, j, image[j][i]);
       fprintf(gnu,"\n");
  }
  fprintf(gnu,"e\n");
  fflush(gnu);
  return;
}
