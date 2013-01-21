#include <stdio.h>
#include <math.h>

extern float re_min, im_min;
extern float scaleRe, scaleIm;
extern int maxIterations;

#ifdef _DOUBLE
#define _DOUBLE_ double
#else
#define _DOUBLE_ float
#endif

struct Complex
{
  _DOUBLE_ re;
  _DOUBLE_ im;
  Complex(_DOUBLE_ r, _DOUBLE_ i) {re = r; im = i; }
  Complex() {}
  _DOUBLE_ abs() {return sqrt(re*re + im*im); }
};

int ComputeMandelbrotPoint(Complex p)
{
	
  Complex z(0.0, 0.0);
  int iter = 0;
  _DOUBLE_ temp;
  
  while(z.re * z.re + z.im * z.im < 4.0 && iter < maxIterations) {
    temp = z.re * z.re - z.im * z.im + p.re;
    z.im = 2*z.re*z.im + p.im;
    z.re = temp;
    iter++;   
  }  
  return iter;
}

int ComputeMandelbrotPoint(int x, int y, int dimX, int dimY)
{
  Complex p;
  p.re = re_min + ((float) x * scaleRe);	
  p.im = im_min + ((float) y * scaleIm);
 
  return ComputeMandelbrotPoint(p);
}

void ComputeMandelbrot(int** points, int dimX, int dimY )
{
  for(int i=0; i<dimX; i++)
      for(int j=0; j<dimY; j++)		        		
	points[j][i] = ComputeMandelbrotPoint(i, j, dimX, dimY);               
}
