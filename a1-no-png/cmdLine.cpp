#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

int cmdLine(int argc, char** argv, int &numThreads, bool &runSerial, int &imgX, int &imgY, float& re_min, float& im_min, float& re_max, float& im_max, int &doPlot, int &chunkSize, char *filename, bool &verify, int &maxIterations)
{
  for( int i=1; i < argc; ++i) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 't':
	if( i+1 >= argc) {
	  cout << "missing argument with -t. Usage: mdb -t numThreads" << endl;
	  return 0;
	}
	numThreads = atoi(argv[i+1]);
	i++;
	break;

      case 'i':
	if( i+1 >= argc) {
	  cout << "missing argument with -i. Usage: mdb -i maxIterations" << endl;
	  return 0;
	}
	maxIterations = atoi(argv[i+1]);
	i++;
	break;

      case 'p':
	doPlot = true;
	i++;
	break;

      case 'x':
	if( i+1 >= argc) {
	  cout << "missing argument with -x. Usage: mdb -x xRange" << endl;
	  return 0;
	}
	imgX = atoi(argv[i+1]);
	i++;
	break;

      case 'y':
	if( i+1 >= argc) {
	  cout << "missing argument with -y. Usage: mdb -y yRange" << endl;
	  return 0;
	}
	imgY = atoi(argv[i+1]);
	i++;
	break;

      case 'b':
	if( i+1 >= argc) {
	  cout << "missing argument with -b. Usage: mdb -b real_min real_max imaginary_min imaginary_max" << endl;
	  return 0;
	}
	re_min = atof(argv[++i]);
	if( i+1 >= argc) {
	  cout << "missing argument with -b. Usage: mdb -b real_min real_max imaginary_min imaginary_max" << endl;
	  return 0;
	}
	re_max = atof(argv[++i]);
	if( i+1 >= argc) {
	  cout << "missing argument with -b. Usage: mdb -b real_min real_max imaginary_min imaginary_max" << endl;
	  return 0;
	}
	im_min = atof(argv[++i]);
	if( i+1 >= argc) {
	  cout << "missing argument with -b. Usage: mdb -b real_min real_max imaginary_min imaginary_max" << endl;
	  return 0;
	}
	im_max = atof(argv[++i]);
	break;

	// run the serial version
      case 's':
	runSerial = true;
	break;

      case 'c':
	chunkSize = atoi(argv[i+1]);
	i++;
	break;
	//output to a png file
      case 'o':
	if( i+1 >= argc) {
	  cout << "missing filename with argument -o. Usage: mdb -o \"filename\" " << endl;
	  return 0;
	}
	strcpy(filename, argv[i+1]);
	i++;
	break;

      case 'v':
	verify = true;
	break;
      }
    }
  }
  if(numThreads > 1 && runSerial) {
    cout << "incompatible options: -t -s cannot be used together. " << endl;
    return 0;
  }
  return 1;
}
