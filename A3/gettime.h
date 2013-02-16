#ifndef _GETTIME_H
#define _GETTIME_H

#include < time.h >
#include <windows.h> 

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
int gettimeofday(struct timeval *tv, struct timezone *tz);

#endif