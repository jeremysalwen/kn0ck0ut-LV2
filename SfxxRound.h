// Copyright (c) 1999-2003 Robin Davies.
// Placed into the public domain. 
// All uses permitted without attribution.
#ifndef SFXXROUND_H
#define SFXXROUND_H


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


#ifdef _M_IX86

#include "SfxxInline.h"

const double _TWOPOW32 = ((double)(65536.0 * 65536.0));

inline double _TWOPOW(int N) 
{
	return (((N) < 32) ? ((double)(1UL << (N))) : 
                        ((double)(1UL << (N - 32)) * _TWOPOW32))
							;
}

const double g_MagicRoundingConstant = _TWOPOW(52) + _TWOPOW(51);

/* 
   Declared as a class to ensure correct 8-byte alignment of 
   mydtemp on x86. Don't use this off the stack!

*/
class CImplementsRounding {
public:
  FORCEINLINE int Round(double x) const {
	((CImplementsRounding *)this)->mydtemp = (double)(x) + g_MagicRoundingConstant;
	return *(int*)&mydtemp;
  };
  FORCEINLINE int Floor(double x)const  {
	((CImplementsRounding*)this)->mydtemp = (double)(x) + g_MagicRoundingConstant;
	int val = *(int*)&mydtemp;
    if (val > x) --val;
    return val;
  } ;
  FORCEINLINE int Ceil(double x) const  {
	((CImplementsRounding*)this)->mydtemp = (double)(x) + g_MagicRoundingConstant;
	int val = *(int*)&mydtemp;
    if (val < x) ++val;
    return val;
  };


protected:
  double mydtemp;
};

#else
//#error Only works on x86.
// but then again, sensible platforms should have to do this
// kind of crap to convert a float to an integer.
// Try this. Should work.
class CImplementsRounding {
public:
    static int Round(double d) { return (int)d; }
};

#endif

#endif
