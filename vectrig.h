#ifndef VECTRIG_H
#define VECTRIG_H

//All this code is borrowed from vectrig, which I also wrote.
#define FLTEND(x) x##f
#define FLTTYPE float

#define HALFPI  1.570796327
#define PI 3.141592653589793238462643
#define TWOPI (2*PI)


static const float A=0xA2F9836Cp-33;

//Modified to
static inline float bringtorange(float value) {
	value-=(value>0)*TWOPI;
	value+=PI;
	value=(value>HALFPI)?(PI-value):value;
	return (value<-HALFPI)?value+PI:-value;
}


#define B91  FLTEND(-1.6666658777408987333e-1)
#define B92  FLTEND(8.3330585579660222884e-3)
#define B93  FLTEND(-1.9809580202098392708e-4)
#define B94  FLTEND(2.6065544878032086552e-6)

static inline float b9series(float value) {
   value=bringtorange(value);
   FLTTYPE sqr=value*value;
   return value *(1 + sqr*(B91+sqr*(B92 +sqr*(B93 +sqr*B94))));
}

#endif