/* kn0ck0ut v0.5 by st3pan0va 2004 */

#ifndef __AGAIN_H
#define __AGAIN_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <lv2plugin.hpp>
#include <fftw3.h>
#include "float.h"
#include "QuickTrig.h"
#include "vectrig.h"
#include "kn0ck0ut.peg"


using namespace LV2;


class AKnockout:public Plugin<AKnockout> {
	public:
	AKnockout(double srate);
	~AKnockout();
	void run(uint32_t numsamples);
	void activate();
	void deactivate();
	CQuickTrig myQT;
	
	
private:
	void clearBuffers();
	void AllocateNewBuffers(unsigned int fftsize);
	void FreeOldBuffers();
	unsigned int goverlap;
	unsigned int gfftSize;
	double sampleRate;
	
	void do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract);
	void makelookup(int fftFrameSize);
	
	float*  __restrict gInFIFO ;
	float* __restrict gOutputAccum;
	float* __restrict FFTRealBuffer;
	float* __restrict gAnaFreq;
	float* __restrict gAnaMagn;
	float* __restrict gInFIFO2;
	float* __restrict gAnaMagn2;
	float* __restrict gDecay;
	float* __restrict window;

	long gRover;
	long samples_needed_in_buffer;
	long outAccumIndex;
	long copiesremaining;
	
	fftwf_complex * gFFTworksp2;
	fftwf_complex * gFFTworksp;
	fftwf_plan forward_sp1;
	fftwf_plan forward_sp2;
	fftwf_plan backwards;

};



#endif
