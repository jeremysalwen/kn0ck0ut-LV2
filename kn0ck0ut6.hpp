/* kn0ck0ut v0.5 by st3pan0va 2004 */

#ifndef __AGAIN_H
#define __AGAIN_H

#include <lv2plugin.hpp>
#include <fftw3.h>
#include "float.h"
#include "QuickTrig.h"


using namespace LV2;

#define MAX_FRAME_LENGTH 16384

enum
{
kCentre, kIn, kLoCut, kHiCut, kDecay, kBlur, kOut, kNumParams
};


class AKnockout;

class AKnockoutProgram
{
friend class AKnockout;
public:
	AKnockoutProgram();
	~AKnockoutProgram() {}
private:	
	float fCentre, fIn, fLoCut,fHiCut, fOut, fDecay, fBlur;
};


class AKnockout:public Plugin<AKnockout> {
	public:
	AKnockout(double srate);
	~AKnockout();
	void run(uint32_t numsamples);
	void suspend();
	CQuickTrig myQT;
	
protected:
	float fGain;
	char programName[32];
	
	
	
private:
	AKnockoutProgram programs;
	float *buffer;
	float fCentre, fIn, fLoCut, fOut, fHiCut, fDecay, fBlur;
	
	void do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, long gInit, float fGain, float fInGain, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract);
	void makelookup(int fftFrameSize);
	
	float* gInFIFO;
	float* gOutFIFO;
	float* gOutputAccum;
	float* gOutputBuffer;
	float* gAnaFreq;
	float* gAnaMagn;
	float* gInFIFO2;
	float* gAnaMagn2;
	float* gDecay;
	double* window;
	static long gRover;
	long gInit;

	
	fftwf_complex * gFFTworksp2;
	fftwf_complex * gFFTworksp;
	fftwf_plan forward_sp1;
	fftwf_plan forward_sp2;
	fftwf_plan backwards;

};



#endif
