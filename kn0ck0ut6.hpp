/* kn0ck0ut v0.5 by st3pan0va 2004 */

#ifndef __AGAIN_H
#define __AGAIN_H

#include "float.h"
#include "QuickTrig.h"

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
	char name[24];
};


class AKnockout : public AudioEffectX
{
public:
	AKnockout(audioMasterCallback audioMaster);
	~AKnockout();

	virtual void process(float **inputs, float **outputs, long sampleFrames);
	virtual void processReplacing(float **inputs, float **outputs, long sampleFrames);
	virtual void setProgram(long program);
	virtual void setProgramName(char *name);
	virtual void getProgramName(char *name);
	virtual void setParameter(long index, float value);
	virtual float getParameter(long index);
	virtual void getParameterLabel(long index, char *label);
	virtual void getParameterDisplay(long index, char *text);
	virtual void getParameterName(long index, char *text);
	virtual float getVu();
	virtual void suspend();
	virtual void resume();
	CQuickTrig myQT;
	
protected:
	float fGain;
	char programName[32];
	
	
	
private:
	void setDelay(float delay);
	
	AKnockoutProgram *programs;
	float *buffer;
	float fCentre, fIn, fLoCut, fOut, fHiCut, fDecay, fBlur;
	float vu;
	long delay;
	
	void do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, long gInit, float fGain, float fInGain, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract);
	void makelookup(int fftFrameSize);
	inline void smsFft(float *fftBuffer, long fftFrameSize, long sign);
	
	float* gInFIFO;
	float* gOutFIFO;
	float* gFFTworksp;
	float* gOutputAccum;
	float* gAnaFreq;
	float* gAnaMagn;
	float* gInFIFO2;
	float* gFFTworksp2;
	float* gAnaMagn2;
	float* gDecay;
	double* window;
	static long gRover;
	long gInit;
	
	
};



#endif
