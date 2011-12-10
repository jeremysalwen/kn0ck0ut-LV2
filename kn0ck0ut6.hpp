/*
	Copyright 2004 St3pan0va, 2011 Jeremy Salwen
	
	This file is part of Kn0ck0ut-LV2
	
    Kn0ck0ut-LV2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __AGAIN_H
#define __AGAIN_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <lv2plugin.hpp>
#include <fftw3.h>
#include "float.h"
#include "QuickTrig.h"
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
	inline void knockout(float origmag,float origphase, fftwf_complex * __restrict workspace, float* __restrict magnbuff, float* __restrict decaybuff, float coef, int k, float fDecayRate, int iBlur);
	void do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata1,float *outdata2, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract, bool consider_phase);
	inline void sumIntoCircularBuffer(float* __restrict outaccum,float dOutFactor,long outAccumIndex,long stepSize,long fftFrameSize);
	void makelookup(int fftFrameSize);
	inline float phaseToFrequency(float phase, int k,float dOversampbytwopi, float expct);
	
	float*  __restrict gInFIFO;
	float* __restrict gInFIFO2;
	float* __restrict gOutputAccum;
	float* __restrict gOutputAccum2;
	float* __restrict FFTRealBuffer;
	float* __restrict gAnaPhase1;
	float* __restrict gAnaPhase2;
	float* __restrict gAnaMagn;
	float* __restrict gAnaMagn2;
	float* __restrict gDecay;
	float* __restrict gDecay2;
	float* __restrict window;

	long gRover;
	long samples_needed_in_buffer;
	long outAccumIndex;
	long copiesremaining;
	
	fftwf_complex * gFFTworksp2;
	fftwf_complex * gFFTworksp;
	fftwf_plan forward_sp1;
	fftwf_plan forward_sp2;
	fftwf_plan backward_sp1;
	fftwf_plan backward_sp2;

};



#endif
