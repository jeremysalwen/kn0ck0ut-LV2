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
#include "kn0ck0ut6.hpp"



void CQuickTrigConsts::Initialize() // function to initialise quick sin & cos
{

	int i;
	for (i = 0; i <= kMsTableSize; ++i) {
		double phi = (2*PI*i/kMsTableSize);
		mMsBitsTable[i].msin = (float)sin(phi);
		mMsBitsTable[i].mcos = (float)cos(phi);
	}
	for (i = 0; i <= kLsTableSize; ++i) {
		double phi = (2*PI*i/(kMsTableSize*1.0*kLsTableSize));
		mLsBitsTable[i].msin = (float)sin(phi);
		mLsBitsTable[i].mcos = (float)cos(phi);
	}
}

// Initialize QuickTrig constants with a global constructor.

class CQuickTrigInitialize: CQuickTrigConsts {
	public:
		CQuickTrigInitialize() {
			CQuickTrigConsts::Initialize();
		}
};

CQuickTrigInitialize gQuickTrigInitialize;


//-----------------------------------------------------------------------------
AKnockout::AKnockout(double rate) : Plugin<AKnockout>(p_n_ports)
{
	sampleRate=rate;
	gfftSize=8192;
	goverlap=8;
}
void AKnockout::activate() {
	AllocateNewBuffers(gfftSize);
	clearBuffers();
}

void AKnockout::deactivate() {// delete buffers when deactivated
	FreeOldBuffers();
}
//-----------------------------------------------------------------------------------------
AKnockout::~AKnockout() 
{
}

void AKnockout::AllocateNewBuffers(unsigned int fftSize) {
	unsigned int fftSize2=fftSize/2+1;
	gInFIFO = new float [fftSize];
	FFTRealBuffer=(float*)fftwf_malloc(sizeof(float)*fftSize);
	gFFTworksp = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fftSize2);
	gOutputAccum = new float [fftSize];
	gOutputAccum2 = new float [fftSize];
	gAnaPhase1 = new float [fftSize2];
	gAnaPhase2 = new float [fftSize2];
	gAnaMagn = new float [fftSize2];
	gInFIFO2 = new float [fftSize];
	gFFTworksp2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fftSize2);
	gAnaMagn2 = new float [fftSize2];
	gDecay = new float [fftSize2];
	gDecay2 = new float [fftSize2];
	window = new float [fftSize];

	forward_sp1= fftwf_plan_dft_r2c_1d(fftSize, FFTRealBuffer , gFFTworksp,
                                    FFTW_ESTIMATE);
	forward_sp2= fftwf_plan_dft_r2c_1d(fftSize,FFTRealBuffer, gFFTworksp2,
                                    FFTW_ESTIMATE);	
	backward_sp1=fftwf_plan_dft_c2r_1d(fftSize, gFFTworksp, FFTRealBuffer,
                                    FFTW_ESTIMATE);
	backward_sp2=fftwf_plan_dft_c2r_1d(fftSize, gFFTworksp2, FFTRealBuffer,
                                    FFTW_ESTIMATE);
	makelookup(fftSize);
}
void AKnockout::FreeOldBuffers() {
	delete[] gInFIFO;
	fftwf_free(FFTRealBuffer);
	fftwf_free(gFFTworksp);
	delete[] gOutputAccum;
	delete[] gOutputAccum2;
	delete[] gAnaPhase1;
	delete[] gAnaPhase2;
	delete[] gAnaMagn;
	delete[] gAnaMagn2;
	delete[] gInFIFO2;
	fftwf_free(gFFTworksp2);
	delete[] gDecay;
	delete[] gDecay2;
	delete[] window;
}


//-----------------------------------------------------------------------------------------
void AKnockout::clearBuffers()
{
	unsigned int fftSize=gfftSize;
	unsigned int fftSize2=fftSize/2+1;
	memset(gInFIFO, 0, fftSize*sizeof(float));
	memset(gFFTworksp, 0, fftSize2*sizeof(fftwf_complex));
	memset(gOutputAccum, 0, fftSize*sizeof(float));
	memset(gOutputAccum2, 0, fftSize*sizeof(float));
	memset(gAnaPhase1, 0, fftSize2*sizeof(float));
	memset(gAnaPhase2, 0, fftSize2*sizeof(float));
	memset(gAnaMagn, 0, fftSize2*sizeof(float));
	memset(gInFIFO2, 0, fftSize*sizeof(float));
	memset(gFFTworksp2, 0, fftSize2*sizeof(fftwf_complex));
	memset(gAnaMagn2, 0, fftSize2*sizeof(float));
	memset(gDecay,0,fftSize2*sizeof(float));
	memset(gDecay2,0,fftSize2*sizeof(float));
	
	long stepSize=fftSize/goverlap;
	gRover=0;
	samples_needed_in_buffer=stepSize;
	outAccumIndex=stepSize;
	copiesremaining=0;
}
//assumes desired and fftsize are multiples of 4.
int calcOsampFromFFTSize(unsigned long desired, unsigned long fftsize) {
	if(fftsize%desired==0) {
		return desired;
	}
	if(desired>fftsize) {
		return fftsize;
	}
	if(desired<4) {
		return 4;
	}
	unsigned int i;
	for(i=desired; fftsize%i!=0; i++);
	if(i==fftsize) {
		for(i=desired; fftsize%i!=0; i--);
	}
	return i;
}

unsigned long findBestFFTSize(unsigned long desired) {
	//Next highest multiple of four.
	unsigned long lastbits=desired & 0x0000000000000003;
	if(lastbits) {
		return ((desired >>2)+1)<<2;
	} else {
		return desired;
	}
}
#define CLAMP_PORT(TYPE, VAR,PORT_SYMBOL)\
	TYPE VAR = (TYPE) (*p(PORT_SYMBOL));\
	if(VAR<p_ports[PORT_SYMBOL].min) {\
		VAR = p_ports[PORT_SYMBOL].min;\
	}\
	if(VAR>p_ports[PORT_SYMBOL].max) {\
		VAR = p_ports[PORT_SYMBOL].max;\
	}\
//-----------------------------------------------------------------------------------------
void AKnockout::run(uint32_t sampleFrames)
{
	CLAMP_PORT(int,loCut,p_lowcut)
	

	int centre = (*p(p_mode)>0);

	//Fill in fft size here.
	
	int hiCut = int (*p(p_highcut) *gfftSize/2);
	if(hiCut<p_ports[p_highcut].min) {
		hiCut=p_ports[p_highcut].min;
	}
	if(hiCut>p_ports[p_highcut].max) {
		hiCut=p_ports[p_highcut].max;
	}
	CLAMP_PORT(long,windowsize,p_windowsize)
	unsigned long newwin=findBestFFTSize((unsigned long)windowsize);
	bool resetbuffers=false;
	if(newwin!=gfftSize) {
		gfftSize=newwin;
		FreeOldBuffers();
		AllocateNewBuffers(newwin);
		resetbuffers=true;
	}
	unsigned int iOsamp = 4*((unsigned int)(*p(p_overlapf)));

	iOsamp=calcOsampFromFFTSize(iOsamp,gfftSize);
	if(iOsamp!=goverlap) {
		goverlap=iOsamp;
		resetbuffers=true;
	}
	if(resetbuffers) {
		clearBuffers();
	}
	//Latency works out to be exactly the fft size.
	*p(p_latency)=gfftSize;
	CLAMP_PORT(int,iBlur,p_blur)
	CLAMP_PORT(float,fDecay,p_decay)
	bool consider_phase = *p(p_phase)>0;
	// arguments are number of samples to process, fft window size, sample overlap (4-32), input buffer, output buffer, init flag, gain, R input gain, decay, l cut, hi cut

	do_rebuild(sampleFrames, gfftSize, goverlap, sampleRate, p(p_left), p(p_right), p(p_outl), p(p_outr), fDecay, iBlur, loCut, hiCut, centre, consider_phase);
}
#define DEFINE_CIRC_COPY(NAME, OUTBUFFER, INBUFFER) \
	static inline int NAME(int circularsize,\
	float *  flat,\
	float *  circular,\
	int start,	int length,\
	float* __restrict window) {\
		int leftover=start+length-circularsize;\
		if(leftover>0) {\
			for(int circularindex=start; circularindex<circularsize; circularindex++) {\
				int flatindex=circularindex-start;\
				OUTBUFFER=INBUFFER;\
			}\
			for(int circularindex=0; circularindex<leftover; circularindex++) {\
				int flatindex=circularindex+(circularsize-start);\
				OUTBUFFER=INBUFFER;\
			}\
			return leftover;\
		} else {\
			for(int flatindex=0; flatindex<length; flatindex++) {\
				int circularindex=flatindex+start;\
				OUTBUFFER=INBUFFER;\
			}\
			return length+start;\
		}\
}
DEFINE_CIRC_COPY(copy_to_circular_buffer,circular[circularindex],flat[flatindex])
DEFINE_CIRC_COPY(copy_from_circular_buffer,flat[flatindex],circular[circularindex])

DEFINE_CIRC_COPY(copy_to_circular_buffer_extractcentre,circular[circularindex],flat[flatindex]-window[flatindex])
DEFINE_CIRC_COPY(copy_from_circular_buffer_window,flat[flatindex],circular[circularindex]*window[flatindex])

#define COPY_TO_FIFOS(amount)\
if(centreExtract>0) {\
	copy_to_circular_buffer_extractcentre(fftFrameSize, indata2, tInFIFO2, gRover, amount,indata);\
} else {\
	copy_to_circular_buffer(fftFrameSize, indata2, tInFIFO2, gRover, amount,NULL);\
}\
	gRover=copy_to_circular_buffer(fftFrameSize, indata, tInFIFO, gRover, amount,NULL);
	
// -----------------------------------------------------------------------------------------------------------------

float AKnockout::phaseToFrequency(float phase, int k,float dOversampbytwopi, float expct) {
	phase-=k*expct;
	phase+=PI;
	//now bring into range [0, 2*pi)
	long qpd = phase/(2*PI);
	phase -= 2*PI*qpd;
	phase+=2*PI*(phase<0);
	//return to [-pi,pi)
	phase-=PI;

	return ((double)k + phase*dOversampbytwopi);
}

inline void AKnockout::knockout(float origmag,float origphase, fftwf_complex * __restrict workspace, float* __restrict magnbuff, float* __restrict decaybuff, float coef, int k, float fDecayRate, int iBlur) {
	/* decay control */
	if (magnbuff[k]>decaybuff[k]) {
		decaybuff[k]=magnbuff[k];
	} else {
		decaybuff[k]=(decaybuff[k]-fDecayRate);
		decaybuff[k]=decaybuff[k]*(decaybuff[k]>1);
	}

	if (fDecayRate==0) {
		decaybuff[k]=magnbuff[k];/* if decay is off, set this value to right channel magn*/
	}
	/* spectral blur control */

	for (long m=-iBlur; m<iBlur; m++) {
		if (magnbuff[k+m]>decaybuff[k]) {
			decaybuff[k]=magnbuff[k+m];
		}
	}

	/* this is the 'knockout' process */

	origmag -= coef* decaybuff[k]; /* subtract right channel magnitudes from left, with decay*/
	origmag*= (origmag>0); /* zero -ve partials*/
	/* get real and imag part and re-interleave */
	myQT.QuickSinCos(origphase,workspace[k],workspace[k]+1);
	workspace[k][0] *=origmag;
	workspace[k][1] *=origmag;
}

inline void AKnockout::sumIntoCircularBuffer(float* __restrict outaccum,float dOutfactor,long outAccumIndex,long stepSize,long fftFrameSize){
	long inindex=0;
	/* do windowing and add to output accumulator */ 

	/*
	 *	For this step, we observe that tOutputAccum is a circular buffer
	 *  of size fftFrameSize.  The previous stepSize frames of data we
	 *  can discard, since it was just written out.  However, the remaining
	 *  Data we want to accumulate to.  Thus, we add to the next fftFrameSize - stepSize
	 *  frames of data, and overwrite the last stepSize frames.  However,
	 *  This is a little more complicated because we have to wrap around the circular buffers.
	 */
	long lastind=outAccumIndex-stepSize;
	if(lastind<0) {
		lastind+=fftFrameSize;
		//Fill up the part which we are still adding to.
		for(long k=outAccumIndex; k < lastind; k++) {
			outaccum[k] += window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
		//start to overwrite the old data at the end of the buffer.
		for(long k=lastind; k<fftFrameSize; k++) {
			outaccum[k] = window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
		//finish overwriting the remaining data at the wraparound.
		for(long k=0; k<outAccumIndex; k++) {
			outaccum[k] = window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
	} else {
		//Start accumulating at the end of the buffer.
		for(long k=outAccumIndex; k < fftFrameSize; k++) {
			outaccum[k] += window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
		//continue accumulating at the beginning.
		for(long k=0; k < lastind; k++) {
			outaccum[k] += window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
		//overwrite the last bit of the data.
		for(long k=lastind; k<outAccumIndex; k++) {
			outaccum[k] = window[inindex]*FFTRealBuffer[inindex]*(dOutfactor);
			inindex++;
		}
	}
}

void AKnockout::do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp,
		float sampleRate, float *indata, float *indata2, float *outdata1, float* outdata2,
		float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract, bool consider_phase) {
	//Declare these new local variables, so the compiler knows none of them are aliased.
	float* __restrict tInFIFO=gInFIFO;
	float* __restrict tOutputAccum=gOutputAccum;
	float* __restrict tOutputAccum2=gOutputAccum2;
	float* __restrict tFFTRealBuffer=FFTRealBuffer;
	float* __restrict tAnaPhase1=gAnaPhase1;
	float* __restrict tAnaPhase2=gAnaPhase2;
	float* __restrict tAnaMagn=gAnaMagn;
	float* __restrict tInFIFO2=gInFIFO2;
	float* __restrict tAnaMagn2=gAnaMagn2;
	float* __restrict twindow=window;

	/* set up some handy variables */
	long fftFrameSize2 = fftFrameSize/2;
	long stepSize = fftFrameSize/osamp;
	//float dOversampbytwopi = osamp/(PI*2);	
	//float freqPerBin = sampleRate/fftFrameSize;
	float dOutfactor = 1/(fftFrameSize2*osamp);
	fDecayRate=(fDecayRate>0)*(4.00001-(fDecayRate*fDecayRate*4));

	//float expct = 2.*PI*(double)stepSize/(double)fftFrameSize;
	{
		int numpro=copiesremaining;
		if(numSampsToProcess<copiesremaining) {
			numpro=numSampsToProcess;
			copiesremaining-=numSampsToProcess;
		} else {
			copiesremaining=0;
		}
		copy_from_circular_buffer(fftFrameSize,outdata1,tOutputAccum,outAccumIndex,numpro,NULL);
		outAccumIndex=copy_from_circular_buffer(fftFrameSize,outdata2,tOutputAccum2,outAccumIndex,numpro,NULL);
		outdata1+=numpro;
		outdata2+=numpro;
	}

	while(numSampsToProcess>=samples_needed_in_buffer) {
		
		COPY_TO_FIFOS(samples_needed_in_buffer)
		
		copy_from_circular_buffer_window(fftFrameSize,tFFTRealBuffer,tInFIFO,gRover,fftFrameSize,twindow);

		/* do transform */
		fftwf_execute(forward_sp1);

		copy_from_circular_buffer_window(fftFrameSize,tFFTRealBuffer,tInFIFO2,gRover,fftFrameSize,twindow);

		fftwf_execute(forward_sp2);

		/* frequency analysis */
		for (long k = 0; k <= fftFrameSize2; k++) {

			/* de-interlace FFT buffer */
			float real = gFFTworksp[k][0];
			float imag = gFFTworksp[k][1];

			/* compute magnitude and phase */
			tAnaMagn[k] = 2.*sqrt(real*real + imag*imag);
			tAnaPhase1[k]=atan2(imag,real);
		}

		/* this is the processing section */

		// lo cut
		for (long k = 0; k <= loCut+iBlur; k++) {  
			gFFTworksp[k][0]=0;
			gFFTworksp[k][1]=0;
			gFFTworksp2[k][0]=0;
			gFFTworksp2[k][1]=0;
		}

		// hi cut
		for (long k = fftFrameSize2-HiCut-iBlur; k <= fftFrameSize2; k++) {  
			gFFTworksp[k][0]=0;
			gFFTworksp[k][1]=0;		
			gFFTworksp2[k][0]=0;
			gFFTworksp2[k][1]=0;		
		}

		/* get R input magnitudes */

		for (long k = loCut; k <= fftFrameSize2-HiCut; k++) {
			float real=gFFTworksp2[k][0];
			float imag=gFFTworksp2[k][1];
			tAnaMagn2[k]=(2.*sqrt(real*real+imag*imag));
			tAnaPhase2[k]=atan2(imag,real);
		}
		
		for (long k = loCut+iBlur; k <= fftFrameSize2-HiCut-iBlur; k++) {
			/*we subtract out the part of the right channel *parallel to* the left*/
			float diff=fabs(tAnaPhase1[k]-tAnaPhase2[k]);
			diff*=consider_phase;
			float coef=cosf(diff);
			if(coef<0) {
				coef=0;
			}
			knockout(gAnaMagn[k],gAnaPhase1[k],gFFTworksp, gAnaMagn2,gDecay,coef,k,fDecayRate,iBlur);
			if(!centreExtract) {
				knockout(gAnaMagn2[k],gAnaPhase2[k],gFFTworksp2, gAnaMagn,gDecay2,coef,k,fDecayRate,iBlur);
			}
		}
		/* do inverse transform */
		fftwf_execute(backward_sp1);
		sumIntoCircularBuffer(tOutputAccum,dOutfactor,outAccumIndex,stepSize,fftFrameSize);
		if(!centreExtract) {
			fftwf_execute(backward_sp2);
		}
		sumIntoCircularBuffer(tOutputAccum2,dOutfactor,outAccumIndex,stepSize,fftFrameSize);

		//We output as much of the buffer as we can now, and use copiesremaining to notify later run() calls to empty the
		//rest of the buffer when it can.
		long overflow=stepSize-numSampsToProcess;
		long numcopy=stepSize;
		if(overflow>0) {
			copiesremaining=overflow;
			numcopy=numSampsToProcess;
		}
		copy_from_circular_buffer(fftFrameSize,outdata1,tOutputAccum,outAccumIndex,numcopy,NULL);
		outAccumIndex=copy_from_circular_buffer(fftFrameSize,outdata2,tOutputAccum2,outAccumIndex,numcopy,NULL);
		outdata1+=numcopy;
		outdata2+=numcopy;
		
		numSampsToProcess-=samples_needed_in_buffer;
		indata+=samples_needed_in_buffer;
		indata2+=samples_needed_in_buffer;
		samples_needed_in_buffer=stepSize;
	}
	COPY_TO_FIFOS(numSampsToProcess)
	samples_needed_in_buffer-=numSampsToProcess;
}

// -----------------------------------------------------------------------------------------------------------------

/* make lookup table for window function

   The windowing function used is a raised cosine.  This has very nice
 overlapping properties for any multiple of four ratio of overlapping.  He
 actually windows twice, once before analysis, and once after.  Thus the actual
 window used is the square of this window, i.e. the original window is 
 0.5 - 0.5*cos, so the square is (0.5-0.5*cos)^2=0.25*(1+cos^2-2*cos).  Now, let
 us look at the case of the 1/4 overlapping window.  Shifting the window by 1/4,
 the cosine becomes sine.  Thus we get a total of 0.25*(1+cos^2-2*cos+1+sin^2-2*sin).
 The sin^2 and cos^2 sum to 1, leaving 0.25*(2-2*sin-2*cos).  However, we have 
 to consider the sum of four consecutive windows, as this is the period by which
 it will repeat.  Since we have figured out the sum of the first two windows, we
 can simply shift this sum by 1/2 a window length to get the third and fourth
 windows.  The sines and cosines are both flipped in sign when we shift them
 half a period, so the sum is just going to be the sum of the constant parts.
 That is, it will be 0.25*2 +0.25*2 = 2.  This will be the factor by which we
 need to scale the windows.  However, if we overlap by some factor which is a
 multiple of 4, i.e. N*4, we can see that we will have N times as many windows
 covering the same area, and thus we will need to scale down by a factor of N*2.

	The variables used for this information are:
	 
	osamp: The overlapping factor.  Has to be a multiple of 4.
	
	stepSize: The number of samples in the overlapping period.
	
	freqPerBin: When we take the DFT, each bin represents a pure sine wave of a
		certain frequency.  Given a bin n, the sine wave frequency will be
		n*freqPerBin Hertz.
	
*/

void AKnockout::makelookup(int fftFrameSize)
{
	for (int k=0; k<fftFrameSize; k++) {
		window[k] = -.5*cos(2*PI*(double)k/(double)fftFrameSize)+.5;
	}
}

static int _ = AKnockout::register_class(p_uri);
