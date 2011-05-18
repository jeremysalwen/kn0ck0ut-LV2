/* kn0ck0ut v0.4 by st3pan0va 2004 */

#include "kn0ck0ut6.hpp"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FFTWINDOW 8192
#define SAMPLERATE 44100

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

AKnockoutProgram::AKnockoutProgram ()
{
	fCentre = 0;
	fIn=0.25;
	fLoCut = 0;
	fOut = 0.25;
	fHiCut = 0;
	fDecay = 0;

}


//-----------------------------------------------------------------------------
AKnockout::AKnockout(double rate) : Plugin<AKnockout>(8)
{
		fLoCut = fHiCut = 0;


	gInFIFO = new float [MAX_FRAME_LENGTH];
	gOutFIFO = new float [MAX_FRAME_LENGTH];
	gFFTworksp = (fftwf_complex*)fftw_malloc(sizeof(fftw_complex) * MAX_FRAME_LENGTH);
	gOutputBuffer=new float [2*MAX_FRAME_LENGTH];
	gOutputAccum = new float [2*MAX_FRAME_LENGTH];
	gAnaFreq = new float [MAX_FRAME_LENGTH];
	gAnaMagn = new float [MAX_FRAME_LENGTH];
	gInFIFO2 = new float [MAX_FRAME_LENGTH];
	gFFTworksp2 = (fftwf_complex*)fftw_malloc(sizeof(fftw_complex) * MAX_FRAME_LENGTH);
	gAnaMagn2 = new float [MAX_FRAME_LENGTH];
	gDecay = new float [MAX_FRAME_LENGTH];
	window = new double [FFTWINDOW];

	forward_sp1= fftwf_plan_dft_r2c_1d(FFTWINDOW, gInFIFO , gFFTworksp,
                                    FFTW_ESTIMATE);
	forward_sp2= fftwf_plan_dft_r2c_1d(FFTWINDOW, gInFIFO2, gFFTworksp2,
                                    FFTW_ESTIMATE);	
	backwards=fftwf_plan_dft_c2r_1d(FFTWINDOW, gFFTworksp, gOutputBuffer,
                                    FFTW_ESTIMATE);

	makelookup(FFTWINDOW);

	suspend ();		// flush buffer

}

//-----------------------------------------------------------------------------------------
AKnockout::~AKnockout() // delete buffers in destructor
{
	delete gInFIFO;
	delete gOutFIFO;
	fftw_free(gFFTworksp);
	delete gOutputAccum;
	delete gOutputBuffer;
	delete gAnaFreq;
	delete gAnaMagn;
	delete gAnaMagn2;
	delete gInFIFO2;
	fftw_free(gFFTworksp2);
	delete gDecay;
	delete window;

}


//-----------------------------------------------------------------------------------------
void AKnockout::suspend ()
{
	memset(gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gOutFIFO, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gFFTworksp, 0, 2*MAX_FRAME_LENGTH*sizeof(float));
	memset(gOutputAccum, 0, 2*MAX_FRAME_LENGTH*sizeof(float));
	memset(gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gInFIFO2, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gFFTworksp2, 0, 2*MAX_FRAME_LENGTH*sizeof(float));
	memset(gAnaMagn2, 0, MAX_FRAME_LENGTH*sizeof(float));
}

//-----------------------------------------------------------------------------------------
void AKnockout::run(uint32_t sampleFrames)
{

	int loCut = int (fLoCut*128); 
	int hiCut = int (fHiCut*FFTWINDOW/2);
	int centre = (fCentre>0.5);

	int iOsamp = 8;

	int iBlur = int (fBlur*24);

	// arguments are number of samples to process, fft window size, sample overlap (4-32), input buffer, output buffer, init flag, gain, R input gain, decay, l cut, hi cut

	do_rebuild(sampleFrames, FFTWINDOW, iOsamp, SAMPLERATE, p(0), p(1), p(2),1, fOut*4, fIn*4, fDecay, iBlur, loCut, hiCut, centre);
}

// -----------------------------------------------------------------------------------------------------------------


void AKnockout::do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, long gInit, float fGain, float fInGain, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract)

{

	static long gRover=false;

	double magn, phase, tmp, real, imag;
	double freqPerBin, expct;
	long i,k, qpd, inFifoLatency, stepSize, fftFrameSize2, m;
	double dOversampbytwopi, dFreqfactor, dOutfactor;

	/* set up some handy variables */
	fftFrameSize2 = fftFrameSize/2;
	stepSize = fftFrameSize/osamp;
	dOversampbytwopi = (double)osamp/PI/2;	
	freqPerBin = sampleRate/(double)fftFrameSize;
	dFreqfactor = PI/(double)osamp/freqPerBin*2;
	dOutfactor = (double)fftFrameSize2*(double)osamp;
	fDecayRate=(fDecayRate>0)*(4.00001-(fDecayRate*fDecayRate*4));

	expct = 2.*PI*(double)stepSize/(double)fftFrameSize;
	inFifoLatency = fftFrameSize-stepSize;
	if (gRover == false) gRover = inFifoLatency;

	/* main processing loop */
	for (i = 0; i < numSampsToProcess; i++){

		/* As long as we have not yet collected enough data just read in */
		gInFIFO[gRover] = indata[i];
		gInFIFO2[gRover] = (indata2[i]-(indata[i]*centreExtract)) * fInGain; // R input gain factor here
		outdata[i] = gOutFIFO[gRover-inFifoLatency] * fGain; // output, gain control here

		gRover++;

		/* now we have enough data for processing */
		if (gRover >= fftFrameSize) {
			gRover = inFifoLatency;


			/* do windowing and re,im interleave */
			for (k = 0; k < fftFrameSize;k++) {
				gInFIFO[k] *= window[k];
				gInFIFO2[k] *= window[k];
			}

			/* do transform */
			fftwf_execute(forward_sp1);

			/* frequency analysis */
			for (k = 0; k <= fftFrameSize2; k++) {

				/* de-interlace FFT buffer */
				real = gFFTworksp[2*k][0];
				imag = gFFTworksp[2*k][1];

				/* compute magnitude and phase */
				gAnaMagn[k] = 2.*sqrt(real*real + imag*imag);
				phase = atan2(imag,real);

				tmp = phase-(double)k*expct;
				qpd = tmp/PI;
				if (qpd >= 0) qpd += qpd&1;
				else qpd -= qpd&1;
				tmp -= PI*(double)qpd;
				tmp *= dOversampbytwopi;

				/* store frequency in analysis array */

				gAnaFreq[k] = ((double)k + tmp)*freqPerBin;

			}

			/* do transform */
			fftwf_execute(forward_sp2);


			/* this is the processing section */

			// lo cut
			for (k = 0; k <= loCut+iBlur; k++) {  
				gFFTworksp[2*k][0]=0;
				gFFTworksp[2*k][1]=0;
			}

			// hi cut
			for (k = fftFrameSize2-HiCut-iBlur; k <= fftFrameSize2; k++) {  
				gFFTworksp[2*k][0]=0;
				gFFTworksp[2*k][1]=0;		
			}

			/* get R input magnitudes */

			for (k = loCut; k <= fftFrameSize2-HiCut; k++) {
				gAnaMagn2[k]=(2.*sqrt(gFFTworksp2[2*k][0]*gFFTworksp2[2*k][0] + gFFTworksp2[2*k][1]*gFFTworksp2[2*k][1]));
			}


			for (k = loCut+iBlur; k <= fftFrameSize2-HiCut-iBlur; k++) {


				/* decay control */

				if (gAnaMagn2[k]>gDecay[k]) gDecay[k]=gAnaMagn2[k];
				else {
					gDecay[k]=(gDecay[k]-fDecayRate);
					gDecay[k]=gDecay[k]*(gDecay[k]>1);
				}
				if (fDecayRate==0) gDecay[k]=gAnaMagn2[k]; // if decay is off, set this value to right channel magn

				/* spectral blur control */

				for (m=-iBlur; m<iBlur; m++) {
					if (gAnaMagn2[k+m]>gDecay[k]) gDecay[k]=gAnaMagn2[k+m];
				}					   

				/* this is the 'knockout' process */

				magn = gAnaMagn[k] - gDecay[k]; // subtract right channel magnitudes from left, with decay
				magn = magn * (magn>0); // zero -ve partials

				/* correct the frequency - sm sprenger method */
				tmp = gAnaFreq[k];
				tmp -= (double)k*freqPerBin;
				tmp *= dFreqfactor;
				tmp += (double)k*expct;

				/* get real and imag part and re-interleave */
				gFFTworksp[2*k][0] = magn*myQT.QuickCos(tmp);
				gFFTworksp[2*k][1] = magn*myQT.QuickSin(tmp);

			} 

			/* do inverse transform */
			fftwf_execute(backwards);

			/* do windowing and add to output accumulator */ 
			for(k=0; k < fftFrameSize; k++) {
				gOutputAccum[k] += window[k]*gOutputBuffer[k]/(dOutfactor);
			}

			/* transfer output accum to output buffer */
			memmove (gOutFIFO, gOutputAccum, stepSize*sizeof(float));

			/* shift accumulator */
			memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(float));

			memmove (gInFIFO, gInFIFO+stepSize, inFifoLatency*sizeof(float));
			memmove (gInFIFO2, gInFIFO2+stepSize, inFifoLatency*sizeof(float));

		}
	}
}

// -----------------------------------------------------------------------------------------------------------------

/* make lookup table for window function */

void AKnockout::makelookup(int fftFrameSize)
{
	for (int k=0; k<fftFrameSize; k++) {
		window[k] = -.5*cos(2*PI*(double)k/(double)fftFrameSize)+.5;
	}
}