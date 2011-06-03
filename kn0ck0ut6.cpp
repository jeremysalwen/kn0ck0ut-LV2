/* kn0ck0ut v0.4 by st3pan0va 2004 */

#include "kn0ck0ut6.hpp"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FFTWINDOW 8192

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
	gInFIFO = new float [MAX_FRAME_LENGTH];
	gOutBuffer = new float [MAX_FRAME_LENGTH];
	FFTRealBuffer=new float[MAX_FRAME_LENGTH];
	gFFTworksp = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MAX_FRAME_LENGTH);
	gOutputAccum = new float [2*MAX_FRAME_LENGTH];
	gAnaFreq = new float [MAX_FRAME_LENGTH];
	gAnaMagn = new float [MAX_FRAME_LENGTH];
	gInFIFO2 = new float [MAX_FRAME_LENGTH];
	gFFTworksp2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MAX_FRAME_LENGTH);
	gAnaMagn2 = new float [MAX_FRAME_LENGTH];
	gDecay = new float [MAX_FRAME_LENGTH];
	memset(gDecay,0,MAX_FRAME_LENGTH*sizeof(float));
	window = new double [FFTWINDOW];

	forward_sp1= fftwf_plan_dft_r2c_1d(FFTWINDOW, FFTRealBuffer , gFFTworksp,
                                    FFTW_ESTIMATE);
	forward_sp2= fftwf_plan_dft_r2c_1d(FFTWINDOW,FFTRealBuffer, gFFTworksp2,
                                    FFTW_ESTIMATE);	
	backwards=fftwf_plan_dft_c2r_1d(FFTWINDOW, gFFTworksp, FFTRealBuffer,
                                    FFTW_ESTIMATE);

	makelookup(FFTWINDOW);

	suspend ();		// flush buffer

}

//-----------------------------------------------------------------------------------------
AKnockout::~AKnockout() // delete buffers in destructor
{
	delete[] gInFIFO;
	delete[] FFTRealBuffer;
	delete[] gOutBuffer;
	fftwf_free(gFFTworksp);
	delete[] gOutputAccum;
	delete[] gAnaFreq;
	delete[] gAnaMagn;
	delete[] gAnaMagn2;
	delete[] gInFIFO2;
	fftwf_free(gFFTworksp2);
	delete[] gDecay;
	delete[] window;

}


//-----------------------------------------------------------------------------------------
void AKnockout::suspend ()
{
	memset(gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(float));
	memset(gOutBuffer, 0, MAX_FRAME_LENGTH*sizeof(float));
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

	int loCut = int (*p(p_lowcut)); 
	int hiCut = int (*p(p_highcut) *FFTWINDOW/2);
	int centre = (*p(p_mode)>0);

	int iOsamp = 8;

	int iBlur = int (*p(p_blur));
	float fDecay= *p(p_decay);
	// arguments are number of samples to process, fft window size, sample overlap (4-32), input buffer, output buffer, init flag, gain, R input gain, decay, l cut, hi cut

	do_rebuild(sampleFrames, FFTWINDOW, iOsamp, sampleRate, p(p_left), p(p_right), p(p_out),1, fDecay, iBlur, loCut, hiCut, centre);
}

// -----------------------------------------------------------------------------------------------------------------


void AKnockout::do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, long gInit, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract)

{
	static long gRover=false;

	/* set up some handy variables */
	long fftFrameSize2 = fftFrameSize/2;
	long stepSize = fftFrameSize/osamp;
	double dOversampbytwopi = (double)osamp/PI/2;	
	double freqPerBin = sampleRate/(double)fftFrameSize;
	double dFreqfactor = PI/(double)osamp/freqPerBin*2;
	double dOutfactor = (double)fftFrameSize2*(double)osamp;
	fDecayRate=(fDecayRate>0)*(4.00001-(fDecayRate*fDecayRate*4));

	double expct = 2.*PI*(double)stepSize/(double)fftFrameSize;
	long inFifoLatency = fftFrameSize-stepSize;
	if (gRover == false) {
		gRover = inFifoLatency;
	}
	/* main processing loop */
	for (long i = 0; i < numSampsToProcess; i++){

		/* As long as we have not yet collected enough data just read in */
		gInFIFO[gRover] = indata[i];
		gInFIFO2[gRover] = (indata2[i]-(indata[i]*centreExtract)); // R input gain factor here
		outdata[i] = gOutBuffer[gRover-inFifoLatency]; // output, gain control here

		gRover++;

		/* now we have enough data for processing */
		if (gRover >= fftFrameSize) {
			gRover = inFifoLatency;


			/* do windowing  */
			for (long k = 0; k < fftFrameSize;k++) {
				FFTRealBuffer[k]=gInFIFO[k] * window[k];
			}

			/* do transform */
			fftwf_execute(forward_sp1);
			
			for (long k = 0; k < fftFrameSize;k++) {
				FFTRealBuffer[k]=gInFIFO2[k] * window[k];
			}
			
			fftwf_execute(forward_sp2);

			/* frequency analysis */
			for (long k = 0; k <= fftFrameSize2; k++) {

				/* de-interlace FFT buffer */
				double real = gFFTworksp[k][0];
				double imag = gFFTworksp[k][1];

				/* compute magnitude and phase */
				gAnaMagn[k] = 2.*sqrt(real*real + imag*imag);
				double phase = atan2(imag,real);

				double tmp = phase-(double)k*expct;
				long qpd = tmp/PI;
				if (qpd >= 0) qpd += qpd&1;
				else qpd -= qpd&1;
				tmp -= PI*(double)qpd;
				tmp *= dOversampbytwopi;

				/* store frequency in analysis array */

				gAnaFreq[k] = ((double)k + tmp)*freqPerBin;

			}

			/* this is the processing section */

			// lo cut
			for (long k = 0; k <= loCut+iBlur; k++) {  
				gFFTworksp[k][0]=0;
				gFFTworksp[k][1]=0;
			}

			// hi cut
			for (long k = fftFrameSize2-HiCut-iBlur; k <= fftFrameSize2; k++) {  
				gFFTworksp[k][0]=0;
				gFFTworksp[k][1]=0;		
			}

			/* get R input magnitudes */

			for (long k = loCut; k <= fftFrameSize2-HiCut; k++) {
			    float real=gFFTworksp2[k][0];
			    float imag=gFFTworksp2[k][1];
				gAnaMagn2[k]=(2.*sqrt(real*real+imag*imag));
			}


			for (long k = loCut+iBlur; k <= fftFrameSize2-HiCut-iBlur; k++) {


				/* decay control */

				if (gAnaMagn2[k]>gDecay[k]) {
					gDecay[k]=gAnaMagn2[k];
				} else {
					gDecay[k]=(gDecay[k]-fDecayRate);
					gDecay[k]=gDecay[k]*(gDecay[k]>1);
				}
				if (fDecayRate==0) {
					gDecay[k]=gAnaMagn2[k]; // if decay is off, set this value to right channel magn
				}
				/* spectral blur control */

				for (long m=-iBlur; m<iBlur; m++) {
					if (gAnaMagn2[k+m]>gDecay[k]) {
						gDecay[k]=gAnaMagn2[k+m];
					}
				}					   

				/* this is the 'knockout' process */

				double magn = gAnaMagn[k] - gDecay[k]; // subtract right channel magnitudes from left, with decay
				magn = magn * (magn>0); // zero -ve partials

				//(Note by Jeremy): The phase has not been modified at all.
				//This exactly undoes the transformation of the phase of the left
				//channel which we did during the analysis phase.
				
				/* correct the frequency - sm sprenger method */
				double tmp = gAnaFreq[k];
				tmp -= (double)k*freqPerBin;
				tmp *= dFreqfactor;
				tmp += (double)k*expct;

				/* get real and imag part and re-interleave */
				gFFTworksp[k][0] = magn*myQT.QuickCos(tmp);
				gFFTworksp[k][1] = magn*myQT.QuickSin(tmp);

			} 

			/* do inverse transform */
			fftwf_execute(backwards);

			/* do windowing and add to output accumulator */ 
			for(long k=0; k < fftFrameSize; k++) {
				gOutputAccum[k] += window[k]*FFTRealBuffer[k]/(dOutfactor);
			}

			/* transfer output accum to output buffer */
			memmove (gOutBuffer, gOutputAccum, stepSize*sizeof(float));

			/* shift accumulator */
			memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(float));

			memmove (gInFIFO, gInFIFO+stepSize, inFifoLatency*sizeof(float));
			memmove (gInFIFO2, gInFIFO2+stepSize, inFifoLatency*sizeof(float));

		}
	}
}

// -----------------------------------------------------------------------------------------------------------------

/* make lookup table for window function

   The windowing function used is a raised cosine.  This has very nice
 overlapping properties for any even ratio of overlapping.  He actually windows
 twice, once before analysis, and once after.  Thus the actual window used is
 the square of this window, i.e. the original window is 0.5 - 0.5*cos, so the
 square is (0.5-0.5*cos)^2=0.25*(1+cos^2-2*cos).  Now, let us look at the case
 of the 1/4 overlapping window.  Shifting the window by 1/4, the cosine becomes
 sine.  Thus we get a total of 0.25*(1+cos^2-2*cos+1+sin^2-2*sin).  The sin^2
 and cos^2 sum to 1, leaving 0.25*(2-2*sin-2*cos).  However, we have to consider
 the sum of four consecutive windows, as this is the period by which it will
 repeat.  Since we have figured out the sum of the first two windows, we can
 simply shift this sum by 1/2 a window length to get the third and fourth
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
