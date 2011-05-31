/* kn0ck0ut v0.4 by st3pan0va 2004 */

#include "kn0ck0ut6.hpp"
#include "AEffEditor.hpp"
#include <stdio.h> 
#include <string.h> 
#include <math.h> 

#define FFTWINDOW 8192
#define SAMPLERATE 44100 

//-----------------------------------------------------------------------------------------
void AKnockout::processReplacing(float **inputs, float **outputs, long sampleFrames)
{
    float *in1  =  inputs[0];
    float *in2  =  inputs[1];   
    float *out1 = outputs[0];
  	
  	int loCut = int (fLoCut*128); 
  	int hiCut = int (fHiCut*FFTWINDOW/2);
  	int centre = (fCentre>0.5);
  	
  	int iOsamp = 8;
  	
  	int iBlur = int (fBlur*24);
  	
    // arguments are number of samples to process, fft window size, sample overlap (4-32), input buffer, output buffer, init flag, gain, R input gain, decay, l cut, hi cut
  
  	do_rebuild(sampleFrames, FFTWINDOW, iOsamp, SAMPLERATE, in1, in2, out1, fOut*4, fIn*4, fDecay, iBlur, loCut, hiCut, centre);
	
    vu=(*out1);
    
}

// -----------------------------------------------------------------------------------------------------------------

	
void AKnockout::do_rebuild(long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *indata2, float *outdata, float fGain, float fInGain, float fDecayRate, int iBlur, int loCut, int HiCut, int centreExtract)

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
				gFFTworksp[2*k] = gInFIFO[k] * window[k];
				gFFTworksp2[2*k] = gInFIFO2[k] * window[k];
				gFFTworksp[2*k+1] = gFFTworksp2[2*k+1] = 0.;
				}
				
			/* do transform */
			smsFft(gFFTworksp, fftFrameSize, -1);

			/* frequency analysis */
			for (k = 0; k <= fftFrameSize2; k++) {

				/* de-interlace FFT buffer */
				real = gFFTworksp[2*k];
				imag = gFFTworksp[2*k+1];

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
			smsFft(gFFTworksp2, fftFrameSize, -1);


			/* this is the processing section */
			
			// lo cut
			for (k = 0; k <= loCut+iBlur; k++) {  
				gFFTworksp[2*k]=0;
				gFFTworksp[2*k+1]=0;
				}
				
			// hi cut
			for (k = fftFrameSize2-HiCut-iBlur; k <= fftFrameSize2; k++) {  
				gFFTworksp[2*k]=0;
				gFFTworksp[2*k+1]=0;		
				}
				
				/* get R input magnitudes */
				
			for (k = loCut; k <= fftFrameSize2-HiCut; k++) {
					gAnaMagn2[k]=(2.*sqrt(gFFTworksp2[2*k]*gFFTworksp2[2*k] + gFFTworksp2[2*k+1]*gFFTworksp2[2*k+1]));
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
				gFFTworksp[2*k] = magn*myQT.QuickCos(tmp);
				gFFTworksp[2*k+1] = magn*myQT.QuickSin(tmp);
	
				} 


			/* zero negative frequencies */
			for (k = fftFrameSize+2; k < 2*fftFrameSize; k++) gFFTworksp[k] = 0.;
			
			/* do inverse transform */
			smsFft(gFFTworksp, fftFrameSize, 1);

			/* do windowing and add to output accumulator */ 
			for(k=0; k < fftFrameSize; k++) {
				gOutputAccum[k] += window[k]*gFFTworksp[2*k]/(dOutfactor);
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

for (int k=0; k<fftFrameSize; k++) window[k] = -.5*cos(2*PI*(double)k/(double)fftFrameSize)+.5;

}
// -----------------------------------------------------------------------------------------------------------------