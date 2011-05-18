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
	
	strcpy (name, "Init");
}


//-----------------------------------------------------------------------------
AKnockout::AKnockout(audioMasterCallback audioMaster) : AudioEffectX(audioMaster, 16, kNumParams)	// 1 program, 1 parameter only
{
	
	programs = new AKnockoutProgram[numPrograms];
	fLoCut = fHiCut = vu = 0;
	
	if (programs)
		setProgram (0);
	
	hasVu ();	
	setNumInputs(2);		// stereo in
	setNumOutputs(2);		// mono out but cant tell wavelab that
	setUniqueID('k0ut');	// identify
//	canMono();				// doesn't make sense to feed both inputs with the same signal
	canProcessReplacing();	// supports only replacing output
	strcpy(programName, "Default");	// default program name
	
	gInFIFO = new float [MAX_FRAME_LENGTH];
	gOutFIFO = new float [MAX_FRAME_LENGTH];
	gFFTworksp = new float [2*MAX_FRAME_LENGTH];
	gOutputAccum = new float [2*MAX_FRAME_LENGTH];
	gAnaFreq = new float [MAX_FRAME_LENGTH];
	gAnaMagn = new float [MAX_FRAME_LENGTH];
	gInFIFO2 = new float [MAX_FRAME_LENGTH];
	gFFTworksp2 = new float [2*MAX_FRAME_LENGTH];
	gAnaMagn2 = new float [MAX_FRAME_LENGTH];
	gDecay = new float [MAX_FRAME_LENGTH];
	window = new double [FFTWINDOW];
	
	makelookup(FFTWINDOW);
	
	suspend ();		// flush buffer
	
}

//-----------------------------------------------------------------------------------------
AKnockout::~AKnockout() // delete buffers in destructor
{
	delete gInFIFO;
	delete gOutFIFO;
	delete gFFTworksp;
	delete gOutputAccum;
	delete gAnaFreq;
	delete gAnaMagn;
	delete gAnaMagn2;
	delete gInFIFO2;
	delete gFFTworksp2;
	delete gDecay;
	delete window;
	
}

//-----------------------------------------------------------------------------------------
void AKnockout::setProgram (long program)
{
	AKnockoutProgram * ap = &programs[program];

	curProgram = program;
	setParameter (kCentre, ap->fCentre);
	setParameter (kIn, ap->fIn);	
	setParameter (kLoCut, ap->fLoCut);
	setParameter (kOut, ap->fOut);
	setParameter (kHiCut, ap->fHiCut);
	setParameter (kDecay, ap->fDecay);
	setParameter (kBlur, ap->fBlur);
	
}

//------------------------------------------------------------------------
void AKnockout::setDelay (float fdelay)
{	
}
//------------------------------------------------------------------------

void AKnockout::setProgramName(char *name)
{
	strcpy(programName, name);
}

//-----------------------------------------------------------------------------------------
void AKnockout::getProgramName (char *name)
{
	if (!strcmp (programs[curProgram].name, "Init"))
		sprintf (name, "%s %d", programs[curProgram].name, curProgram + 1);
	else
		strcpy (name, programs[curProgram].name);
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

void AKnockout::resume ()
{
	return;
}
//------------------------------------------------------------------------
float AKnockout::getVu ()
{
	float cvu = vu;
	
	vu = 0;
	return cvu;
}

//------------------------------------------------------------------------
void AKnockout::setParameter (long index, float value)
{
	AKnockoutProgram * ap = &programs[curProgram];

	switch (index)
	{
		case kCentre :       fCentre = ap->fCentre = value; break;
		case kIn :       fIn = ap->fIn = value; break;
		case kLoCut : fLoCut = ap->fLoCut = value; break;
		case kOut :      fOut = ap->fOut = value; break;
		case kHiCut : 	 fHiCut = ap->fHiCut = value; break;
		case kDecay :	 fDecay = ap->fDecay = value; break;
		case kBlur :	 fBlur = ap->fBlur = value; break;
		
	}
	if (editor)
		editor->postUpdate ();
}


//-----------------------------------------------------------------------------------------
float AKnockout::getParameter (long index)
{
	float v = 0;

	switch (index)
	{
		case kCentre :		 v = (fCentre>0.5); break;
		case kIn :		 v = fIn; break;
		case kLoCut : v = fLoCut; break;
		case kOut :      v = fOut; break;
		case kHiCut :    v = fHiCut; break;
		case kDecay :	 v = fDecay; break;
		case kBlur :	 v = fBlur; break;
		
	}
	return v;
}
//-----------------------------------------------------------------------------------------

void AKnockout::getParameterName(long index, char *label)
{
	switch (index)
	{
		case kCentre :    strcpy (label, "Extract centre"); break;
		case kIn :    strcpy (label, "R Input Gain"); break;
		case kLoCut : strcpy (label, "Lo freq cut"); break;
		case kOut :      strcpy (label, "Output Gain "); break;
		case kHiCut : strcpy (label, "Hi freq cut"); break;
		case kDecay : strcpy(label, "Decay"); break;
		case kBlur : strcpy(label, "Blur"); break;
	}
}

//-----------------------------------------------------------------------------------------
void AKnockout::getParameterDisplay(long index, char *text)
{
	switch (index)
	{
		case kCentre :   if (fCentre>0.5) strcpy (text, "   On   ");
						if (fCentre<=0.5) strcpy (text, "   Off   ");	break;
		case kIn :    dB2string (fIn*4, text); break;
		case kLoCut : float2string (fLoCut, text);	break;
		case kOut :      dB2string (fOut*4, text); break;
		case kHiCut : float2string(fHiCut, text); break;
		case kDecay : float2string(fDecay, text); break;
		case kBlur : float2string(int(fBlur*24), text); break;
	}
}
//-----------------------------------------------------------------------------------------
void AKnockout::getParameterLabel(long index, char *label)
{
	switch (index)
	{
		case kCentre :   strcpy (label, "        "); 	break;
		case kIn :		strcpy (label, "   dB   ");	break;		
		case kLoCut : strcpy (label, " threshold ");	break;
		case kOut :      strcpy (label, "   dB   ");	break;
		case kHiCut : strcpy (label, " threshold ");	break;
		case kDecay : strcpy(label, " amount "); break;
		case kBlur : strcpy(label, " amount "); break;
	}
}
//-----------------------------------------------------------------------------------------
void AKnockout::process(float **inputs, float **outputs, long sampleFrames)
{
    //no send effect
}

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
  
  	do_rebuild(sampleFrames, FFTWINDOW, iOsamp, SAMPLERATE, in1, in2, out1,1, fOut*4, fIn*4, fDecay, iBlur, loCut, hiCut, centre);
	
    vu=(*out1);
    
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


inline void AKnockout::smsFft(float *fftBuffer, long fftFrameSize, long sign)
/* 
	FFT routine, (C)1996 S.M.Sprenger. Sign = -1 is FFT, 1 is iFFT (inverse)
	Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
	time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
	and returns the cosine and sine parts in an interleaved manner, ie.
	fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
	must be a power of 2. It expects a complex input signal (see footnote 2),
	ie. when working with 'common' audio signals our input signal has to be
	passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
	of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
	float wr, wi, arg, *p1, *p2, temp;
	float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
	long i, bitm, j, le, le2, k;
	long loopCondition = (long) (log(fftFrameSize)/log(2.)); 		// only a 
	long twiceFftFS = fftFrameSize*2;								// few tiny

	for (i = 2; i < twiceFftFS-2; i += 2) {
		for (bitm = 2, j = 0; bitm < twiceFftFS; bitm <<= 1) {
			if (i & bitm) j++;
			j <<= 1;
		}
		if (i < j) {
			p1 = fftBuffer+i; p2 = fftBuffer+j;
			temp = *p1; *(p1++) = *p2;
			*(p2++) = temp; temp = *p1;
			*p1 = *p2; *p2 = temp;
		}
	}
	for (k = 0, le = 2; k < loopCondition; k++) { 					// changes to
		le <<= 1;
		le2 = le>>1;
		ur = 1.0;
		ui = 0.0;
		arg = PI / (le2>>1);
		wr = myQT.QuickCos(arg); 									// sm sprenger's
		wi = sign*myQT.QuickSin(arg);								// fft code
		for (j = 0; j < le2; j += 2) {
			p1r = fftBuffer+j; p1i = p1r+1;
			p2r = p1r+le2; p2i = p2r+1;
			for (i = j; i < twiceFftFS; i += le) {
				tr = *p2r * ur - *p2i * ui;
				ti = *p2r * ui + *p2i * ur;
				*p2r = *p1r - tr; *p2i = *p1i - ti;
				*p1r += tr; *p1i += ti;
				p1r += le; p1i += le;
				p2r += le; p2i += le;
			}
			tr = ur*wr - ui*wi;
			ui = ur*wi + ui*wr;
			ur = tr;
		}
	}
}
