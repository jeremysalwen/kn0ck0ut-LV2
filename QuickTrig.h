// Copyright (c) 1999-2003 Robin Davies.
// Placed into the public domain. 
// All uses permitted without attribution.

#define PI 3.141592653589793238462643383279502

struct TSinCos { 
    float msin, mcos;
};

enum { 
        kMSBits = 10,
        kLSBits = 10,
        kMsTableSize = 1 << kMSBits,
        kLsTableSize = 1 << kLSBits,
        kBits = kMSBits + kLSBits,
		kPowBits = 8,
    };

static TSinCos mMsBitsTable[kMsTableSize+1];
static TSinCos mLsBitsTable[kLsTableSize+1];

class CQuickTrigConsts {
public:
	
	enum { 
		kMaxValidIndex = (2 << 20),
	};

protected:
   

protected:
    static void Initialize();
    
};

class CQuickTrig
: public CQuickTrigConsts 
{
public:

    inline double QuickSinQ(int nIndex) const {
        // Based on the identity sin(u+v) = sinu cosv + cosu sinv
        TSinCos *pscu = mMsBitsTable +( (nIndex >> kLSBits) & (kMsTableSize-1));
        TSinCos *pscv = mLsBitsTable + ( (nIndex) & (kLsTableSize-1));
        return pscu->msin * pscv->mcos + pscu->mcos * pscv->msin;
    };
    inline double QuickSin(double dAngle) const // Returns sin with 20 bits of precision.
    {
        return QuickSinQ((int)(dAngle*(kMsTableSize*kLsTableSize/(2*PI))) );
    }
    inline double QuickCosQ(int nIndex) const {
        // based on the identity cos(u+v) = cosu cosv + sinu sinv
        TSinCos *pscu = mMsBitsTable +( (nIndex >> kLSBits) & (kMsTableSize-1));
        TSinCos *pscv = mLsBitsTable + ( (nIndex) & (kLsTableSize-1));
        return pscu->mcos * pscv->mcos - pscu->msin * pscv->msin;
    };

    inline double QuickCos(double dAngle) const // Returns cos with 20 bits of precision.
    {
        return QuickCosQ((int)(dAngle*(kMsTableSize*kLsTableSize/(2*PI))) );
    }
	inline void QuickSinCos(double dAngle,float * cos, float* sin) const // Returns cos with 20 bits of precision.
    {
        return QuickSinCosQ((int)(dAngle*(kMsTableSize*kLsTableSize/(2*PI))),cos,sin );
    }

	inline void QuickSinCosQ(int nIndex, float* cos, float*sin) const {
	 // Based on the identity sin(u+v) = sinu cosv + cosu sinv
	    // based on the identity cos(u+v) = cosu cosv + sinu sinv
        TSinCos *pscu = mMsBitsTable +( (nIndex >> kLSBits) & (kMsTableSize-1));
        TSinCos *pscv = mLsBitsTable + ( (nIndex) & (kLsTableSize-1));
        *cos= pscu->mcos * pscv->mcos - pscu->msin * pscv->msin;
        *sin= pscu->msin * pscv->mcos + pscu->mcos * pscv->msin;
	}


};
