// Copyright (c) 1999-2003 Robin Davies.
// Placed into the public domain. 
// All uses permitted without attribution.

#include "QuickTrig.h"
#include <math.h> 



void CQuickTrigConsts::Initialize()
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