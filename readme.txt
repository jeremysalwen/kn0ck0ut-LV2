kn0ck0ut-LV2 version 1.0 by Jeremy Salwen. An LV2 port of kn0ck0ut (http://www.freewebs.com/st3pan0va/)

Kn0ck0ut-LV2 is an LV2 plugin to perform spectral subtraction.  It can be used 
to achieve a wide variety of effects, most notably removing or extracting the
center of a two channel audio file.  As Kn0ck0ut is only a plugin, you will
need a host for LV2 plugins in order to use it, such as Ardour, Qtractor, Igen,
lv2_jack_host, or lv2file.

It requires lv2-c++-tools, lv2core, and fftw3 to build.

to build, run make.  To install, try make install.

kn0ck0ut-LV2 is released under the terms of the GPL version 3. 
(this is done so with the permission of st3pan0va).

In addition to the features of the original Kn0ck0ut, Kn0ck0ut-LV2 features:

* Improved performance through use of FFTW and more efficent buffering code
* Completely variable FFT size and overlap amount controls.
* Experimental "Phase Compensation" option, which will perhaps preserve
  additional fidelity in certain cases.
* Restored Low-cut filter which was removed in later releases.

Below Reads the Original Readme

====================================================================

Kn0ck0ut v0.5

--++/* by St3pan0va */++--


Notes on the source codes:

the code is a mess. sorry. im no programmer.

any newbies out there - don't use this as an example. 

any VST wizards - plz feel free to develop (mend ;) 

so long as yr source + plugin are freely available. 

no commercial use. no use without acknowledgement.

written/compiled on codewarrior 4.04 for windows

.mcp file included is the codewarrior project file.

you will need vst SDK from steinberg website.

contact st3pan0va [at] netscape [dot] net.



----

Acknowledgements / Licence

Kn0ck0ut is freeware and beta. No warranty, no

support. The plugin may be redistributed but only for free, 

and only with this readme file.

Made using the VST SDK and adapts (with thanks) code from S.M.Sprenger

available under the Wide Open Licence at www.dspdimension.com

Uses the QuickTrig class by Robin Davies at www.musicdsp.com

Thanks to TimG, FondueMeltdown, brainy & the ppl of *GYBO*.

