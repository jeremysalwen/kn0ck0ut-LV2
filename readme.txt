kn0ck0ut-LV2 version 1.0 by Jeremy Salwen (jeremysalwen@gmail.com). 

An LV2 port of kn0ck0ut (http://www.freewebs.com/st3pan0va/)

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

See oldreadme.txt for the original readme of Kn0ck0ut by st3pan0va (Note that it is included for historical purposes only: The development and licensing information are both out of date and not applicable).
