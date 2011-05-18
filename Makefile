CXXFLAGS+=`pkg-config lv2-plugin fftw3f --cflags`
LDFLAGS +=`pkg-config lv2-plugin fftw3f --libs`
all: libkn0ck0ut.so

libkn0ck0ut.so: kn0ck0ut6.o QuickTrig.o
