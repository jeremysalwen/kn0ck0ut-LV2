CXXFLAGS+=-fPIC `pkg-config lv2-plugin fftw3f --cflags`
LDFLAGS +=-shared `pkg-config lv2-plugin fftw3f --libs`
all: libkn0ck0ut.so

libkn0ck0ut.so: kn0ck0ut6.o
	$(CXX) $(LDFLAGS) kn0ck0ut6.o -o libkn0ck0ut.so
clean:
	rm -f *.o *.so
.phony: clean all
