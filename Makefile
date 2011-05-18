CXXFLAGS+=-fPIC -DPIC `pkg-config lv2-plugin fftw3f --cflags`
LDFLAGS +=-shared `pkg-config lv2-plugin fftw3f --libs`

BUNDLE = kn0ck0ut.lv2
INSTALL_DIR = $(DESTDIR)/lib/lv2

all: $(BUNDLE)

$(BUNDLE): manifest.ttl kn0ck0ut.ttl libkn0ck0ut.so
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl kn0ck0ut.ttl libkn0ck0ut.so $(BUNDLE)
	
kn0ck0ut6.o: kn0ck0ut6.cpp kn0ck0ut6.hpp kn0ck0ut.peg

install: $(BUNDLE)
	mkdir -p $(INSTALL_DIR)
	rm -rf $(INSTALL_DIR)/$(BUNDLE)
	cp -R $(BUNDLE) $(INSTALL_DIR)

uninstall:
	rm -rf $(INSTALL_DIR)/$(BUNDLE)
	
libkn0ck0ut.so: kn0ck0ut6.o
	$(CXX) $(LDFLAGS) kn0ck0ut6.o -o libkn0ck0ut.so
	
clean:
	rm rm -rf $(BUNDLE) *.so *.o
	
.phony: clean all install
