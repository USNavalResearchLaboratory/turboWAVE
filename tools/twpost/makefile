# makefile for turboWAVE post processing command line tool

PLATFORM = OSX

ALL_HEADERS = definitions.h fft.h twpost.h

ifeq ($(PLATFORM),OSX)
	CCC = g++
	CCFLAGS = -c -O3 
	LIB_PATH = 
	LIBS = -lstdc++
endif

tw3d: twpost.o FFT.o
	$(CCC) -o twpost twpost.o FFT.o $(LIB_PATH) $(LIBS)

twpost.o: twpost.cpp $(ALL_HEADERS)
	$(CCC) $(CCFLAGS) twpost.cpp

FFT.o: FFT.cpp $(ALL_HEADERS)
	$(CCC) $(CCFLAGS) FFT.cpp
