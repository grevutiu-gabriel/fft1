CXXFLAGS=`Magick++-config --cppflags --cxxflags --ldflags --libs` `pkg-config fftw3 --libs` -g -Wall -O0

CXX=g++

all: fourier

fourier: main.o
	$(CXX) -o $@ $< $(CXXFLAGS) 

main.o: main.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<
	
clean:
	rm main.o

install:
	cp fourier /usr/bin/

uninstall:
	rm /usr/bin/fourier
