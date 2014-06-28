/****************************\
| Using libfftw with C++     |
| www.admindojo.com          |
\****************************/
#include <fftw3.h>
#include <Magick++.h>
#include <math.h>
#include <complex>
#include <iostream>
#include <string>
using namespace Magick;
using namespace std;

void swapQuadrants(int squareSize, PixelPacket* pixels)
{
	int half = floor(squareSize / (double)2);
	
	// swap quadrants diagonally
	for(int i = 0; i < half; i++) {
		for(int j = 0; j < half; j++) {
			int upper = j + (squareSize * i);
			int lower = upper + (squareSize * half) + half;
			
			PixelPacket cur0 = *(pixels + upper);
			*(pixels + upper) = *(pixels + lower);
			*(pixels + lower) = cur0;
			
			PixelPacket cur1 = *(pixels + upper + half);
			*(pixels + upper + half) = *(pixels + lower - half);
			*(pixels + lower - half) = cur1;
		}
	}
}

void fft(int squareSize, PixelPacket* pixels, PixelPacket* outMag, PixelPacket* outPhase)
{
	fftw_plan planR, planG, planB;
	fftw_complex *inR, *inG, *inB, *outR, *outG, *outB;
	
	// allocate input arrays
	inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	inG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	inB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	
	// allocate output arrays
	outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	outG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	outB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	
	// create plans
	planR = fftw_plan_dft_2d(squareSize, squareSize, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);
	planG = fftw_plan_dft_2d(squareSize, squareSize, inG, outG, FFTW_FORWARD, FFTW_ESTIMATE);
	planB = fftw_plan_dft_2d(squareSize, squareSize, inB, outB, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// assign values to real parts (values between 0 and MaxRGB)
	for(int i = 0; i < squareSize * squareSize; i++) {
		PixelPacket current = *(pixels + i);
		double red = current.red;
		double green = current.green;
		double blue = current.blue;
		
		// save as real numbers
		inR[i][0] = red;
		inG[i][0] = green;
		inB[i][0] = blue;
	}
	
	// perform FORWARD fft
	fftw_execute(planR);
	fftw_execute(planG);
	fftw_execute(planB);
	
	// transform imaginary number to phase and magnitude and save to output
	for(int i = 0; i < squareSize * squareSize; i++) {
		// normalize values
		double realR = outR[i][0] / (double)(squareSize * squareSize);
		double imagR = outR[i][1] / (double)(squareSize * squareSize);
				
		double realG = outG[i][0] / (double)(squareSize * squareSize);
		double imagG = outG[i][1] / (double)(squareSize * squareSize);
		
		double realB = outB[i][0] / (double)(squareSize * squareSize);
		double imagB = outB[i][1] / (double)(squareSize * squareSize);
		
		// magnitude
		double magR = sqrt((realR * realR) + (imagR * imagR));
		double magG = sqrt((realG * realG) + (imagG * imagG));
		double magB = sqrt((realB * realB) + (imagB * imagB));
		
		// write to output
		(*(outMag + i)).red = magR;
		(*(outMag + i)).green = magG;
		(*(outMag + i)).blue = magB;
		
		// std::complex for arg()
		complex<double> cR(realR, imagR);
		complex<double> cG(realG, imagG);
		complex<double> cB(realB, imagB);
		
		// phase
		double phaseR = arg(cR) + M_PI;
		double phaseG = arg(cG) + M_PI;
		double phaseB = arg(cB) + M_PI;
		
		// scale and write to output
		(*(outPhase + i)).red = (phaseR / (double)(2 * M_PI)) * MaxRGB;
		(*(outPhase + i)).green = (phaseG / (double)(2 * M_PI)) * MaxRGB;
		(*(outPhase + i)).blue = (phaseB / (double)(2 * M_PI)) * MaxRGB;
	}
	
	// move zero frequency to (squareSize/2, squareSize/2)
	swapQuadrants(squareSize, outMag);
	swapQuadrants(squareSize, outPhase);
	
	// free memory
	fftw_destroy_plan(planR);
	fftw_destroy_plan(planG);
	fftw_destroy_plan(planB);
	fftw_free(inR); fftw_free(outR);
	fftw_free(inG); fftw_free(outG);
	fftw_free(inB); fftw_free(outB);
}

void ift(int squareSize, PixelPacket* inMag, PixelPacket* inPhase, PixelPacket* outPixels)
{
	// move zero frequency back to corners
	swapQuadrants(squareSize, inMag);
	swapQuadrants(squareSize, inPhase);
	
	fftw_plan planR, planG, planB;
	fftw_complex *inR, *inG, *inB, *outR, *outG, *outB;
	
	// allocate input arrays
	inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	inG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	inB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	
	// allocate output arrays
	outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	outG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	outB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * squareSize * squareSize);
	
	// create plans
	planR = fftw_plan_dft_2d(squareSize, squareSize, inR, outR, FFTW_BACKWARD, FFTW_ESTIMATE);
	planG = fftw_plan_dft_2d(squareSize, squareSize, inG, outG, FFTW_BACKWARD, FFTW_ESTIMATE);
	planB = fftw_plan_dft_2d(squareSize, squareSize, inB, outB, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	// transform magnitude/phase to real/imaginary
	for(int i = 0; i < squareSize * squareSize; i++) {
		double magR = inMag[i].red;
		double phaseR = ((inPhase[i].red / (double)MaxRGB) * 2 * M_PI) - M_PI;
		inR[i][0] = (magR * cos(phaseR));
		inR[i][1] = (magR * sin(phaseR));
		
		double magB = inMag[i].blue;
		double phaseB = ((inPhase[i].blue / (double)MaxRGB) * 2 * M_PI) - M_PI;
		inB[i][0] = (magB * cos(phaseB));
		inB[i][1] = (magB * sin(phaseB));
		
		double magG = inMag[i].green;
		double phaseG = ((inPhase[i].green / (double)MaxRGB) * 2 * M_PI) - M_PI;
		inG[i][0] = (magG * cos(phaseG));
		inG[i][1] = (magG * sin(phaseG));
	}
	
	// perform ift
	fftw_execute(planR);
	fftw_execute(planG);
	fftw_execute(planB);
	
	// save real parts to output
	for(int i = 0; i < squareSize * squareSize; i++) {
		double magR = outR[i][0];
		double magG = outG[i][0];
		double magB = outB[i][0];
		
		// make sure it's capped at MaxRGB
		(*(outPixels + i)).red = magR > MaxRGB ? MaxRGB : magR;
		(*(outPixels + i)).green = magG > MaxRGB ? MaxRGB : magG;
		(*(outPixels + i)).blue = magB > MaxRGB ? MaxRGB : magB;
	}
	
	// free memory
	fftw_destroy_plan(planR);
	fftw_destroy_plan(planG);
	fftw_destroy_plan(planB);
	fftw_free(inR); fftw_free(outR);
	fftw_free(inG); fftw_free(outG);
	fftw_free(inB); fftw_free(outB);
}

int main(int argc, char** argv)
{
	if(argc < 4) {
		cout << "usage:\t" << argv[0] << " in outMag outPhase" << endl
				<< "\t" << argv[0] << " -ift inMag inPhase out" << endl;
		return 0;
	}
	
	string arg = argv[1];
	
	if(arg != "-ift") {
		// FORWARD FOURIER TRANSFORM
		Image img(argv[1]);
		
		// get the length of the longer side of the image
		int squareSize = img.columns() < img.rows() ? img.rows() : img.columns();
		
		// the geometry of our padded image		
		Geometry padded(squareSize, squareSize);
		padded.aspect(true);
		
		// make image square
		img.extent(padded);
		
		// create templates for magnitude and phase
		Image mag(Geometry(squareSize, squareSize), "Black");
		Image phase(Geometry(squareSize, squareSize), "Black");
	
		// get image pixels
		img.modifyImage();
		Pixels pixelCache(img);
		PixelPacket* pixels;
		pixels = pixelCache.get(0, 0, squareSize, squareSize);
		
		// get magnitude pixels
		mag.modifyImage();
		Pixels pixelCacheMag(mag);
		PixelPacket* pixelsMag;
		pixelsMag = pixelCacheMag.get(0, 0, squareSize, squareSize);
		
		// get phase pixels
		phase.modifyImage();
		Pixels pixelCachePhase(phase);
		PixelPacket* pixelsPhase;
		pixelsPhase = pixelCachePhase.get(0, 0, squareSize, squareSize);
		
		// perform fft
		fft(squareSize, pixels, pixelsMag, pixelsPhase);
		
		// write changes
		pixelCache.sync();
		pixelCacheMag.sync();
		pixelCachePhase.sync();
		
		// save files
		mag.write(argv[2]);
		phase.write(argv[3]);
	}
	else {
		// BACKWARD FOURIER TRANSFORM
		Image mag(argv[2]);
		Image phase(argv[3]);
		
		// get size
		int squareSize = mag.columns();
		Image img(Geometry(squareSize, squareSize), "Black");
		
		// get image pixels
		img.modifyImage();
		Pixels pixelCache(img);
		PixelPacket* pixels;
		pixels = pixelCache.get(0, 0, squareSize, squareSize);
		
		// get magnitude pixels
		mag.modifyImage();
		Pixels pixelCacheMag(mag);
		PixelPacket* pixelsMag;
		pixelsMag = pixelCacheMag.get(0, 0, squareSize, squareSize);
		
		// get phase pixels
		phase.modifyImage();
		Pixels pixelCachePhase(phase);
		PixelPacket* pixelsPhase;
		pixelsPhase = pixelCachePhase.get(0, 0, squareSize, squareSize);
		
		// perform ift
		ift(squareSize, pixelsMag, pixelsPhase, pixels);
		
		// write changes
		pixelCache.sync();
		
		// save file
		img.write(argv[4]);
	}

	return 0;
}
