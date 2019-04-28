

#include "image.h"

int main(int argc, char const *argv[]){
	
	/*
	Channel cha1, cha2;
	cha1.loadFile("../images/star.pgm");
	
	cha1.sobel();
	cha1.saveFile("1-sobel.pgm");
	cha1.reloadFile();

	cha1.roberts();
	cha1.saveFile("3-roberts.pgm");
	cha1.reloadFile();

	cha1.robinson();
	cha1.saveFile("5-robinson.pgm");
	*/

	Image img1;
	img1.loadFile("../images/gabarito.ppm");
	img1.applyFilter("../csv/gauss.csv");
	
	img1.sobel();
	img1.saveFile("2-sobel.ppm");
	img1.threshold(50);
	img1.saveFile("bin-2-sobel.pgm");
	img1.reloadFile();

	img1.roberts();
	img1.saveFile("4-roberts.ppm");
	img1.threshold(50);
	img1.saveFile("bin-4-roberts.pgm");
	img1.reloadFile();

	img1.robinson();
	img1.saveFile("6-robinson.ppm");
	img1.threshold(50);
	img1.saveFile("bin-6-robinson.pgm");
	img1.reloadFile();
	
	Channel *cha2;
	cha2 = img1.grayscale();

	cha2->sobel();
	cha2->saveFile("1-sobel.pgm");
	cha2->threshold(50);
	cha2->saveFile("bin-1-sobel.pgm");
	cha2 = img1.grayscale();

	cha2->roberts();
	cha2->saveFile("3-roberts.pgm");
	cha2->threshold(50);
	cha2->saveFile("bin-3-roberts.pgm");
	cha2 = img1.grayscale();

	cha2->robinson();
	cha2->saveFile("5-robinson.pgm");
	cha2->threshold(50);
	cha2->saveFile("bin-5-robinson.pgm");

	return 0;
}