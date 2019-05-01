#ifndef IMAGE_H
#define IMAGE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "channel.h"

class Image{
public:
	Channel *R, *G, *B;
	int *filter;
	int width, height, range, type; // atributos da imagem
	std::string path, csv_path;

	// Constructors
	Image();
	Image(std::string path);
	~Image();

	// IO
	void loadFile(std::string path);
	void reloadFile();
	void saveFile();
	void saveFile(std::string path);
	void loadCSV(std::string path);
	void reloadCSV();

	// Modifications
	void negative();
	void mirror();
	void plus(int add);
	void plusNoLimit(int add);
	void minus(int sub);
	void threshold(int threshold); // REVER
	void maximizer();

	// Filters
	void applyFilterNoLimit();
	void applyFilter();
	void applyFilterNoLimit(std::string csvfile);
	void applyFilter(std::string csvfile);
	void applyKernel(unsigned int index);
	void applyNoLinear(float div);

	// Distances
	void L1(int refR, int refG, int refB);
	void L2(int refR, int refG, int refB);
	void Mahalanobis(std::string path);
	void Knn(std::string path, int K);

	// Edge Detection
	void sobel();
	void roberts();
	void robinson();

	// Fill
	void fill(int threshold);
	void flood(int index, unsigned int *groups, unsigned int &number, int &th);
	void fill2(int threshold);
	void flood2(int index, unsigned int *groups, unsigned int &number, int &th);
	void fill3(int threshold);
	void flood3(int index, unsigned int *groups, unsigned int &number, int &th);
	bool Try(int index, unsigned int *groups, int &th, int dir);

	// Others
	Channel *grayscale();
	// Mahalanobis subsection
	// Knn subsection

	//debug
	void showData(int a);
};

#endif // TEST_H
