#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class Channel{
public:
	unsigned char *data, *copy;
	int *filter;
	int width, height, range, type, min, max; // atributos da imagem
	int dimension, offset, weight; // atributos do filtro carregado
	std::string path, csv_path;

	// Constructors
	Channel();
	Channel(std::string path);
	Channel(int o_width, int o_height);
	Channel(unsigned char *o_data, int o_width, int o_height);
	Channel(int *o_data, int o_width, int o_height);
	~Channel();

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
	void minusNoLimit(int add);
	void threshold(int threshold);
	void thresholdRange(int low, int high);
	void maximizer();

	// Filters
	void applyFilterNoLimit();
	void applyFilter();
	void applyFilterNoLimit(std::string csvfile);
	void applyFilter(std::string csvfile);
	int* returnFiltered();
	void applyKernel(unsigned int index);
	void applyNoLinear(float div);

	// Edge Detection
	void sobel_old();
	void sobel();
	void roberts();
	void robinson();

	// Others
	unsigned char *expandedCopy();
	void update();
	unsigned char *to255(int *in);
	unsigned char *int2char(int *in);

	// Debug
	void showData(int a);
};

#endif // CHANNEL_H
