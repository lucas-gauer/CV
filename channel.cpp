#include "channel.h"

void negativeFilter(int *in, int length);

Channel::Channel()
:data(nullptr),copy(nullptr),filter(nullptr),
width(0),height(0),range(255),type(0),min(255),max(0),
dimension(0),offset(0),weight(0)
{}

Channel::Channel(std::string path)
:copy(nullptr),filter(nullptr),
min(255),max(0),
dimension(0),offset(0),weight(0)
{ // load genérico
	loadFile(path);
}

Channel::Channel(int o_width, int o_height)
:copy(nullptr),filter(nullptr),
width(o_width),height(o_height),range(255),type(2),min(255),max(0),
dimension(0),offset(0),weight(0)
{ // para quando souber apenas o tamanho desejado	
	data = new unsigned char[width * height]{0};
}

Channel::Channel(unsigned char *o_data, int o_width, int o_height)
:data(o_data),copy(nullptr),filter(nullptr),
width(o_width),height(o_height),range(255),type(2),min(255),max(0),
dimension(0),offset(0),weight(0)
{} // para quando souber o tamanho e a informação

Channel::Channel(int *o_data, int o_width, int o_height)
:copy(nullptr),filter(nullptr),
width(o_width),height(o_height),range(255),type(2),min(255),max(0),
dimension(0),offset(0),weight(0)
{
	data = to255(o_data);
}

Channel::~Channel(){
	delete[] data;
	//delete[] copy;
	delete[] filter;
}

// -----------------------------------IO-----------------------------------

void Channel::loadFile(std::string o_path){ // le o arquivo e define os parâmetros do objeto
	std::ifstream inFile;
	std::string strbuff;

	path = std::string(o_path);

	inFile.open(path);
	if(!inFile){
		std::cout << "Unable to open file\n";
		exit(1);
	}

	inFile >> strbuff;
	if(strbuff == "P2"){
		type = 2;
	}
	else if(strbuff == "P3"){
		std::cout << "This class is for single channel images (PGM) :(\n";
		exit(1);
	}

	inFile >> width;
	inFile >> height;
	inFile >> range;

	delete[] data;
	data = new unsigned char[width * height];

	int i = 0;
	while(i < width * height){ // linhas de conteúdo
		inFile >> strbuff;
		data[i] = stoi(strbuff);
		if(data[i] > max){
			max = data[i];
		}
		if(data[i] < min){
			min = data[i];
		}
		i++;
	}

	inFile.close();
}

void Channel::reloadFile(){
	loadFile(path);
}

void Channel::saveFile(){ // salva o estado atual
	std::ofstream outFile;

	outFile.open("output.pgm");
	outFile << "P" << type << "\n";
	outFile << width << " " << height << " " << range << "\n";
	
	for(int k = 0; k < height; k++){
		for(int l = 0; l < width; l++){
			outFile << +data[k*width+l] << "\n";
		}
	}

	outFile.close();
}

void Channel::saveFile(std::string o_path){ // salva em um caminho específico
	std::ofstream outFile;

	outFile.open(o_path);
	outFile << "P" << type << "\n";
	outFile << width << " " << height << " 255\n";
	
	for(int k = 0; k < height; k++){
		for(int l = 0; l < width; l++){
			outFile << +data[k*width+l] << "\n";
		}
	}

	outFile.close();
}

void Channel::loadCSV(std::string o_path){ // le arquivos de 'matriz' para aplicar filtros
	csv_path = std::string(o_path);
	std::ifstream inFile;
	std::string strbuff;
	dimension = 1;
	int length = 0;
	char buffer[10];
	weight = 0;

	inFile.open(o_path);
	if(!inFile){
		std::cout << "Unable to open .csv file\n";
		exit(1);
	}

	inFile >> strbuff;
	length = strbuff.length();
	for(int i = 0; i < length; i++){
		if(strbuff[i] == ','){
			dimension++;
		}
	}
	offset = dimension/2;
	
	for(int i = 0; i < length; ++i){
		inFile.unget();
	}

	delete[] filter;
	filter = new int[dimension*dimension];

	int j;
	for(int i = 0; i < dimension; i++){
		for(j = 0; j < dimension-1; j++){
			inFile.get(buffer, 10, ',');
			filter[i*dimension+j] = stoi(std::string(buffer));
			weight += (filter[i*dimension+j]);
			inFile.get();
		}
		inFile.get(buffer, 10, '\n');
		filter[i*dimension+j] = stoi(std::string(buffer));
		weight += filter[i*dimension+j];
		inFile.get();
	}

	if(weight == 0) weight = 1;

	inFile.close();
}

void Channel::reloadCSV(){
	if(!csv_path.empty()){
		loadCSV(csv_path);
	}
}

// -----------------------------------MODIFICATIONS-----------------------------------

void Channel::negative(){ // inverte os valores
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			data[i*width+j] = 255 - data[i*width+j];
		}
	}
}

void Channel::mirror(){ // espelha a imagem em relação a um eixo vertical
	copy = new unsigned char[width*height];

	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			copy[i*width+j] = data[i*width+j];
		}
	}

	for(int i = 0; i < height; i++){
		for(int j = 1; j <= width; j++){
			data[(i+1)*width-j] = copy[i*width+j-1];
		}
	}

	delete[] copy;
}

void Channel::plus(int add){ // adiciona 'add' a todos os pixels
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(data[i*width+j] + add > 255){
				data[i*width+j] = 255;
			}
			else{
				data[i*width+j] += add;
			}
		}
	}
}

void Channel::plusNoLimit(int add){ // adiciona 'add' sem controle de ser maior que 255
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			data[i*width+j] += add;
		}
	}
}

void Channel::minus(int sub){ // subtrai 'sub' de todos os pixels
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(data[i*width+j] - sub < 0){
				data[i*width+j] = 0;
			}
			else{
				data[i*width+j] -= sub;
			}
		}
	}
}

void Channel::minusNoLimit(int add){ // adiciona 'add' sem controle de ser maior que 255
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			data[i*width+j] -= add;
		}
	}
}

void Channel::threshold(int threshold){ // transforma em imagem binária com o turnin point de 'threshold'
	for(int k = 0; k < height; k++){
		for(int l = 0; l < width; l++){
			if(data[k*width+l] > threshold){
				data[k*width+l] = 255;
			}
			else{
				data[k*width+l] = 0;
			}
		}
	}
}

void Channel::thresholdRange(int low, int high){ // zera os valores que estiverem fora do range low-high
	if(low > high){
		std::cout << "Invalid input";
		return;
	}
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(data[i*width+j] > high){
				data[i*width+j] = 0;
			}
			else if(data[i*width+j] < low){
				data[i*width+j] = 0;
			}
		}
	}
}

void Channel::maximizer(){ // expande o range atual para 0-255
	for (int i = 0; i < width * height; i++){
		data[i] = (unsigned char) ((((float)data[i]-min) / (float)(max-min)) * 255);
	}
}

// -----------------------------------FILTERS-----------------------------------

void Channel::applyFilterNoLimit(){ // aplica convolução com o filtro csv carregado no momento
	int sum;
	int truewidth = width + 2*offset;

	copy = expandedCopy();

	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filter[(k+offset)*dimension+(l+offset)];
				}
			}
			data[(i-offset)*width+(j-offset)] = sum/weight;
		}
	}

	delete[] copy;
}

void Channel::applyFilter(){ // trunca em 0 ou 255
	int sum;
	int truewidth = width + 2*offset;
	float out;

	copy = expandedCopy();

	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filter[(k+offset)*dimension+(l+offset)];
				}
			}
			out = sum/weight;
			if(out > 255) out = 255;
			else if(out < 0) out = 0;
			data[(i-offset)*width+(j-offset)] = (unsigned char)out;
		}
	}

	delete[] copy;
}

void Channel::applyFilterNoLimit(std::string csvfile){ // carrega um filtro e aplica
	loadCSV(csvfile);
	applyFilterNoLimit();
}

void Channel::applyFilter(std::string csvfile){ // carrega um filtro e aplica
	loadCSV(csvfile);
	applyFilter();
}

int* Channel::returnFiltered(){
	int sum;
	int truewidth = width + 2*offset;

	copy = expandedCopy();
	int *out = new int[width * height];
	
	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filter[(k+offset)*dimension+(l+offset)];
				}
			}
			out[(i-offset)*width+(j-offset)] = sum/weight;
		}
	}

	delete[] copy;

	return out;
}

void Channel::applyKernel(unsigned int index){ // faz convolução no pixel indicado
	unsigned int i = (index / width) + offset;
	unsigned int j = (index % width) + offset;

	int truewidth = width + 2*offset;

	int sum = 0;
	float out;

	for(int k = 0-offset; k < dimension-offset; k++){
		for(int l = 0-offset; l < dimension-offset; l++){
			sum += copy[(i+k)*truewidth+(j+l)]*
			filter[(k+offset)*dimension+(l+offset)];
		}
	}

	out = sum/weight;
	if(out > 255) out = 255;
	else if(out < 0) out = 0;
	data[(i - offset) * width + (j - offset)] = (unsigned char)out;
	
}

void Channel::applyNoLinear(float div){ // precisa de um csv carregado
	copy = expandedCopy(); // cria uma cópia expandida para os cálculos

	int sum;
	float mean, var, hotspot;
	int truewidth = width + 2*offset;

	//float min = 99999;
	//float max = -min;

	int filter4 = filter[4]; // só funciona com kernel de 3x3 assim
	int f_weight = weight;

	//Channel sd(width, height); // imagem para printar o desvio padrão pontual

	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)];
				}
			}
			mean = sum/(dimension*dimension); // média

			var = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					var += pow(copy[(i+k)*truewidth+(j+l)] - mean,2);
				}
			}

			var /= (dimension*dimension);
			var = sqrt(var); // desvio padrão

			//sd.set((i-offset)*width+(j-offset),var); // imagem desvio padrão

			//std::cout << var << std::endl;
			//if(var > 30) var = 30;
			//if(var > 20) var = 20;

			//hotspot = exp(var/(div*M_E)); // fórmula alternativa

			hotspot = exp(var/div); // e^(var/div)  // FÓRMULA PRINCIPAL

			if(hotspot > 500) hotspot = 500; // valores muito altos não fazem diferença

			weight = f_weight;
			weight += (hotspot - filter4); // ajusta o peso do kernel
			//std::cout << hotspot << std::endl;

			filter[4] = hotspot; // altera o valor do centro do kernel

			applyKernel((i-offset)*width+(j-offset)); // aplica a convolução

			//if(var < min) min = var;
			//if(var > max) max = var;
		}
	}

	//std::cout << "min: " << min << "\n" << "max: " << max << std::endl;
	
	//sd.saveFile("sd.pgm");

	delete[] copy;
	reloadCSV();
}

// -----------------------------------EDGE-DETECTION-----------------------------------

void Channel::sobel_old(){ // detecção de linhas pelo método de sobel
	int sum;
	float avrg1, avrg2;

	dimension = 3;
	offset = 1;
	weight = 4;
	int truewidth = width + 2*offset;
	int out;

	char filterX[9] = {
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1
	};
	char filterY[9] = {
		-1,-2,-1,
		 0, 0, 0,
		 1, 2, 1
	};

	copy = expandedCopy();

	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filterX[(k+offset)*dimension+(l+offset)];
				}
			}
			avrg1 = sum/weight;
			avrg1 *= avrg1;
			sum = 0;
			for(int k = 0-offset; k < dimension-offset; k++){
				for(int l = 0-offset; l < dimension-offset; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filterY[(k+offset)*dimension+(l+offset)];
				}
			}
			avrg2 = sum/weight;
			avrg2 *= avrg2;
			out = sqrt(avrg1+avrg2);
			if(out > 255) out = 255;
			data[(i-offset)*width+(j-offset)] = out;
		}
	}

	delete[] copy;
}

void Channel::sobel(){
	dimension = 3;
	offset = 1;
	weight = 1;

	int filterX[9] = {
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1
	};
	int filterY[9] = {
		-1,-2,-1,
		 0, 0, 0,
		 1, 2, 1
	};

	int *v[2];
	delete[] filter;
	filter = filterX;
	v[0] = returnFiltered();

	filter = filterY;
	v[1] = returnFiltered();

	/*
	data = to255(v[0]);
	saveFile("x.pgm");

	data = to255(v[1]);
	saveFile("y.pgm");
	*/

	int *result = new int[width * height];
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			result[i*width+j] = sqrt(
				pow(v[0][i*width+j],2) +
				pow(v[1][i*width+j],2)
				);
		}
	}

	delete[] data;
	data = to255(result);
	delete[] result;

	delete[] v[0];
	delete[] v[1];
	filter = nullptr;
}

void Channel::roberts(){
	dimension = 3;
	offset = 1;
	weight = 1;

	int filterX[4] = {
		 0, 1,
		-1, 0
	};
	int filterY[4] = {
		1,  0,
		0, -1
	};

	int *v[2];
	v[0] = new int[width * height];
	v[1] = new int[width * height];

	copy = expandedCopy();

	int sum;
	int truewidth = width + 2*offset;
	unsigned int index;
	int *result = new int[width * height];

	for(int i = offset; i < height+offset; i++){
		for(int j = offset; j < width+offset; j++){
			index = (i-offset)*width+(j-offset);
			sum = 0;
			for(int k = 0; k <= 1; k++){
				for(int l = 0; l <= 1; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filterX[k*2+l];
				}
			}
			v[0][index] = sum;
			
			sum = 0;
			for(int k = 0; k <= 1; k++){
				for(int l = 0; l <= 1; l++){
					sum += copy[(i+k)*truewidth+(j+l)]*filterY[k*2+l];
				}
			}
			v[1][index] = sum;
			
			result[index] = sqrt(pow(v[0][index], 2) + pow(v[1][index], 2));
		}
	}

	/*	
	data = to255(v[0]);
	saveFile("d1.pgm");

	data = to255(v[1]);
	saveFile("d2.pgm");
	*/

	delete[] data;
	data = to255(result);
	delete[] result;

	delete[] copy;
	delete[] v[0];
	delete[] v[1];
}

void Channel::robinson(){
	dimension = 3;
	offset = 1;
	weight = 1;

	int *filters[4];
	int filterX[9] = {
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1
	};
	int filterY[9] = {
		-1,-2,-1,
		 0, 0, 0,
		 1, 2, 1
	};
	int filterD1[9] = {
		-2,-1, 0,
		-1, 0, 1,
		 0, 1, 2
	};
	int filterD2[9] = {
		 0,-1,-2,
		 1, 0,-1,
		 2, 1, 0
	};

	filters[0] = filterX;
	filters[1] = filterY;
	filters[2] = filterD1;
	filters[3] = filterD2;

	int *v[8];
	delete[] filter;

	for(unsigned i = 0; i < 4; i++){
		filter = filters[i];
		v[i] = returnFiltered();
	}

	negativeFilter(filterX,9);
	negativeFilter(filterY,9);
	negativeFilter(filterD1,9);
	negativeFilter(filterD2,9);

	for(unsigned i = 0; i < 4; i++){
		filter = filters[i];
		v[i+4] = returnFiltered();
	}

	/*for(unsigned i = 0; i < 8; i++){
		data = to255(v[i]);
		saveFile(std::to_string(i)+".pgm");
	}*/

	int sum;
	int *result = new int[width * height];
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			sum = 0;
			for(int k = 0; k < 8; k++){
				sum += pow(v[k][i*width+j],2);
			}
			result[i*width+j] = sqrt(sum);
		}
	}

	delete[] data;
	data = to255(result);
	delete[] result;

	for(unsigned i = 0; i < 8; i++){
		delete[] v[i];
	}
	filter = nullptr;
}

// -----------------------------------OTHERS-----------------------------------

unsigned char* Channel::expandedCopy(){ // retorna um vetor com a imagem original mais as linhas de fora replicadas
	// é feito com base no filtro carregado no momento pois cada tamanho de filtro demanda 'offset' linhas em cada direção
	int truewidth = width + (2*offset);
	int trueheight = height + (2*offset);

	unsigned char *expCopy = new unsigned char[truewidth*trueheight] {0};

	int i, j;
	for(i = offset; i < height+offset; i++){
		for(j = offset; j < width+offset; j++)
		{
			expCopy[i*truewidth+j] = data[(i-offset)*width+(j-offset)];
		}
	}
	for(i = 0; i < offset; i++){					// primeiras linhas
		for(j = offset; j < truewidth-offset; j++){
			expCopy[i*truewidth+j] = data[j-offset];
		}
	}
	for(i = height+offset; i < trueheight; i++){	// últimas linhas
		for(j = offset; j < truewidth-offset; j++){
			expCopy[i*truewidth+j] = data[(height-1)*width+j-offset];
		}
	}
	for(i = offset; i < height+offset; i++){		// primeiras colunas
		for(j = 0; j < offset; j++){
			expCopy[i*truewidth+j] = data[(i-offset)*width];
		}
	}
	for(i = offset; i < height+offset; i++){		// últimas colunas
		for(j = width+offset; j < truewidth; j++){
			expCopy[i*truewidth+j] = data[(i-offset)*width+(width-1)];
		}
	}
	
	for(i = 0; i < 2; i++){							// 4 cantos
		for(j = 0; j < 2; j++){
			for(int k = 0; k < offset; k++){
				for(int l = 0; l < offset; l++){
					expCopy[(((i*(height+offset))+k)*truewidth) + (j*(width+offset) + l)] = 
					data[0 + ((i*(height-1))*width) + (j*(width-1))];
				}
			}
		}
	}
	return expCopy;
}

void Channel::update(){ // atualizada parâmetros min e max
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			if(data[i*width+j] > max){
				max = data[i*width+j];
			}
			if(data[i*width+j] < min){
				min = data[i*width+j];
			}
		}
	}
}

unsigned char* Channel::to255(int *in){
	int min = 99999;
	int max = -min;
	for(int index = 0; index < width*height; index++){
		if(in[index] > max) max = in[index]; 
        if(in[index] < min) min = in[index]; 
	}

	//std::cout << "max e min: " << max << " " << min << std::endl;
	unsigned char *out = new unsigned char[width * height];
	for(int index = 0; index < width*height; index++){
		out[index] = (unsigned char) ((((float)in[index]-min) / (float)(max-min)) * 255);
	}

	return out;
}

unsigned char* Channel::int2char(int *in){
	unsigned char *out = new unsigned char[width * height];
	for(int index = 0; index < width*height; index++){
		out[index] = (unsigned char) in[index];
	}

	return out;
}

void negativeFilter(int *in, int length){
	for(int i = 0; i < length; i++){
		in[i] = -in[i];
	}
}

// -----------------------------------DEBUG-----------------------------------

void Channel::showData(int a){ // printa diversos parâmetros para debug
	// a = 1 base
	// a = 2 printa o vetor da imagem
	// a = 3 printa o vetor do filtro
	// a = 4 printa ambos
	// sempre são printados parâmetros 'pequenos'
	std::cout << std::endl;
	if(a % 2 == 0){
		for (int i = 0; i < height; i++){
			std::cout << "Line " << i << ": ";
			for (int j = 0; j < width; j++){
				std::cout << +data[i*width+j] << " ";
			}
			std::cout << std::endl;
		}
	}

	if(a >= 3){
		for (int i = 0; i < dimension; i++){
			std::cout << "Line " << i << ": ";
			for (int j = 0; j < dimension; j++){
				std::cout << filter[i*dimension+j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "Weight: " << weight << std::endl;
	}

	std::cout << "Path:\t" 	<< path 	<< std::endl;
	std::cout << "Width:\t" << width 	<< std::endl;
	std::cout << "Height:\t"<< height 	<< std::endl;
	std::cout << "Range:\t" << range 	<< std::endl;
	std::cout << "Type:\t" 	<< type 	<< std::endl;

	std::cout << "CSVPath:" << csv_path	<< std::endl;
	std::cout << "Offset:\t"<< offset 	<< std::endl;
	std::cout << "Weight:\t"<< weight 	<< std::endl;

	std::cout << "Min:\t" 	<< min 		<< std::endl;
	std::cout << "Max:\t" 	<< max 		<< std::endl;
	//std::cout << "" <<  << std::endl;
}