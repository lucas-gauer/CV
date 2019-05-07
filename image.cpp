#include "image.h"

double covariance(Channel *X, Channel *Y, double mX, double mY);
void inverse(double* in, double* out);
void matrixMultiplier(double* a, double* b, double* result);
void sort(float *v, int start, int end);

// -----------------------------------CONSTRUCTORS-----------------------------------

Image::Image()
:width(0),height(0),range(255),type(0)
{}

Image::Image(std::string o_path)
:width(0),height(0),range(255),type(0)
{ // load genérico
	loadFile(o_path);
}

Image::~Image(){
	delete R;
	delete G;
	delete B;
}

// -----------------------------------IO-----------------------------------

void Image::loadFile(std::string o_path){ // cria 3 sub-imagens Channel para R, G e B
	unsigned char *r, *g, *b;
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
		std::cout << "This class is for color pictures only :(\n";
		exit(1);
	}
	else if(strbuff == "P3") type = 3;

	inFile >> width;
	inFile >> height;
	inFile >> range;

	r = new unsigned char[width * height];
	g = new unsigned char[width * height];
	b = new unsigned char[width * height];

	int i = 0;
	while(i < width * height){ // linhas de conteúdo
		inFile >> strbuff;
		r[i] = stoi(strbuff);
		inFile >> strbuff;
		g[i] = stoi(strbuff);
		inFile >> strbuff;
		b[i] = stoi(strbuff);
		++i;
	}

	delete R;
	delete G;
	delete B;
	R = new Channel(r,width,height);
	G = new Channel(g,width,height);
	B = new Channel(b,width,height);

	inFile.close();

	// name part:
	if(o_path[0] == '.') i = 2;
	else i = 0;
	while(o_path[i] != '.'){
		++i;
	}
	int j = i;
	while(o_path[i] != '/' && i > 0){
		--i;
	}
	if(o_path[0] == '.') ++i;

	name = o_path.substr(i, j - i);
}

void Image::reloadFile(){
	loadFile(path);
	reloadCSV();
}

void Image::saveFile(){ // salva o estado atual
	saveFile(name + ".ppm");
}

void Image::saveFile(std::string o_path){ // salva em um caminho específico
	std::ofstream outFile;

	outFile.open(o_path);
	outFile << "P" << type << "\n";
	outFile << width << " " << height << " " << range << "\n";
	
	for(int k = 0; k < height; ++k){
		for(int l = 0; l < width; ++l){
			outFile << +R->data[k*width + l] << "\n";
			outFile << +G->data[k*width + l] << "\n";
			outFile << +B->data[k*width + l] << "\n";
		}
	}

	//R->saveFile("R.pgm");
	//G->saveFile("G.pgm");
	//B->saveFile("B.pgm");

	outFile.close();
}

void Image::loadCSV(std::string o_path){
	csv_path = o_path;
	R->loadCSV(o_path);
	G->loadCSV(o_path);
	B->loadCSV(o_path);
}

void Image::reloadCSV(){
	if(!csv_path.empty()) loadCSV(csv_path);
}

// -----------------------------------MODIFICATIONS-----------------------------------

void Image::negative(){
	R->negative();
	G->negative();
	B->negative();
}

void Image::mirror(){
	R->mirror();
	G->mirror();
	B->mirror();
}

void Image::plus(int add){
	R->plus(add);
	G->plus(add);
	B->plus(add);
}

void Image::plusNoLimit(int add){
	R->plusNoLimit(add);
	G->plusNoLimit(add);
	B->plusNoLimit(add);
}

void Image::minus(int sub){
	R->minus(sub);
	G->minus(sub);
	B->minus(sub);
}

void Image::minusNoLimit(int sub){
	R->minusNoLimit(sub);
	G->minusNoLimit(sub);
	B->minusNoLimit(sub);
}

void Image::threshold(int threshold){
	R->threshold(threshold);
	G->threshold(threshold);
	B->threshold(threshold);
}

void Image::maximizer(){ // TODO

}

// -----------------------------------FILTERS-----------------------------------

void Image::applyFilterNoLimit(){
	R->applyFilterNoLimit();
	G->applyFilterNoLimit();
	B->applyFilterNoLimit();
}

void Image::applyFilter(){
	R->applyFilter();
	G->applyFilter();
	B->applyFilter();
}

void Image::applyFilterNoLimit(std::string csvfile){
	loadCSV(csvfile);
	applyFilterNoLimit();
}

void Image::applyFilter(std::string csvfile){
	loadCSV(csvfile);
	applyFilter();
}

/*void Image::applyKernel(unsigned int index){
	R->applyKernel(index);
	G->applyKernel(index);
	B->applyKernel(index);
}*/

void Image::applyNoLinear(float div){
	R->applyNoLinear(div);
	G->applyNoLinear(div);
	B->applyNoLinear(div);
}

// -----------------------------------DISTANCES-----------------------------------

float Image::distance(unsigned int a, unsigned int b){ //linearDistance
	return sqrt(
		pow(R->data[a] - R->data[b],2) + 
		pow(G->data[a] - G->data[b],2) + 
		pow(B->data[a] - B->data[b],2)
		);
}

void Image::L1(int refR, int refG, int refB){ // distância L1
	Channel result(width, height);

	int aux;
	for(int i = 0; i < height; ++i){
		for(int j = 0; j < width; ++j){
			aux = 	abs(R->data[i*width + j] - refR) +
					abs(G->data[i*width + j] - refG) +
					abs(B->data[i*width + j] - refB);
			if(aux > 255) aux = 255;
			result.data[i*width + j] = aux;
		}
	}

	result.saveFile("L1.pgm");
	result.threshold(50);
	result.saveFile("Bin-L1.pgm");
}

void Image::L2(int refR, int refG, int refB){ // distância L2
	Channel result(width, height);

	int aux;
	for(int i = 0; i < height; ++i){
		for(int j = 0; j < width; ++j){
			aux = sqrt(
				pow(R->data[i*width + j] - refR, 2) +
				pow(G->data[i*width + j] - refG, 2) +
				pow(B->data[i*width + j] - refB, 2)
				);
			if(aux > 255) aux = 255;
			result.data[i*width + j] = aux;
		}
	}

	result.saveFile("L2.pgm");
	result.threshold(50);
	result.saveFile("Bin-L2.pgm");
}

void Image::Mahalanobis(std::string o_path){ // distância de Mahalanobis
	Image ref(o_path); // 'imagem' com amostras de referências
	double n_samples = ref.width * ref.height;

	unsigned int add[3] = {0};
	for(int i = 0; i < n_samples; ++i){ // media das ref
		add[0] += ref.R->data[i]; // get(Cor,index)
		add[1] += ref.G->data[i]; // 0 = R, 1 = G, 2 = B
		add[2] += ref.B->data[i];
	}

	double mean[3];
	mean[0] = add[0]/n_samples;
	mean[1] = add[1]/n_samples;
	mean[2] = add[2]/n_samples;

	double cov[9]; // vetor 'matriz' de covariâncias
	cov[0] = covariance(ref.R, ref.R, mean[0], mean[0]);// XX
	cov[1] = covariance(ref.R, ref.G, mean[0], mean[1]);// XY
	cov[2] = covariance(ref.R, ref.B, mean[0], mean[2]);// XZ
	cov[3] = cov[1];
	cov[6] = cov[2];
	cov[4] = covariance(ref.G, ref.G, mean[1], mean[1]);// YY
	cov[5] = covariance(ref.G, ref.B, mean[1], mean[2]);// YZ
	cov[7] = cov[5];
	cov[8] = covariance(ref.B, ref.B, mean[2], mean[2]);// ZZ

	double invcov[9]; // cov^-1
	inverse(cov, invcov); // (in, out)

	Channel result(width, height);
	Channel bin(width, height);

	double diff[3];
	double q[3]; // resultado da multiplicação
	double dist; // distancia
	for(int a = 0; a < width * height; ++a){
		diff[0] = (double)R->data[a] - mean[0];
		diff[1] = (double)G->data[a] - mean[1];
		diff[2] = (double)B->data[a] - mean[2];

		matrixMultiplier(invcov, diff, q);

		dist = 0;
		for(int i = 0; i < 3; i++){
			dist += (diff[i] * q[i]);
		}

		dist = sqrt(dist);

		if(dist > 255) dist = 255;
		result.data[a] = round(dist);

		if(dist > 3) bin.data[a] = 255;
		else bin.data[a] = 0;
	}

	result.saveFile("Mahalanobis.pgm");
	result.maximizer();
	result.saveFile("Max-Mahalanobis.pgm");
	
	bin.saveFile("Bin-Mahalanobis.pgm");
}

void Image::Knn(std::string o_path, int K){
	Image ref(o_path); // 'imagem' com amostras de referências
	int n_samples = ref.width * ref.height;

	float aux, min_dist, mean_dist = 0;
	float *NND = new float[n_samples];
	for(int i = 0; i < n_samples; ++i){
		min_dist = 9999;
		for(int j = 0; j < n_samples; ++j){
			aux = ref.distance(i, j);
			if(aux < min_dist && aux != 0) min_dist = aux;
		}
		mean_dist += NND[i] = min_dist;
	}

	delete[] NND;
	mean_dist = mean_dist/(float)n_samples; // distância média na amostra

	Channel result(width, height);
	Channel bin(width, height);

	float knn_mean_dist;
	float *dist = new float[n_samples];
	for(int i = 0; i < width * height; ++i){
		min_dist = 9999;
		for(int s = 0; s < n_samples; ++s){
			dist[s] = sqrt(
				pow(R->data[i] - ref.R->data[s],2) +
				pow(G->data[i] - ref.G->data[s],2) +
				pow(B->data[i] - ref.B->data[s],2)
				);
			if(dist[s] < min_dist){
				min_dist = dist[s];
			}
		}
		sort(dist, 0, n_samples);

		knn_mean_dist = 0;
		for(int k = 0; k < K; ++k){
			knn_mean_dist += dist[k];
		}
		knn_mean_dist /= K;

		if(min_dist > 255) min_dist = 255;
		result.data[i] = min_dist;

		if((knn_mean_dist/mean_dist) > 3) bin.data[i] = 255;
		else bin.data[i] = 0;
	}
	
	delete[] dist;

	result.saveFile("Knn.pgm");
	result.maximizer();
	result.saveFile("Max-Knn.pgm");

	bin.saveFile("Bin-Knn.pgm");
}

// -----------------------------------EDGE-DETECTION-----------------------------------

void Image::sobel(){
	R->sobel();
	G->sobel();
	B->sobel();
}

void Image::roberts(){
	R->roberts();
	G->roberts();
	B->roberts();
}

void Image::robinson(){
	R->robinson();
	G->robinson();
	B->robinson();
}

void Image::canny(int maxDist){
	R->canny(maxDist);
	G->canny(maxDist);
	B->canny(maxDist);
}

// -----------------------------------FILL-----------------------------------

void Image::fill(int threshold, bool adapt){
	unsigned int size = width * height;

	unsigned char *father = new unsigned char[size]{0};
	int *groups = new int[size]{0};

	unsigned int index  = 0;
	unsigned int id 	= 0;

	while(index < size){
		++id;

		flood(index, groups, father, id, threshold, adapt);

		while(index < size && groups[index] != 0){
			++index;
		}
	}
	std::cout << "Groups: " << id << std::endl;

	unsigned int *r 	= new unsigned int[id]{0};
	unsigned int *g 	= new unsigned int[id]{0};
	unsigned int *b 	= new unsigned int[id]{0};
	unsigned int *count = new unsigned int[id]{0};

	for(unsigned i = 0; i < size; ++i){
		r[groups[i]-1] += R->data[i];
		g[groups[i]-1] += G->data[i];
		b[groups[i]-1] += B->data[i];
		++count[groups[i]-1];
	}

	for(unsigned i = 0; i < id; ++i){
		if(count[i] == 0) continue;
		r[i] /= count[i];
		g[i] /= count[i];
		b[i] /= count[i];
	}

	for(unsigned i = 0; i < size; ++i){
		R->data[i] = r[groups[i]-1];
		G->data[i] = g[groups[i]-1];
		B->data[i] = b[groups[i]-1];
	}

	segmentEdges(groups, id+1, false);

	delete[] groups;
	delete[] father;
	delete[] r;
	delete[] g;
	delete[] b;
	delete[] count;
}

void Image::flood(int index, int *groups, unsigned char *father, int id, int th, bool adapt){

	groups[index] = id;
	int reference = index;

	bool moving = true;
	bool flag   = true;
	while(flag){
		moving = false;

		while(Try(reference, index, groups, th, LEFT)){
			index -= 1;
		    moving = true;
		    groups[index] = id;
		    father[index] = RIGHT;
		    if(adapt) reference = index;
		}

		while(Try(reference, index, groups, th, UP)){
		    index -= width;
		    moving = true;
		    groups[index] = id;
		    father[index] = DOWN;
		    if(adapt) reference = index;
		}

		while(Try(reference, index, groups, th, RIGHT)){
		    index += 1;
		    moving = true;
		    groups[index] = id;
		    father[index] = LEFT;
		    if(adapt) reference = index;
		}

		while(Try(reference, index, groups, th, DOWN)){
		    index += width;
		    moving = true;
		    groups[index] = id;
		    father[index] = UP;
		    if(adapt) reference = index;
		}

		if(!moving){
			switch(father[index]){
			case LEFT:
				index -= 1;
				break;
			case UP:
				index -= width;
				break;
			case RIGHT:
				index += 1;
				break;
			case DOWN:
				index += width;
				break;
			case 0:
				flag = false;
			}
			if(adapt) reference = index;
		}
	}
}

/*bool Image::Try(int reference, int index, int *groups, int &th, DIR dir){
	switch(dir){
	case LEFT:
		if(index % width != 0){ // se existe
			if(groups[index - 1] == 0){ // se não foi agrupado ainda
				if(distance(reference, index - 1) < th){ // se pertence
					return true;
				}
			}
		}
		break;

	case UP:
		if(index - width >= 0){
			if(groups[index - width] == 0){
				if(distance(reference, index - width) < th){
					return true;
				}
			}
		}
		break;

	case RIGHT:
		if(index % width != width - 1){
			if(groups[index + 1] == 0){
				if(distance(reference, index + 1) < th){
					return true;
				}
			}
		}
		break;

	case DOWN:
		if((index + width) < (height * width)){
			if(groups[index + width] == 0){
				if(distance(reference, index + width) < th){
					return true;
				}
			}
		}
		break;
	}

	return false;
}*/

bool Image::Try(int reference, int index, int *groups, int &th, DIR dir){
	switch(dir){
	case LEFT:
		if(index % width != 0 && groups[index - 1] == 0 
		&& distance(reference, index-1) < th){
			return true;
		}
		break;

	case UP:
		if(index - width >= 0 && groups[index - width] == 0 
		&& distance(reference, index - width) < th){
			return true;
		}
		break;

	case RIGHT:
		if(index % width != width - 1 && groups[index + 1] == 0
		&& distance(reference, index + 1) < th){
			return true;
		}
		break;

	case DOWN:
		if(index + width < height * width && groups[index + width] == 0
		&& distance(reference, index + width) < th){
			return true;
		}
		break;
	}

	return false;
}

void Image::segmentEdges(int *groups, unsigned int n_groups, bool colors){
	Channel result(width, height);

	unsigned char *color = new unsigned char[n_groups]{};
	if(colors){
		for(unsigned i = 0; i < n_groups; ++i){
			color[i] = (rand() % 25) * 10;
		}
	}
	else{
		for(unsigned i = 0; i < n_groups; ++i){
			color[i] = 255;
		}
	}

	for(int index = 0; index < width * height; ++index){
		if(index % width != width - 1){ 			// RIGHT
			if(groups[index] != groups[index + 1]){
				result.data[index] = color[groups[index]];
				result.data[index + 1] = color[groups[index + 1]];
			}
		}

		if(index + width < height * width){			// DOWN
			if(groups[index] != groups[index + width]){
				result.data[index] = color[groups[index]];
				result.data[index + width] = color[groups[index + width]];
			}
		}
	}

	delete[] color;

	if(colors) result.maximizer();
	result.saveFile(name + "-" + std::to_string(n_groups-1) + "-edges.pgm");
}

void Image::floodFrom(int x, int y, int threshold, bool adapt){
	unsigned int size  = width * height;

	unsigned char *father = new unsigned char[size]{0};
	int *groups = new int [size]{0};

	unsigned int index = y*width + x;
	unsigned int id = 1;

	flood(index, groups, father, id, threshold, adapt);

	unsigned int r, g, b, count;
	r = g = b = count = 0;

	for(unsigned i = 0; i < size; ++i){
		if(groups[i] != 1) continue;
		r += R->data[i];
		g += G->data[i];
		b += B->data[i];
		count++;
	}

	r /= count;
	g /= count;
	b /= count;
	
	for(unsigned i = 0; i < size; ++i){
		if(groups[i] != 1) continue;
		R->data[i] = r;
		G->data[i] = g;
		B->data[i] = b;
	}

	segmentEdges(groups, 2, true);

	delete[] groups;
	delete[] father;
}

// -----------------------------------OTHERS-----------------------------------

Channel* Image::grayscale(){ // converte para escala de cinza
	Channel *gray = new Channel(width, height);

	int aux;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			aux = 0;
			aux += R->data[i*width + j];
			aux += G->data[i*width + j];
			aux += B->data[i*width + j];
			aux = aux/3;
			gray->data[i*width + j] = aux;
		}
	}

	return gray;
}

// Mahanalobis subsection
double covariance(Channel *X, Channel *Y, double mX, double mY){
	double n_samples = X->width * X->height;

	unsigned int aux = 0;
	int i;
	for(i = 0; i < n_samples; i++){
		aux += (X->data[i] - mX) * (Y->data[i] - mY);
	}
	
	return (double)aux/n_samples;
}

void inverse(double* in, double* out){
	double det = 0;
	for(int i = 0; i < 3; i++){
		det += (in[i] * in[3+((i+1)%3)] * in[6+((i+2)%3)]) - 
			   (in[i] * in[3+((i+2)%3)] * in[6+((i+1)%3)]);
	}
	
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++)
			out[i*3+j] = ((in[3*((j+1)%3)+((i+1)%3)] * in[3*((j+2)%3)+((i+2)%3)]) - 
						  (in[3*((j+1)%3)+((i+2)%3)] * in[3*((j+2)%3)+((i+1)%3)]))/det;
	}
}

void matrixMultiplier(double* a, double* b, double* result){
	double add;
	for(int i = 0; i < 3; i++){
		add = 0;
		for(int j = 0; j < 3; j++){
			add += a[i*3 + j] * b[j];
		}
		result[i] = add;
	}
}

// Knn subsection
void intercalate(float* v, int start, int mid, int end){
	int i = start, j = mid , k = 0;
    float* aux = new float[end - start];
    while(i < mid && j < end){
        if(v[i] <= v[j]){
            aux[k++] = v[i++];
        }
        else{
            aux[k++] = v[j++];
        }
    }
    while(i < mid){
        aux[k++] = v[i++];
    }
    while(j < end){
        aux[k++] = v[j++];
    }
    for(i = start; i < end; i++){
        v[i] = aux[i-start];
    }
    delete[] aux;
}

void sort(float* v, int start, int end){ // mergesort
    int mid;
    if(start < end-1){
        mid = (start + end)/2;
        sort(v, start, mid);
        sort(v, mid, end);
        intercalate(v, start, mid, end);
    }
}

// -----------------------------------DEBUG-----------------------------------

void Image::showData(int a){
	
	std::cout << std::endl;
	if(a % 2 == 0){
		for (int i = 320; i < 322; i++){
			std::cout << "Line " << i << ": ";
			for (int j = 0; j < width; j++){
				std::cout << +R->data[i*width + j] << "-";
				std::cout << +G->data[i*width + j] << "-";
				std::cout << +B->data[i*width + j] << " ";
			}
			std::cout << std::endl;
		}
	}

	std::cout << "Path:\t" 	<< path 	<< std::endl;
	std::cout << "Width:\t" << width 	<< std::endl;
	std::cout << "Height:\t"<< height 	<< std::endl;
	std::cout << "Range:\t" << range 	<< std::endl;
	std::cout << "Type:\t" 	<< type 	<< std::endl;

	std::cout << "CSVPath:" << csv_path	<< std::endl;
	//std::cout << "" <<  << std::endl;
}
