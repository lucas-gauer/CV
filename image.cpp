#include "image.h"

double covariance(Channel *X, Channel *Y, double mX, double mY);
void inverse(double* in, double* out);
void matrixMultiplier(double* a, double* b, double* result);
void sort(float *v, int start, int end);

// -----------------------------------CONSTRUCTORS-----------------------------------

Image::Image()
:width(0),height(0),range(255),type(0)
{}

Image::Image(std::string path)
:width(0),height(0),range(255),type(0)
{ // load genérico
	loadFile(path);
}

Image::~Image(){

}

// -----------------------------------IO-----------------------------------

void Image::loadFile(std::string o_path){ // cria 3 sub-imagens da classe Channel para R, G e B
	unsigned char *r, *g, *b;
	std::ifstream inFile;
	std::string hue;
	int i = 0;

	path = std::string(o_path);

	inFile.open(path);
	if(!inFile){
		std::cout << "Unable to open file\n";
		exit(1);
	}
	else{
		//std::cout << "Loaded\n";
	}

	inFile >> hue;
	if(hue == "P2"){
		type = 2;
		std::cout << "This class is for color pictures only :(\n";
		exit(1);
	}
	else if(hue == "P3"){
		type = 3;
	}

	inFile >> width;
	inFile >> height;
	inFile >> range;

	//std::cout << width << std::endl;
	//std::cout << height << std::endl;
	//std::cout << range << std::endl;

	r = new unsigned char[width * height];
	g = new unsigned char[width * height];
	b = new unsigned char[width * height];

	i = 0;
	while(inFile){ // linhas de conteúdo
		inFile >> hue;
		r[i] = stoi(hue);
		inFile >> hue;
		g[i] = stoi(hue);
		inFile >> hue;
		b[i] = stoi(hue);
		i++;
	}

	R = new Channel(r,width,height);
	G = new Channel(g,width,height);
	B = new Channel(b,width,height);

	inFile.close();

}

void Image::reloadFile(){
	delete R;
	delete G;
	delete B;
	loadFile(path);
	reloadCSV();
}

void Image::saveFile(){ // salva o estado atual
	std::ofstream outFile;

	outFile.open("Output.ppm");
	outFile << "P" << type << "\n";
	outFile << width << " " << height << " 255\n";
	
	for(int k = 0; k < height; k++){
		for(int l = 0; l < width; l++){
			outFile << +R->data[k*width+l] << "\n";
			outFile << +G->data[k*width+l] << "\n";
			outFile << +B->data[k*width+l] << "\n";
		}
	}
	//R->saveFile("R.pgm");
	//G->saveFile("G.pgm");
	//B->saveFile("B.pgm");

	outFile.close();
}

void Image::saveFile(std::string o_path){ // salva em um caminho específico
	std::ofstream outFile;

	outFile.open(o_path);
	outFile << "P" << type << "\n";
	outFile << width << " " << height << " 255\n";
	
	for(int k = 0; k < height; k++){
		for(int l = 0; l < width; l++){
			outFile << +R->data[k*width+l] << "\n";
			outFile << +G->data[k*width+l] << "\n";
			outFile << +B->data[k*width+l] << "\n";
		}
	}
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

void Image::L1(int refR, int refG, int refB){ // distância L1
	Channel result(width, height);
	int aux;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){ // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			aux = 	abs(R->data[i*width+j] - refR) +
					abs(G->data[i*width+j] - refG) +
					abs(B->data[i*width+j] - refB);
			if(aux > 255) aux = 255;
			result.data[i*width+j] = aux;
		}
	}
	result.saveFile("L1.pgm");
	result.threshold(50);
	result.saveFile("Bin-L1.pgm");
}

void Image::L2(int refR, int refG, int refB){ // distância L2
	Channel result(width, height);
	int aux;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			aux = sqrt((R->data[i*width+j] - refR) * (R->data[i*width+j] - refR) +
					   (G->data[i*width+j] - refG) * (G->data[i*width+j] - refG) +
					   (B->data[i*width+j] - refB) * (B->data[i*width+j] - refB));
			if(aux > 255) aux = 255;
			result.data[i*width+j] = aux;
		}
	}
	result.saveFile("L2.pgm");
	result.threshold(50);
	result.saveFile("Bin-L2.pgm");
}

void Image::Mahalanobis(std::string path){ // distância de Mahalanobis
	//applyFilter("../csv/gauss.csv");
	Image ref(path); // 'imagem' com amostras de referências
	double n_samples = ref.width * ref.height;
	//std::cout << "n_samples: " << n_samples << std::endl;
	double mean[3];
	unsigned int add[3] = {0};
	int i;
	for(i = 0; i < n_samples; i++){ // media das ref
		add[0] += ref.R->data[i]; // get(Cor,index)
		add[1] += ref.G->data[i];	// 0 = R, 1 = G, 2 = B
		add[2] += ref.B->data[i];
	}

	mean[0] = add[0]/n_samples;
	mean[1] = add[1]/n_samples;
	mean[2] = add[2]/n_samples;

	double cov[9]; // vetor 'matriz' de covariâncias
	double invcov[9]; // cov^-1
	cov[0] = covariance(ref.R, ref.R, mean[0], mean[0]);// XX
	cov[1] = covariance(ref.R, ref.G, mean[0], mean[1]);// XY
	cov[2] = covariance(ref.R, ref.B, mean[0], mean[2]);// XZ
	cov[3] = cov[1];
	cov[6] = cov[2];
	cov[4] = covariance(ref.G, ref.G, mean[1], mean[1]);// YY
	cov[5] = covariance(ref.G, ref.B, mean[1], mean[2]);// YZ
	cov[7] = cov[5];
	cov[8] = covariance(ref.B, ref.B, mean[2], mean[2]);// ZZ

	inverse(cov, invcov);

	Channel result(width, height);
	Channel bin(width, height);
	double diff[3];
	double* q = new double[3];
	double dist;

	for(int a = 0; a < width*height; a++){
		dist = 0;

		diff[0] = (double)R->data[a] - mean[0];
		diff[1] = (double)G->data[a] - mean[1];
		diff[2] = (double)B->data[a] - mean[2];

		matrixMultiplier(invcov, diff, q);

		for(int i = 0; i < 3; i++){
			dist += (diff[i]*q[i]);
		}

		dist = sqrt(dist);

		if(dist > 255) dist = 255;

		result.data[a] = round(dist);

		if(dist > 3) bin.data[a] = 255;
		else bin.data[a] = 0;
	}

	result.saveFile("Mahalanobis.pgm");
	result.update();
	result.maximizer();
	result.saveFile("Max-Mahalanobis.pgm");
	//result.threshold(3);
	bin.saveFile("Bin-Mahalanobis.pgm");
}

void Image::Knn(std::string path, int K){
	Image ref(path); // 'imagem' com amostras de referências
	int n_samples = ref.width * ref.height;
	//std::cout << "n_samples: " << n_samples << std::endl;

	float aux, min_dist, mean_dist = 0;
	float NND[n_samples];
	for(int i = 0; i < n_samples; i++){
		min_dist = 9999;
		for(int j = 0; j < n_samples; j++){
			aux = sqrt(
				pow(ref.R->data[i] - ref.R->data[j],2) +
				pow(ref.G->data[i] - ref.G->data[j],2) +
				pow(ref.B->data[i] - ref.B->data[j],2)
				);
			//if(aux == 0) continue;
			if(aux < min_dist && aux != 0) min_dist = aux;
		}
		//std::cout << min_dist << std::endl;
		mean_dist += NND[i] = min_dist;
	}

	mean_dist = mean_dist/(float)n_samples;
	//std::cout << "mean_dist: " << mean_dist << std::endl;

	Channel result(width, height);
	Channel bin(width, height);
	float knn_mean_dist;
	float dist[n_samples];
	for(int i = 0; i < width*height; i++){
		min_dist = 9999;
		for(int s = 0; s < n_samples; s++){
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
		for(int k = 0; k < K; k++){
			knn_mean_dist += dist[k];
			//std::cout << "dist[k]: " << dist[k] << std::endl;
		}
		knn_mean_dist /= K;
		//std::cout << "knn_mean_dist:  " << knn_mean_dist << std::endl;
		//std::cout << "min_dist:  " << min_dist << std::endl;
		//std::cout << "NND[near]: " << NND[near] << std::endl;

		if(min_dist > 255) min_dist = 255;
		result.data[i] = min_dist;
		if((knn_mean_dist/mean_dist) > 3) bin.data[i] = 255;
		else bin.data[i] = 0;
	}
	result.saveFile("Knn.pgm");
	result.update();
	result.maximizer();
	result.saveFile("Max-Knn.pgm");
	//result.threshold(mean_dist*3);
	//result.saveFile(folder + "/Bin-L3-2.pgm");
	bin.saveFile("Bin-Knn.pgm");
}

// -----------------------------------EDGE-DETECTION-----------------------------------

void Image::sobel(){
	R->sobel();
	G->sobel();
	B->sobel();
	//R->saveFile("1.pgm");
	//G->saveFile("2.pgm");
	//B->saveFile("3.pgm");
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

// -----------------------------------OTHERS-----------------------------------

Channel* Image::grayscale(){ // converte para escala de cinza
	Channel *gray = new Channel(width, height);
	int aux;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			aux = 0;
			aux += R->data[i*width+j];
			aux += G->data[i*width+j];
			aux += B->data[i*width+j];
			aux = aux/3;
			gray->data[i*width+j] = aux;
		}
	}
	return gray;
}

// Mahanalobis subsection
double covariance(Channel *X, Channel *Y, double mX, double mY){
	double n_samples = X->width * X->height;
	//std::cout << "n_samples: " << n_samples << std::endl;
	unsigned int aux = 0;
	int i;
	for(i = 0; i < n_samples; i++){
		aux += (X->data[i] - mX) * (Y->data[i] - mY);
	}
	
	//std::cout << "aux: " << aux << std::endl;
	return (double)aux/n_samples;
}

void inverse(double* in, double* out){
	double det = 0;
	for(int i = 0; i < 3; i++){
		det += (in[i] * in[3+((i+1)%3)] * in[6+((i+2)%3)]) - 
			   (in[i] * in[3+((i+2)%3)] * in[6+((i+1)%3)]);
	}
	//std::cout<<"\n\ndet: "<<det;
	//std::cout<<"\n\n1/det: "<<1/det;
	
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
			add += a[i*3+j]*b[j];
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
	/*
	std::cout << std::endl;
	if(a % 2 == 0){
		for (int i = 0; i < height; i++){
			std::cout << "Line " << i << ": ";
			for (int j = 0; j < width; j++){
				std::cout << +data[i*width+j] << " ";
			}
			std::cout << std::endl;
		}
	}*/

	std::cout << "Path:\t" 	<< path 	<< std::endl;
	std::cout << "Width:\t" << width 	<< std::endl;
	std::cout << "Height:\t"<< height 	<< std::endl;
	std::cout << "Range:\t" << range 	<< std::endl;
	std::cout << "Type:\t" 	<< type 	<< std::endl;

	std::cout << "CSVPath:" << csv_path	<< std::endl;
	//std::cout << "" <<  << std::endl;
}