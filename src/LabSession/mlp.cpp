/*****************************************************************************/
/* Multi-Layer Perceptron                                                    */
/*****************************************************************************/
/* Written by Benoit Huet in March 2001                                      */
/* Institut Eurecom                                                          */
/* benoit.huet@eurecom.fr                               www.eurecom.fr/~huet */
/*****************************************************************************/

//#include "stdafx.h"
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "Image.h"

using namespace std;

//Variables you may wish to change the value of.
//
double learning_rate =			0.5; // 0.5 for the sigmoid and 0.01 for tanh are usually good choices
int MAXCYCLES =					100; //defines the maximum number of training cycles.
double RMS_ERROR_THRESHOLD =	0.0001;  //defines the stopping threshold during training.

//
//Variables you should not modify
//
int nb_layers;
struct Layer{
	int size; // number of neuron in that unit (includes the neuron employed as bias except for the output layer which does not require an extra neuron)
	double ** w; // the connection weights arriving at this layer (note that the input layer does not have weights)
	double * o; // contains the output value of this layer's neurons (after activation function)
	double * u; // contain the net input (activation) value of a neuron 
	double * dE_du; //contains the partial derivative of the error with respect to the activation of the corresponding neuron (only used during backpropagation)
};
struct Layer * layers;
double ** pattern=NULL;
int nb_pattern=0;
double rms_error; //the root mean square error of the networks

//
// Some Memory Allocation Functions
//
#define Alloc2D(Variable,Type,n1,n2) { Variable=(Type **)malloc(sizeof(Type *)*n1); Variable[0]=(Type *)malloc(sizeof(Type)*n1*n2); for(unsigned int A2D=1;A2D<n1;A2D++) Variable[A2D]=Variable[A2D-1]+n2; }
#define Free2D(Variable) { delete Variable[0]; delete Variable;}
#define Resize2D(Variable,Type,n1,n2)	{ Variable=(Type **)realloc(Variable,sizeof(Type *)*n1); Variable[0]=(Type *)realloc(Variable[0],sizeof(Type)*n1*n2); for(int A2D=1;A2D<n1;A2D++) Variable[A2D]=Variable[A2D-1]+n2; }
#define SQR(x) ((x)*(x))

//
// Utility function which returns the number of the output layer
//
int outputlayer()
{
	return nb_layers-1;
}

//
// the transfer function of the neurons (sigmoid in this case)
//
double transfer_f(double value)
{
	return( (double)(1.0 / (1.0 + exp(- (value)))) );
}

//
// the derivative of the sigmoid transfer function  (used during backpropagation)
//
double  deriv_transfer_f(double value)
{
	double tmp=transfer_f(value);
	return( (double) (tmp * (1.0 - tmp)) );
}

//
// Some Random Number Generators 
//
void  Randomise(void)
{
	/* Seed the random-number generator with current time so that
* the numbers will be different every time we run.    */
	srand( (unsigned)time( NULL ) );
}
double Random(double j) {
	return  (double) ((double)rand()/(RAND_MAX+1.0f)*j);
}

//
//Read the MLP description file which contains info about
//   number of layers
//   number of neurons per layers
//and Allocates the necessary structures
// 
void read_NNdescription(char * filename)
{
	ifstream input_file(filename);
	if(!input_file) {
		cerr << filename << " cannot be opened... Aborting" << endl;

		getchar();
		exit(-1);
	}

	//read nb layers

	input_file >> nb_layers;
	layers = (Layer *) calloc(nb_layers,sizeof(Layer));
	
	//read the number of neurons in each layer

	for (int i=0; i<nb_layers; i++) {
		input_file >> layers[i].size;
	}

	for (int i=0; i<nb_layers; i++) {
		if(i!=outputlayer()) layers[i].size+=1; //add an extra neuron for the bias
		layers[i].o = (double *) calloc(layers[i].size,sizeof(double));
		if(i) {
			Alloc2D(layers[i].w,double,(unsigned)layers[i-1].size,(unsigned)(layers[i].size));
			layers[i].u = (double *) calloc(layers[i].size,sizeof(double));
		}
	}
	input_file.close();
}

//
// Reads the NN weights from a file
// weigths are stores layer by layer from input to output
// each line correspond to weights from a single neuron to neurons in the next layer.
// the first line of each layer corresponds to the bias neuron
//
void read_NNweights(char *filename)
{
	ifstream input_file(filename);
	if(!input_file) {
		cerr << filename << " cannot be opened... Aborting" << endl;

		getchar();
		exit(-1);
	}

	//for each layer read the weights
	for (int i=1; i<nb_layers; i++) {
		for (int j=0; j<layers[i-1].size; j++) {
			// deals with the bias neuron on each layer but the outputlayer()
			int k;
			if(i!=outputlayer()) k=1;
			else k=0;
			
			for (; k<(layers[i].size); k++) {
				input_file >> layers[i].w[j][k];
			}
		}
	}
	input_file.close();
}

//
// Saves the NN weights from a file
// weigths are stores layer by layer from input to output
// each line correspond to weights from a single neuron to neurons in the next layer.
// the first line of each layer corresponds to the bias neuron
//
void save_NNweights(char *filename)
{
	ofstream OutFile(filename);
	if(!OutFile) {
		cerr << filename << " cannot be created for output... Aborting" << endl;

		getchar();
		exit(-1);
	}

	for (int i=1; i<nb_layers; i++) {
		for (int j=0; j<layers[i-1].size; j++) {
			// deals with the bias neuron on each layer but the outputlayer()
			int k;
			if(i!=outputlayer()) k=1;
			else k=0;
			
			for (; k<(layers[i].size); k++) {
				OutFile << layers[i].w[j][k] << " ";
			}
			OutFile << endl;
		}
	}

	OutFile.close();
}

//
//Read the trainning or recall patterns from a file
//and store them in the pattern array which contains 
//both input and ouput information of each examplar
//
void read_patterns(char *filename)
{
	ifstream input_file(filename);
	if(!input_file) {
		cerr << filename << " cannot be opened... Aborting" << endl;

		getchar();
		exit(-1);
	}

	//for each lines in the input file read inputs and outputs.

	if(pattern) {
		cout << "patterns already exist...freeing patterns from memory" << endl;
		Free2D(pattern);
		pattern=NULL;
	}
	
	input_file >> nb_pattern;
	Alloc2D(pattern,double,(unsigned)nb_pattern,(unsigned)((layers[0].size-1)+layers[outputlayer()].size));

	int i=0;
	int j=0;
	while(input_file.good() && !input_file.eof() && i<nb_pattern) {
		for (j=0; j<layers[0].size-1; j++) {
			input_file >> pattern[i][j];
		}
		for (int k=0; k<layers[outputlayer()].size; k++) {
			input_file >> pattern[i][j+k];
		}
		i++;
	}
	if(i!=nb_pattern) {
		cout << "Problem! the " << filename << "seems corrupted." << endl;
		cout << "Continuing with " << i << " patterns only" << endl;
		nb_pattern=i;
	}
	input_file.close();
}

//
//Computes the output of the MLP for a particular input pattern (pat)
//
void feedforward(int pat)
{

	//copy pattern into input layer
	
	for(int i=0; i < layers[0].size; i++) {
		if(i==0) layers[0].o[i] = 1.0; //bias neuron
		else layers[0].o[i] = pattern[pat][i-1]/128.0; //input normalisation required (between 0 and 1)
	}

	// propagate the data through layers 
	for(int i=1; i<nb_layers; i++) {
		
		//at each layer deal with the bias neuron (unit 0)
		int k;
		if(i!=outputlayer()) {
			k=1; 
			layers[i].o[0]=1.0; //bias neuron
		}                     
		else k=0;           // except the output layer since it doesn't have one!

		for(; k<layers[i].size; k++) {
			//
			// BEGIN IMPLEMENTATION
			//
			// compute the activation of neuron k in layer i (weighted sum)
			//initialize the value of activation
			layers[i].u[k] = 0;
			//compute activation
			for(int j=0; j<layers[i-1].size;j++)
			{
				layers[i].u[k]+=layers[i].w[j][k]*layers[i-1].o[j];
			}

			// then compute the output of neuron k in layer i (using the activation function transfer_f())
			layers[i].o[k]=transfer_f(layers[i].u[k]);
			// that's it for the feedforward procedure
			// END IMPLEMENTATION
		}
	}
}

//
//Implementation of the BackPropagation Algorithm 
// compute and performs the weigth'change according to the current input 
// pattern pat
//
void backpropagation(int pat)
{
	//compute the error and dE_du at the output layer (outputlayer())
	for(int i=0; i < layers[outputlayer()].size; i++) {
		//
		// BEGIN IMPLEMENTATION
		//
		
		// compute the difference between the network ouput and the desired output (pattern[pat][])
		// compute the local contribution to the RMS error (rms_error+=square of the difference computed above)
			rms_error = layers[outputlayer()].o[i]-pattern[pat][i+2];
		// compute the partial derivative of the error with respect to the activation at the output layer (layers[outputlayer()].dE_du[]) using the derivative of the activation function (deriv_transfer_f(layers[].u[]))
			layers[outputlayer()].dE_du[i]=2*rms_error*deriv_transfer_f(transfer_f(layers[outputlayer()].u[i]));
			                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         	
		// END IMPLEMENTATION
	}

	//compute the deltas for the remaining layers if any!
	for(int i=nb_layers-2; i>0; i--) {
		for(int j=0; j<layers[i].size; j++) {
			//
			// BEGIN IMPLEMENTATION
			//
			// compute the weighted sum of delta j in hidden layer i   
			int delta_w=0;
			for(int k=0; k<layers[i-1].size; k++)
			{
				delta_w += layers[i+1].dE_du[k]*layers[i+1].w[j][k];
			}

			// compute the partial derivative of the error with respect to the activation at the current layer (layers[i].dE_du[j]) using the derivative of the activation function (deriv_transfer_f(layers[].u[]))
			layers[i].dE_du[j]=delta_w*deriv_transfer_f(transfer_f(layers[i].u[j]));
			// END IMPLEMENTATION

		}
	}

	//compute weight change and update accordingly for all layers
	for(int i=1; i <nb_layers ; i++) {
		for(int j=0; j<layers[i-1].size; j++) {

			int k;
			if(i!=outputlayer()) k=1; //outputlayer doesn't have a bias neuron but the others do
			else k=0;

			for(; k<layers[i].size; k++) {
				//
				// BEGIN IMPLEMENTATION
				// update the weight jk in layer i using the learning rate (learning_rate) and the previously computed layer[].dE_du[] and output layer[].o[]
				layers[i].w[j][k]-=learning_rate*layers[i].dE_du[k]*layers[i].o[j];
				// END IMPLEMENTATION

			}
		}
	}
}

//
//Stores a visual representation of the classification of the 
//patterns for the current network parameters
//
void visualise_results(char* filename)
{
	Image * img= new Image("../../data/class.ppm");
	for(unsigned int i=0; i<img->getWidth(); i++) {
		for(unsigned int j=0; j<img->getWidth(); j++) {

			if((int)img->getR(i,j)==255) {
				img->put(i,j,212,128,128);
			} else {
				if(img->getG(i,j)==255) {
					img->put(i,j,128,212,128);
				} else {
					if(img->getB(i,j)==255) {
						img->put(i,j,128,128,212);
					} else {
						cout << "Invalid Pixel Value! " << img->getR(i,j) << "," << img->getG(i,j) << "," << img->getB(i,j) << endl;
					}
				}
			}
		}
	}

	for(int pat=0; pat<nb_pattern; pat++) {
		feedforward(pat);
		
		if(layers[outputlayer()].o[0]>=.5 && layers[outputlayer()].o[1]<.5 && layers[outputlayer()].o[2]<.5) 
		img->put((unsigned int)pattern[pat][0],127-(unsigned int)pattern[pat][1],255,0,0);
		else if(layers[outputlayer()].o[0]<.5 && layers[outputlayer()].o[1]>=.5 && layers[outputlayer()].o[2]<.5) 
		img->put((unsigned int)pattern[pat][0],127-(unsigned int)pattern[pat][1],0,255,0);
		else if(layers[outputlayer()].o[0]<.5 && layers[outputlayer()].o[1]<.5 && layers[outputlayer()].o[2]>=.5) 
		img->put((unsigned int)pattern[pat][0],127-(unsigned int)pattern[pat][1],0,0,255);
		else img->put((unsigned int)pattern[pat][0],127-(unsigned int)pattern[pat][1],0,0,0);

	}
	img->write(Image::color_ascii,filename);
	delete img;
	cout << "Visual Results saved in " << filename << endl;

}

//
//Stores the visual representation of the classification for all  
//input patterns for the current network parameters
//File stored in output.ppm
void visualise_output()
{
	Image * img= new Image(128,128,Image::color);
	for(int i=0;i<128;i++)
	for(int j=0;j<128;j++) {
		pattern[0][0]=i;
		pattern[0][1]=j;
		feedforward(0);

		if(layers[outputlayer()].o[0]>=.5 && layers[outputlayer()].o[1]<.5 && layers[outputlayer()].o[2]<.5) 
		img->put(i,127-j,255,0,0);
		else if(layers[outputlayer()].o[0]<.5 && layers[outputlayer()].o[1]>=.5 && layers[outputlayer()].o[2]<.5) 
		img->put(i,127-j,0,255,0);
		else if(layers[outputlayer()].o[0]<.5 && layers[outputlayer()].o[1]<.5 && layers[outputlayer()].o[2]>=.5) 
		img->put(i,127-j,0,0,255);
		else img->put(i,127-j,0,0,0);
	}
	img->write(Image::color_ascii,"output.ppm");
	delete img;
}

 
//
//
//
int main(int argc, char ** argv)
{
	if(argc!=5) {
		cout << "Usage: " << endl;
		cout << "[Training Mode] " << argv[0] << " -t neuralnet_description_file training_set neuralnet_weights_file" << endl;
		cout << "    or "  << endl;
		cout << "[Recall   Mode] " << argv[0] << " -r neuralnet_description_file recall_set   neuralnet_weights_file" << endl;

		getchar();
		exit(-1);
	}
	//for(int i=0;i<argc;i++){
	//printf(argv[i]);
	//}
	enum MODE {TRAINING,RECALL};
	MODE mode;

	if (strcmp(argv[1],"-t")==0) {
		cout << "Running in training mode." << endl;
		mode = TRAINING;
	} else if(strcmp(argv[1],"-r")==0) {
		cout << "Running in recall mode." << endl;
		mode = RECALL;
	}

	read_NNdescription(argv[2]);

	if(mode==TRAINING) {

		Randomise(); // init random seed!
		
		FILE *input_file;
			
		if(fopen_s(&input_file,argv[argc-1], "r") != 0) {
			cout << argv[argc-1] << " doesn't exists. Using random weights..." << endl;

			// Randomise the weights with small random values.
			for (int i=1; i<nb_layers; i++) {
				for(int j=0;j<layers[i-1].size;j++) {
					
					int k; // unit (neuron) 0 is used as bias except for the output layer
					if(i!=outputlayer()) k=1;
					else k=0;

					for(;k<layers[i].size;k++) {
						if(Random(1.0)>=.5)
						layers[i].w[j][k]=Random((double)(1.0/(double)(1.0*layers[i-1].size)));
						else layers[i].w[j][k]=-Random((double)(1.0/(double)(1.0*layers[i-1].size)));
					}
				}
			}
		} else {
			cout << argv[argc-1] << " exists. Loading pre-defined weights from " << argv[argc-1] << "." << endl;
			fclose(input_file);
			read_NNweights(argv[argc-1]);
		}
		
		//allocate memory for dE_du (only used in training mode)
		for (int i=1; i<nb_layers; i++) {
			layers[i].dE_du = (double *) calloc(layers[i].size,sizeof(double));
		}

		//read the training patterns
		read_patterns(argv[argc-2]);

		//perform backpropagation

		int nb_training_cycles=0;
		do {
			rms_error=0.0;
			for(int i=0;i<nb_pattern;i++) {
				feedforward(i);
				backpropagation(i);
			}

			rms_error = (double)sqrt(rms_error)/(double)(nb_pattern+layers[outputlayer()].size);
			
			nb_training_cycles++;
			if(nb_training_cycles%10==0) cout << nb_training_cycles << " RMS error=" << rms_error << endl;

		} while(rms_error>RMS_ERROR_THRESHOLD && nb_training_cycles<=MAXCYCLES);

		// save the weights to file
		cout << "saving weights to file " << argv[argc-1] << endl;
		save_NNweights(argv[argc-1]);
	} 
	else {
		cout << "Recall Mode" << endl;
		read_NNweights(argv[argc-1]);
		//perform the feed-forward operation for each pattern in the pattern list
		cout << "reading weights" << endl;
		//read the training patterns
		read_patterns(argv[argc-2]);

		double classification_error=0.0;
		for(int pat=0; pat<nb_pattern; pat++) {

			feedforward(pat);

			cout << "Pattern " << pat << "=> ";
			for(int i=1; i < layers[0].size; i++)
			cout << " " << layers[0].o[i];
			cout << " NN Out=";
			for(int i=0; i<layers[outputlayer()].size; i++) {
				cout << " " << layers[outputlayer()].o[i] << "[" << pattern[pat][i+(layers[0].size-1)] << "]";
				if(layers[outputlayer()].o[i]>=.5 && pattern[pat][i+(layers[0].size-1)]==0) {
					classification_error++;
					cout << "!";
				}
				else if(layers[outputlayer()].o[i]<.5 && pattern[pat][i+(layers[0].size-1)]==1) {
					classification_error++;
					cout << "! ";
				}
			}
			cout << endl;
		}
		cout << "Classification Errors (by output)=" << 100*classification_error/(double)(nb_pattern*layers[outputlayer()].size) << "%" << endl;
		
		visualise_results("results.ppm");
		visualise_output();
	}

	cout << "Finished. Press enter to exit";
	cin.get();
	return(1);
};
