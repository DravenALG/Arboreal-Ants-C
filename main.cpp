#include <iostream>
#include <fstream>

using namespace std;

/* Basic config */
const char filename[] = "data/gnp_w_loop_w_leakage_0.txt"; // file to save the graph
const int number_of_vertics = 100; // numbers of the vertics in the graph
const double leakage_val = 0.1; // the leakage of each vertics
const double inc_rate = 1; // factor by which flow increases at each time step
const int max_iter = 1000;
const double start_flow = 1;
const double end_flow = 1;

/* Load graph */
double **load_graph()
{
	double **graph = new double*[number_of_vertics];
	
    // open file
    ifstream file;
    file.open(filename);
	
	for (int i = 0; i < number_of_vertics; i++)
	{
		graph[i] = new double[number_of_vertics];
		for(int j = 0; j < number_of_vertics; j++)
			file >> graph[i][j];
	}
 
	return graph;
}

/* Init pheromone */
double **init_pher(double **graph){
	double **init_pher = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++)
	{
		init_pher[i] = new double[number_of_vertics];
		for(int j = 0; j < number_of_vertics; j++)
			init_pher[i][j] = graph[i][j];
	}
 
	return init_pher;
}

/* Init leakage */
double *init_leakage(){
	double *init_leakage = new double[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		init_leakage[i] = leakage_val;
	}
	return init_leakage;
}


/* Init flow */
double *init_flow(bool is_forward){
	double *flow = new double[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		flow[i] = 0;
	}
	if(is_forward){
		flow[0] = start_flow;
	}
	else{
		flow[number_of_vertics-1] = end_flow;
	} 
	return flow;
}


/* Normalize the pheromean*/
int normalize(double **pher, bool is_row){
	double **normalize_pher = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		normalize_pher[i] = new double[number_of_vertics];
	}
	if(is_row){
		for (int i = 0; i < number_of_vertics; i++){
			double sum = 0;
			for(int j = 0; j < number_of_vertics; j++){
				sum += pher[i][j];
			}
			for(int j = 0; j < number_of_vertics; j++){
				normalize_pher[i][j] = pher[i][j] / sum;
			}
		}
	}
	else{
		for(int j = 0; j < number_of_vertics; j++){
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++){
				sum += pher[i][j];
			}
			for (int i = 0; i < number_of_vertics; i++){
				normalize_pher[i][j] = pher[i][j] / sum;
			}
		}
	}
}
	

/* Main */
int main() {
    double **graph = load_graph();
    cout << "-----------load graph, done-----------\n";
	double **pher = init_pher(graph);
	cout << "-----------initialize the phermone, done-----------\n";
	double *leakage = init_leakage();
	cout << "-----------initialize the leakage, done-----------\n";
	double *forward_flow = init_flow(true);
	double *backward_flow = init_flow(false);
	cout << "-----------initialize forward and backward flows, done-----------\n";
	
	
	for(int iter=0; iter<max_iter; iter++){
		// normalize the pheromone
		double **norm_pher_forward =  normalize(pher, true);
		double **norm_pher_backward =  normalize(pher, false);
		
		// update pheromone
		
		
		// updtae forward_flow and backward_flow
		
	}
	
    return 0;
}



