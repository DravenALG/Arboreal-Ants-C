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
	}
	
    return 0;
}



