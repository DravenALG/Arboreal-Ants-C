#include <iostream>
#include <fstream>

using namespace std;

/* Basic config */
const char filename[] = "data/gnpLoc_w_loop_w_leakage_0.txt"; // file to save the graph
const int number_of_vertics = 100; // numbers of the vertics in the graph
const double leakage_val = 0.1; // the leakage of each vertics
const double inc_rate = 1; // factor by which flow increases at each time step
const int max_iter = 1000;
const double start_flow = 1;
const double end_flow = 1;
const double decay = 1;
const double epx = 1e-10;


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
double **init_pher(double ** graph){
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
double **normalize(double **pher, bool is_row){
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
				normalize_pher[i][j] = pher[i][j] / (sum + epx);
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
				normalize_pher[i][j] = pher[i][j] / (sum + epx);
			}
		}
	}
	return normalize_pher;
}


/* Update the pheromean*/
double **update(double **pher, double **forward_flow_matrix, double **backward_flow_matrix){
	for (int i = 0; i < number_of_vertics; i++){
		for(int j = 0; j < number_of_vertics; j++){
			pher[i][j] = decay * (pher[i][j] + forward_flow_matrix[i][j] + backward_flow_matrix[i][j]);
		}
	}	
}


/* Flow matrix calculation */ 
double **flow_matrix(double **matrix, double *vertor, bool is_col){
	double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		result[i] = new double[number_of_vertics];
	}
	if(is_col){
		for (int i = 0; i < number_of_vertics; i++){
			for (int j = 0; j < number_of_vertics; j++){
				result[i][j] = matrix[i][j] * vertor[i];
			}
		}
	}
	else{
		for (int j = 0; j < number_of_vertics; j++){
			for (int i = 0; i < number_of_vertics; i++){
				result[i][j] = matrix[i][j] * vertor[j];
			}
		}
	}
	return result;
}


/* New flow calculation */ 
double *new_flow(double **matrix, bool is_col){
	double *result = new double[number_of_vertics];
	if(is_col){
		for (int i = 0; i < number_of_vertics; i++){
			double sum = 0;
			for (int j = 0; j < number_of_vertics; j++){
				sum += matrix[i][j];
			}
			result[i] = sum * (1 - leakage_val);
		}
	}
	else{
		for (int j = 0; j < number_of_vertics; j++){
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++){
				sum += matrix[i][j];
			}
			result[j] = sum * (1 - leakage_val);
		}
	}
	return result;
}


int decrease(double **matrix, double *vertor1, double *vertor2){
	for (int i = 0; i < number_of_vertics; i++){
		for (int j = 0; j < number_of_vertics; j++){
			matrix[i][j] = matrix[i][j] / inc_rate;
		}
		vertor1[i] = vertor1[i] / inc_rate;
		vertor2[i] = vertor2[i] / inc_rate;
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
		double **forward_flow_matrix = flow_matrix(norm_pher_forward, forward_flow, true);
		double **backward_flow_matrix = flow_matrix(norm_pher_backward, backward_flow, false);
		update(pher, forward_flow_matrix, backward_flow_matrix);
	

		// updtae forward_flow and backward_flow
		forward_flow = new_flow(forward_flow_matrix, true);
		backward_flow = new_flow(backward_flow_matrix, false);
		forward_flow[0] = start_flow * inc_rate;
		backward_flow[number_of_vertics-1] = end_flow * inc_rate;
		
		// decrease the value (just optional)
		decrease(pher, forward_flow, backward_flow);
		
		// next iteration
		iter++;
	}
	
		for (int j = 0; j < number_of_vertics; j++){
			for (int i = 0; i < number_of_vertics; i++){
				cout << pher[i][j] << " ";
			}
			cout << "\n";
		}
	
    return 0;
}


