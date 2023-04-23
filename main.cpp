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


double **matrix_scalar_multiply(double **matrix, double scalar) {
    double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		result[i] = new double[number_of_vertics];
	}
    for (int i = 0; i < number_of_vertics; i++) {
        for (int j = 0; j < number_of_vertics; j++) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }
    return result;
}

double **matrix_multiply(double **matrix1, double **matrix2) {
    double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		result[i] = new double[number_of_vertics];
	}

    for (int i = 0; i < number_of_vertics; i++) {
        for (int j = 0; j < number_of_vertics; j++) {
            double sum = 0;
            for (int k = 0; k < number_of_vertics; k++) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}



double **matrix_add(double **matrix1, double **matrix2) {
    double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		result[i] = new double[number_of_vertics];
	}

    for (int i = 0; i < number_of_vertics; i++) {
        for (int j = 0; j < number_of_vertics; j++) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}

double **transpose(double **matrix) {
    double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++){
		result[i] = new double[number_of_vertics];
	}
    for (int i = 0; i < number_of_vertics; i++) {
        for (int j = 0; j < number_of_vertics; j++) {
            result[i][j] = matrix[j][i];
        }}
    return result;
    }


double **matrix_scalar_divide(double **matrix, double scalar) {
    double **result = new double*[number_of_vertics];


    for (int i = 0; i < number_of_vertics; i++) {
        for (int j = 0; j < number_of_vertics; j++) {
            result[i][j] = matrix[i][j] / scalar;
        }
    }
    return result;
}



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
		vector<vector<double>> ff_new = matrix_multiply(ff, f_norm_pher);
        ff_new = matrix_scalar_multiply(ff_new, passage);
        bf_new = matrix_multiply(bf, b_norm_pher);
        bf_new = matrix_scalar_multiply(bf_new, passage);
        ff_new[0][0] = ff_start * inc_rate;
        bf_new[n-1][0] = bf_end * inc_rate;
        ff_new[n-1][0] = 0;
        bf_new[0][0] = 0;

        pher = matrix_add(pher, matrix_add(matrix_scalar_multiply(f_norm_pher, transpose(ff)), matrix_scalar_multiply(transpose(b_norm_pher), transpose(bf))));
        pher = matrix_scalar_multiply(pher, decay);
        ff = matrix_scalar_divide(ff_new, inc_rate);
        bf = matrix_scalar_divide(bf_new, inc_rate);
        pher = matrix_scalar_divide(pher, inc_rate);
        if (!non_lin_params.empty() && non_lin_params.find("k") != non_lin_params.end()) {
            non_lin_params["k"] = non_lin_params["k"] / inc_rate;
        }

        iter++;
		
		// updtae forward_flow and backward_flow
		
	}
	
    return 0;
}


