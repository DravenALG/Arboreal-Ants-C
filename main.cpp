#include <iostream>
#include <fstream>
#include <string>
#include <io.h>
#include <vector>
#include <limits>
#include<algorithm>
using namespace std;
/* Basic config */
const int number_of_vertics = 100; // numbers of the vertics in the graph
const double leakage_val = 0.1; // the leakage of each vertics
const double inc_rate = 1; // factor by which flow increases at each time step
const int max_iter = 2000;
const double start_flow = 1;
const double end_flow = 1;
const double decay = 0.9;
const bool bidirectional_flow = true;
const bool rand_initial = false;
const bool nonlinear = false; 

bool has_short = true;
bool loop = false;
bool dead_end = false;


/* Load All Files */
void getAllFiles(string path, vector<string>& files) {
    long hFile = 0;
    struct _finddata_t fileinfo;  
    string p;  
    if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(),&fileinfo)) != -1) {
        do {
            if ((fileinfo.attrib & _A_SUBDIR)) {
               if (strcmp(fileinfo.name,".") != 0 && strcmp(fileinfo.name,"..") != 0) {
                   files.push_back(p.assign(path).append("\\").append(fileinfo.name));
                   getAllFiles(p.assign(path).append("\\").append(fileinfo.name), files);
               }
           } else {
               files.push_back(p.assign(path).append("\\").append(fileinfo.name));
           }
       } while (_findnext(hFile, &fileinfo) == 0);
       _findclose(hFile);
   }
}


/* Load graph */
double **load_graph(const char *path) {
	double **graph = new double*[number_of_vertics];
	
	// open file
	ifstream file;
	file.open(path);

	for (int i = 0; i < number_of_vertics; i++) {
		graph[i] = new double[number_of_vertics];
		for(int j = 0; j < number_of_vertics; j++)
			file >> graph[i][j];
	}

	return graph;
}


/* Release space */
int release(double **matrix){
	for (int i = 0; i < number_of_vertics; i++){
		delete[] matrix[i];
	}		
	delete[] matrix;
	return 0; 
} 


/* Dijkstra */
int findMinDist(double* shortestDist, bool* visited, int n) {
    int minDist = INT_MAX, minIndex = -1;
    for (int i = 0; i < n; i++) {
        if (!visited[i] && shortestDist[i] < minDist) {
            minDist = shortestDist[i];
            minIndex = i;
        }
    }
    return minIndex;
}

double dijkstra(const char *path, int start, int n) {
	double **weight = load_graph(path);
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			if (weight[i][j] == 0){
				weight[i][j] = 10000000;
			}
		}
	}
	
    double* shortestDist = new double[n];
    bool* visited = new bool[n];
    for (int i = 0; i < n; i++) {
        shortestDist[i] = weight[start][i];
        visited[i] = false;
    }
    shortestDist[start] = 0;
    visited[start] = true;
    for (int i = 1; i < n; i++) {
        int k = findMinDist(shortestDist, visited, n);
        visited[k] = true;
        for (int j = 0; j < n; j++) {
            if (!visited[j] && shortestDist[k] + weight[k][j] < shortestDist[j]) {
                shortestDist[j] = shortestDist[k] + weight[k][j];
            }
        }
    }
    double short_path_length = shortestDist[99];
    delete[] visited;
    delete[] shortestDist;
    release(weight);
    return short_path_length;
}



/* find the path based on pheromone */
vector<int> find_best_path(double **pher) {
    vector<int> path;
    path.push_back(0); 
    int current_node = 0;
    int counter = 0;
    while (current_node != number_of_vertics - 1) {
        double max_pher = -1;
        int next_node = 0;
        for (int i = 0; i < number_of_vertics; i++) {
            if (i != current_node && pher[current_node][i] > max_pher) {
                max_pher = pher[current_node][i];
                next_node = i;
            }
        }
        
        // some outputs for debug
//        cout << current_node << " " << max_pher << "\n"; 
        if(path.size() > 100){
        	has_short = false;
        	loop = true;
            break;
		}
		if(max_pher == 0){
        	has_short = false;
        	dead_end = true;
            break;			
		}
		path.push_back(next_node);
        current_node = next_node;
        counter++;
    }
    return path;
}


/* Init pheromone */
double **init_pher(double ** graph) {
	double **init_pher = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		init_pher[i] = new double[number_of_vertics];
		for(int j = 0; j < number_of_vertics; j++)
			if(rand_initial){
				init_pher[i][j] = graph[i][j] * rand() / RAND_MAX;
			}
			else{
				init_pher[i][j] = graph[i][j];
			}
	}

	return init_pher;
}



/* Init flow */
double *init_flow(bool is_forward) {
	double *flow = new double[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		flow[i] = 0;
	}
	if(is_forward) {
		flow[0] = start_flow;
	} else {
		if(bidirectional_flow){
			flow[number_of_vertics-1] = end_flow;
		}
	}
	return flow;
}


/* Normalize the pheromean*/
double **normalize(double **pher, bool is_row) {
	double **normalize_pher = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		normalize_pher[i] = new double[number_of_vertics];
	}
	if(is_row) {
		for (int i = 0; i < number_of_vertics; i++) {
			double sum = 0;
			for(int j = 0; j < number_of_vertics; j++) {
				if(nonlinear){
					sum += pher[i][j] * pher[i][j];
				}
				else{
					sum += pher[i][j];
				}
			}
			if(sum != 0){
				for(int j = 0; j < number_of_vertics; j++) {
					if(nonlinear){
						normalize_pher[i][j] = pher[i][j] * pher[i][j] / sum;
					}
					else{
						normalize_pher[i][j] = pher[i][j] / sum;
					}
				}
			}
			else{
				for(int j = 0; j < number_of_vertics; j++) {
					normalize_pher[i][j] = 0;
				}
			}
		}
	} else {
		for(int j = 0; j < number_of_vertics; j++) {
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++) {
				if(nonlinear){
					sum += pher[i][j] * pher[i][j];
				}
				else{
					sum += pher[i][j];
				}
			}
			if(sum != 0){
				for (int i = 0; i < number_of_vertics; i++) {
					if(nonlinear){
						normalize_pher[i][j] = pher[i][j] * pher[i][j] / sum;
					}
					else{
						normalize_pher[i][j] = pher[i][j] / sum;
					}
				}
			}
			else{
				for (int i = 0; i < number_of_vertics; i++) {
					normalize_pher[i][j] = 0;
				}
			}
		}
	}
	return normalize_pher;
}


/* Update the pheromean*/
double **update(double **pher, double **forward_flow_matrix, double **backward_flow_matrix) {
	for (int i = 0; i < number_of_vertics; i++) {
		for(int j = 0; j < number_of_vertics; j++) {
			pher[i][j] = decay * (pher[i][j] + forward_flow_matrix[i][j] + backward_flow_matrix[i][j]);
		}
	}
}


/* Flow matrix calculation */
double **flow_matrix(double **matrix, double *vertor, bool is_col) {
	double **result = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		result[i] = new double[number_of_vertics];
	}
	if(is_col) {
		for (int i = 0; i < number_of_vertics; i++) {
			for (int j = 0; j < number_of_vertics; j++) {
				result[i][j] = matrix[i][j] * vertor[i];
			}
		}
	} else {
		for (int j = 0; j < number_of_vertics; j++) {
			for (int i = 0; i < number_of_vertics; i++) {
				result[i][j] = matrix[i][j] * vertor[j];
			}
		}
	}
	return result;
}


/* New flow calculation */
int new_flow(double **matrix, double *flow,  bool is_row) {
	if(is_row) {
		for (int i = 0; i < number_of_vertics; i++) {
			double sum = 0;
			for (int j = 0; j < number_of_vertics; j++) {
				sum += matrix[i][j];
			}
			flow[i] = sum * (1 - leakage_val);
		}
	} else {
		for (int j = 0; j < number_of_vertics; j++) {
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++) {
				sum += matrix[i][j];
			}
			flow[j] = sum * (1 - leakage_val);
		}
	}
	return 0;
}


/* Decrease function, making the imcoming flow increase  equivalent  */
int decrease(double **matrix, double *vertor1, double *vertor2) {
	for (int i = 0; i < number_of_vertics; i++) {
		for (int j = 0; j < number_of_vertics; j++) {
			matrix[i][j] = matrix[i][j] / inc_rate;
		}
		vertor1[i] = vertor1[i] / inc_rate;
		vertor2[i] = vertor2[i] / inc_rate;
	}
	return 0; 
}




/* Main function of  arboreal ants */
int arboreal_ants(const char *path) {
	
	double **graph = load_graph(path);
	//cout << "-----------load graph, done-----------\n";
	double **pher = init_pher(graph);
	//cout << "-----------initialize the phermone, done-----------\n";
	double *forward_flow = init_flow(true);
	double *backward_flow = init_flow(false);
	//cout << "-----------initialize forward and backward flows, done-----------\n";


	for(int iter=0; iter<max_iter; iter++) {
		// normalize the pheromone
		double **norm_pher_forward =  normalize(pher, true);
		double **norm_pher_backward =  normalize(pher, false);
		

		// update pheromone
		double **forward_flow_matrix = flow_matrix(norm_pher_forward, forward_flow, true);
		double **backward_flow_matrix = flow_matrix(norm_pher_backward, backward_flow, false);
		update(pher, forward_flow_matrix, backward_flow_matrix);


		// updtae forward_flow and backward_flow
		new_flow(forward_flow_matrix, forward_flow, false);
		new_flow(backward_flow_matrix, backward_flow, true);
		forward_flow[0] = start_flow * inc_rate;
		forward_flow[number_of_vertics-1] = 0;
		if(bidirectional_flow){
			backward_flow[number_of_vertics-1] = end_flow * inc_rate;
		}
		backward_flow[0] = 0;

		// decrease the value, making the imcoming flow increase  equivalent
		decrease(pher, forward_flow, backward_flow);
		
		// release space
		release(norm_pher_forward);
		release(norm_pher_backward);
		release(forward_flow_matrix);
		release(backward_flow_matrix);
	}
	vector<int> best_path;
	best_path = find_best_path(pher);
    int shortest_path_length = best_path.size() - 1;
    
	// some outputs for debug
//	for (int i = 0; i < number_of_vertics; i++) {
//		for (int j = 0; j < number_of_vertics; j++) {
//			cout << pher[i][j] << " ";
//		}
//		cout << "\n";
//	}
//	
//	for (int i = 0; i < number_of_vertics; i++){
//		cout << forward_flow[i] << " ";
//	}
//	cout << "\n";
//	
//	for (int i = 0; i < number_of_vertics; i++){
//		cout << backward_flow[i] << " ";
//	}
//	cout << "\n";
	
	release(graph);
	release(pher);
	delete[] forward_flow;
	delete[] backward_flow;

	return shortest_path_length;
}


/* Main */
int main() {
	vector<string> files;
	char filePath[] = "data";
	getAllFiles(filePath, files);
	int size = files.size();
	int arboreal_shortest_path_length;
	int dijkstra_shortest_path_length;
	int start_node = 0;
	
	for (int i = 0; i < size; i++) {
		cout<<"-----------Graph Path: "<<files[i]<<" -----------"<<endl;
		string string_path = files[i];
		const char *path=string_path.c_str();
		arboreal_shortest_path_length = arboreal_ants(path);
		dijkstra_shortest_path_length = dijkstra(path, start_node, number_of_vertics);
		if(has_short){
			cout << "####### Arboreal Ants Shortest Path: " << arboreal_shortest_path_length << "#######\n";
			cout << "####### Dijkstra Shortest Path: " << dijkstra_shortest_path_length << " #######\n";
			cout << "-------------------------------------------\n\n";
		}	
		else{
			cout << "####### Arboreal Ants Has No Shortest Path" << " #######\n";
			if(loop){
				cout << "####### Failing Into a Loop" << " #######\n";
			}
			if(dead_end){
				cout << "####### Converge to Dead End" << " #######\n";
			}
			if(dijkstra_shortest_path_length == 10000000){
				cout << "####### Dijkstra Has No Shortest Path" << dijkstra_shortest_path_length << " #######\n";
			}
			else{
				cout << "####### Dijkstra Shortest Path: " << dijkstra_shortest_path_length << " #######\n";
			}
			cout << "-------------------------------------------\n\n";
		}
		has_short = true;
		loop = false;
		dead_end = false;
	}
}
