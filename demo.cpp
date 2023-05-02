#include <iostream>
#include <fstream>
#include <string>
#include <io.h>
#include <vector>
#include <limits>
#include<algorithm>
using namespace std;
/* Basic config */
const char filename[] = "./demo_graph01.txt"; // file to save the graph
//const char filename[] = "data/gnp_w_loop_w_leakage_0.txt"; // file to save the graph
const int number_of_vertics = 4; // numbers of the vertics in the graph
const double leakage_val = 0.1; // the leakage of each vertics
const double inc_rate = 1; // factor by which flow increases at each time step
const int max_iter = 1000;
const double start_flow = 1;
const double end_flow = 1;
const double decay = 1;
const double eps = 1e-10;


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

double* dijkstra(const char *path, int start, int n) {
	double **weight = load_graph(path);
	for (int i = 0; i < number_of_vertics; i++) {
		for (int j = 0; j < number_of_vertics; j++) {
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
    delete[] visited;
    return shortestDist;
}







vector<int> find_best_path(double **pher) {
    vector<int> path;
    path.push_back(0); // 添加第一个节点
    int current_node = 0;
    while (current_node != number_of_vertics - 1) { // 当前节点不是最后一个节点
        double max_pher = -1;
        vector<int> candidate_nodes; // 候选节点列表
        for (int i = 0; i < number_of_vertics; i++) {
            if (i != current_node && pher[current_node][i] > max_pher) {
                max_pher = pher[current_node][i];
                candidate_nodes.clear(); // 清空候选节点列表
                candidate_nodes.push_back(i); // 将当前节点加入候选节点列表
            } else if (i != current_node && pher[current_node][i] == max_pher) {
                candidate_nodes.push_back(i); // 将具有相同信息素值的节点加入候选节点列表
            }
        }
        if (candidate_nodes.empty()) { // 如果找不到候选节点，则路径无法完成
            path.clear();
            break;
        }
        int next_node = candidate_nodes[rand() % candidate_nodes.size()]; // 从候选节点列表中随机选择一个节点作为下一个节点
        path.push_back(next_node);
        current_node = next_node;
    }
    return path;
}




void getAllFiles(string path, vector<string>& files) {
    //鏂囦欢鍙ユ焺
    long hFile = 0;
    //鏂囦欢淇℃伅
    struct _finddata_t fileinfo;  
    string p;  
    if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(),&fileinfo)) != -1) {
        do {
            if ((fileinfo.attrib & _A_SUBDIR)) {  //姣旇緝鏂囦欢绫诲瀷鏄惁鏄枃浠跺す
               if (strcmp(fileinfo.name,".") != 0 && strcmp(fileinfo.name,"..") != 0) {
                   files.push_back(p.assign(path).append("\\").append(fileinfo.name));
                   getAllFiles(p.assign(path).append("\\").append(fileinfo.name), files);
               }
           } else {
               files.push_back(p.assign(path).append("\\").append(fileinfo.name));
           }
       } while (_findnext(hFile, &fileinfo) == 0);  //瀵绘壘涓嬩竴涓紝鎴愬姛杩斿洖0锛屽惁鍒?1
       _findclose(hFile);
   }
}


/* Init pheromone */
double **init_pher(double ** graph) {
	double **init_pher = new double*[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		init_pher[i] = new double[number_of_vertics];
		for(int j = 0; j < number_of_vertics; j++)
			init_pher[i][j] = graph[i][j];
	}

	return init_pher;
}

/* Init leakage */
double *init_leakage() {
	double *init_leakage = new double[number_of_vertics];
	for (int i = 0; i < number_of_vertics; i++) {
		init_leakage[i] = leakage_val;
	}
	return init_leakage;
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
		flow[number_of_vertics-1] = end_flow;
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
				sum += pher[i][j];
			}
			for(int j = 0; j < number_of_vertics; j++) {
				normalize_pher[i][j] = pher[i][j] / (sum + eps);
			}
		}
	} else {
		for(int j = 0; j < number_of_vertics; j++) {
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++) {
				sum += pher[i][j];
			}
			for (int i = 0; i < number_of_vertics; i++) {
				normalize_pher[i][j] = pher[i][j] / (sum + eps);
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
double *new_flow(double **matrix, bool is_col) {
	double *result = new double[number_of_vertics];
	if(is_col) {
		for (int i = 0; i < number_of_vertics; i++) {
			double sum = 0;
			for (int j = 0; j < number_of_vertics; j++) {
				sum += matrix[i][j];
			}
			result[i] = sum * (1 - leakage_val);
		}
	} else {
		for (int j = 0; j < number_of_vertics; j++) {
			double sum = 0;
			for (int i = 0; i < number_of_vertics; i++) {
				sum += matrix[i][j];
			}
			result[j] = sum * (1 - leakage_val);
		}
	}
	return result;
}


int decrease(double **matrix, double *vertor1, double *vertor2) {
	for (int i = 0; i < number_of_vertics; i++) {
		for (int j = 0; j < number_of_vertics; j++) {
			matrix[i][j] = matrix[i][j] / inc_rate;
		}
		vertor1[i] = vertor1[i] / inc_rate;
		vertor2[i] = vertor2[i] / inc_rate;
	}
}



int arboreal_ants(const char *path) {
	
	double **graph = load_graph(path);
	//cout << "-----------load graph, done-----------\n";
	double **pher = init_pher(graph);
	//cout << "-----------initialize the phermone, done-----------\n";
	double *leakage = init_leakage();
	//cout << "-----------initialize the leakage, done-----------\n";
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
		forward_flow = new_flow(forward_flow_matrix, true);
		backward_flow = new_flow(backward_flow_matrix, false);
		forward_flow[0] = start_flow * inc_rate;
		backward_flow[number_of_vertics-1] = end_flow * inc_rate;

		// decrease the value (just optional)
		decrease(pher, forward_flow, backward_flow);

		// next iteration
		iter++;
	}

//	cout << "-----------iter-----------\n";
//	for (int i = 0; i < number_of_vertics; i++) {
//		for (int j = 0; j < number_of_vertics; j++) {
//			cout << pher[i][j] << " ";
//		}
//		cout << "\n";
//	}
	vector<int> best_path;
	best_path = find_best_path(pher);
	//cout<<"\n @@@@@@@@@@@@@@Path Length: "<<best_path.size()<<" @@@@@@@@@@@@@@@@"<<endl;
	//cout << "######### Best Path #########\n";
	//for (int i = 0; i < best_path.size(); i++) {
    //    cout << best_path[i] << " ";
    //}
    //cout << "\n #############################\n";
    int shortest_path_length = best_path.size();


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
	
	
	cout<<"-----------Graph Path: "<<filename<<" -----------"<<endl;
	string string_path = filename;
	const char *path=string_path.c_str();
	arboreal_shortest_path_length = arboreal_ants(path);
	double* shortestDist = dijkstra(path, start_node, number_of_vertics);
	dijkstra_shortest_path_length = shortestDist[number_of_vertics-1];
	cout << "####### Arboreal Ants Shortest Path: " << arboreal_shortest_path_length << "#######\n";
	cout << "####### Dijkstra Shortest Path: " << dijkstra_shortest_path_length << "#######\n";
	cout << "-------------------------------------------\n\n";
		
	
}
