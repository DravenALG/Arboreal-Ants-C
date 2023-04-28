#include <iostream>
#include <limits.h>
#include <fstream>
using namespace std;

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

double* dijkstra(double **weight, int start, int n) {
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

double **load_graph(const char *path) {
	double **graph = new double*[100];
	
	// open file
	ifstream file;
	file.open(path);

	for (int i = 0; i < 100; i++) {
		graph[i] = new double[100];
		for(int j = 0; j < 100; j++)
			file >> graph[i][j];
	}

	return graph;
}


int main() {
	int n = 100;
	const char filename[] = "data/gnpLoc_w_loop_w_leakage_2.txt";
	double **graph = load_graph(filename);
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			if (graph[i][j] == 0){
				graph[i][j] = 10000000;
			}
			cout << graph[i][j] << " ";
		}	
		cout << "\n";
	}
	int start = 0;
	double* shortestDist = dijkstra(graph, start, n);
	for (int i = 0; i < n; i++) {
        cout << "shortest distance from " << start << " to " << i << " is: " << shortestDist[i] << endl;
    }
	for (int i = 0; i < n; i++) {
        delete[] graph[i];
    }
    delete[] graph;
    delete[] shortestDist;
    return 0;
	
}
