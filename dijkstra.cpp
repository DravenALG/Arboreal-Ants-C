#include<iostream>
#include<algorithm>
#include<cstring>
#include <fstream>
using namespace std;
 
 

int n = 100;
int number_of_vertics = 100;
double dist[100];
bool st[100];
 
int dijkstra(double ** g){
    memset(dist, 0x3f, sizeof dist);
    dist[1] = 0;
    
    for (int i = 0; i < n; i ++){
        int t = -1;
        for (int j = 1; j <= n; j ++)
            if (!st[j] && (t == -1 || dist[t] > dist[j]))
                t = j;
                
        for (int j = 1; j <= n; j ++)
            dist[j] = min(dist[j], dist[t] + g[t][j]);
            
        st[t] = true;
    }
    
    if (dist[n] == 0x3f3f3f3f) return -1;
    return dist[n];
}

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

int main(){
	const char filename[] = "data/gnpLoc_w_loop_w_leakage_2.txt";
	double **graph = load_graph(filename);
    
    
    cout << dijkstra(graph) << endl;
    
    return 0;
}
