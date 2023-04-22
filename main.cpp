#include <iostream>
#include <fstream>

using namespace std;

/* Basic config */
const char filename[] = "data/gnp_w_loop_w_leakage_0.txt";
const int number_of_vertics = 100;

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

/* Main */
int main() {
    double **graph = load_graph();
    cout << "-----------Validation of the graph-----------\n";
    for (int i = 0; i < number_of_vertics; i++)
	{
		for(int j = 0; j < number_of_vertics; j++){
			cout << graph[i][j];
		}
		cout << "\n";
	}
    return 0;
}



