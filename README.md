# Arboreal-Ants-C

Non-official implementation of paper "Distributed algorithms from arboreal ants for the shortest path problem" in c++, testing with our own generated data.

[Paper](https://www.pnas.org/doi/abs/10.1073/pnas.2207959120) | [Official Implementation in Python](https://github.com/shivamg13/Arboreal-Ants).

## Usage ##
**Complie and run `main.cpp`.**

**Basic config:**

	const int number_of_vertics = 100; // numbers of the vertics in the graph
	const double leakage_val = 0.1; // the leakage of each vertics
	const double inc_rate = 1; // factor by which flow increases at each time step
	const int max_iter = 2000; // numbers of max iterations
	const bool bidirectional_flow = false; // whether use bidirectional flow
	const bool rand_initial = false; // whether use random initialization of pheromone in each edge
	const bool nonlinear = false; // whether use nonlinear relu