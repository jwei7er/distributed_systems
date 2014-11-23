# Project 1

This program simulates the behavior of many random walks on a large graph. It takes an input file consisting of positive integer edges of an undirected graph, an output file name, and the number of times to simulate. The output file contains the nodes, the number of degrees for each node, and each node's credit for each round in the simulation.

For example, if you wanted to run this program using exampleedges.txt as the input edge file, output.txt as the output file name, and simulate for 5 rounds, you would type:

	./randomwalk ./exampleedges.txt ./output.txt 5

Each node starts with a credit of 1.0. For each subsequent round of simulation, the node's credit is determined by the creit of each neighbor divided by each neighbor's degree.

For example, if node *i* has four neighbors with credit values of 0.1, 0.23, 0.35, and 0.64 and node degrees of 5, 27, 19, and 8 respectively, the credit for node *i* after the round would be:  
`c(t+1, i) = (0.1 / 5) + (0.23 / 27) + (0.35 / 19) + (0.64 / 8) = 0.126940`