# Project 1

This program simulates the behavior of many random walks on a large graph. It takes an input file consisting of positive integer edges of an undirected graph, the output file name, and the number of times to simulate. Each node starts with a credit of 1.0. For each subsequent round of simulation, the node's credit is determined by the creit of each neighbor divided by each neighbor's degree.

For example, if node *i* has four neighbors with credit values of 0.1, 0.23, 0.35, and 0.64 and node degrees of 5, 27, 19, and 8 respectively, the credit for node *i* after the round would be:

`c(t+1, i) = (0.1 / 5) + (0.23 / 27) + (0.35 / 19) + (0.64 / 8) = 0.126940`