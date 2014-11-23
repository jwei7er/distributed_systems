# Project 2

The purpose of this project is to simulate the behavior of many random walks on a large graph using Message Passing Interface (MPI). This program takes an input file listing node partitions and node neighbors (degrees), an input file listing edges in the undirected graph, and an input parameter for the number of rounds to simulate. An output file for each MPI process will be created and will contain the nodes within the process, the number of neighbors for each node, and each node's credit for each round in the simulation. The output file name will be in the form of \#.out, where \# is the rank of the MPI process.

Each MPI process will calculate the credit for the nodes contained within its partition. For nodes outside of its partition, the MPI process will rely on MPI communication to receive the credits of external nodes. Each node starts with a credit of 1.0 before simultation begins. For each subsequent round of simulation, the node's credit is determined by the credit of the each neighbor of the node divided by each neighbor's degree.

For example, if node *i* has four neighbors with credit values of 0.1, 0.23, 0.35, and 0.64 and node degrees of 5, 27, 19, and 8 respectively, the credit for node *i* after the round would be:  
`c(t+1, i) = (0.1 / 5) + (0.23 / 27) + (0.35 / 19) + (0.64 / 8) = 0.126940`