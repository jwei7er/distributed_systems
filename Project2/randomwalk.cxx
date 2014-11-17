/*==============================================================================
 
 Author: Jordan Weiler
 Date:   June 2, 2014
 
 This MPI program simulates the behavior of many random walks on a large graph.
 It takes an input file listing node partitions and degrees, an input file
 listing edges in the undirected graph, and the number of rounds to simulate.
 Each MPI process will calculate the credit for the nodes contained within its
 partition. For nodes outside of its partition, it will rely on MPI message
 passing to receive the credits of external nodes. Each node starts with a 
 credit of 1.0 before simultation begins. For each subsequent round of
 simulation, the node's credit is determined by the credit of the each neighbor
 of the node divided by each neighbor's degree.
 
 It takes an input file consisting of positive integer edges of an undirected
 graph, an output file name, and the number of times to simulate. Each node 
 starts with a credit of 1.0. For each subsequent round of simulation, the
 node's credit is determined by the credit of each neighbor divided by each
 neighbor's degree.
 
 For example, if node i has four neighbors with credit values of 0.1, 0.23,
 0.35, and 0.64 and node degrees of 5, 27, 19, and 8 respectively, the credit
 for node i after the round would be:
 
 c(t+1, i) = (0.1/5) + (0.23/27) + (0.35/19) + (0.64/8)
           = 0.126940
 
==============================================================================*/

#include <iostream>
#include <time.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <mpi.h>
#include <string>

using namespace std;

// *****************************************************************************
// Global variables
// *****************************************************************************
unsigned int maxNode = 0;
double** creditsArr;
unsigned int** partitionArr;
unsigned int** neighborArr;
FILE *fin, *fout;
size_t i, j, k;
int numRounds, rank, number_of_tasks;
string outputFileName;
const int DEGREE = 0;
const int PARTITION = 1;


// *****************************************************************************
// Method: cleanUp
//
// Purpose: Deallocate memory blocks
//
// *****************************************************************************
void cleanUp()
{
    // Credits array is size 2 * maxNode
    for (i = 0; i < 2; i++)
    {
        delete[] creditsArr[i];
    }
    delete[] creditsArr;
    
    // PartitionArr is size maxNode * 2
    // NeighborArr is size maxNode * (degree of node if node is given)
    for (i = 0; i < maxNode; i++)
    {
        if (partitionArr[i][DEGREE] > 0)
        {
            delete[] neighborArr[i];
        }
        delete[] partitionArr[i];
    }
    delete[] neighborArr;
    delete[] partitionArr;
}

// *****************************************************************************
// Method: createOutputFile
//
// Purpose:
//   Create the output file template with the node ids and size of their
//   neighbors. This file will get updated after each round of the simulation.
//   The output file name is created from the node rank and the .out extension.
//
// *****************************************************************************
void createOutputFile()
{
    // Create output file at the same time across all processes
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Create output file name
    ostringstream oss;
    oss << rank;
    outputFileName = oss.str() + string(".out");
    
    // Create a 2D array for old and new credits
    creditsArr = new double*[2];
    for (i = 0; i < 2; i++)
    {
        creditsArr[i] = new double[maxNode];
    }
    
    // Write out each node, then tab, then neighbor size to output file
    int n = 0;
    fout = fopen(outputFileName.c_str(), "w");
    if (fout != NULL)
    {
        for (i = 0; i < maxNode; i++)
        {
            // Only write out node if it is part of the current process's
            // partition and if it has neighbors (i.e. found in input file)
            if (partitionArr[i][PARTITION] == (unsigned)rank && partitionArr[i][DEGREE] > 0)
            {
                // Initialize credit of node to 1
                creditsArr[0][i] = 1.0;
                
                // Write out long unsigned node, then tab, then long unsigned
                // neighbor size
                char buffer [50];
                n = sprintf (buffer, "%lu\t%d\n", (i+1), partitionArr[i][DEGREE]);
                fwrite (buffer, sizeof(char), n, fout);
            }
        }
        fclose(fout);
    }
    else
    {
        cerr << "Error! Can't open output file '" << outputFileName << "'" << endl;
        exit(EXIT_FAILURE);
    }
}

// *****************************************************************************
// Method: getDiffTime
//
// Purpose:
//   Calculate the time difference between the start and stop timers
//
// Arguments:
//   start - the start timer
//
// *****************************************************************************
double getDiffTime(time_t start)
{
    time_t stop;
    time(&stop);
    return difftime(stop, start);
}

// *****************************************************************************
// Method: readEdgeFile
//
// Purpose:
//   Read in the edge input file and populate the array of neighbors for each
//   node
//
// Arguments:
//   edgeFile - name for the edge input file
//
// *****************************************************************************
void readEdgeFile(const char* edgeFile)
{
    // Array containing the current index of each node's neighbor array
    int* neighborIdx = new int[maxNode];
    
    // Read edge input file and add each edge to each node's neighbor array
    fin = fopen(edgeFile, "r");
    if (fin != NULL)
    {
        unsigned int a, b;
        while (fscanf(fin, "%u %u", &a, &b) == 2)
        {
            // Set B as a neighbor of A
            neighborArr[a-1][neighborIdx[a-1]] = b-1;
            neighborIdx[a-1] += 1;
            
            // Set A as a neighbor of B
            neighborArr[b-1][neighborIdx[b-1]] = a-1;
            neighborIdx[b-1] += 1;
        }
        fclose(fin);
    }
    else
    {
        cerr << "Error! Can't open input edge file '" << edgeFile << "'" << endl;
        exit(EXIT_FAILURE);
    }
    
    delete[] neighborIdx;
}

// *****************************************************************************
// Method: readPartitionFile
//
// Purpose:
//   Read in partition file. Create array of nodes, degrees, and partition ids.
//
// Arguments:
//   partitionFile - name for the input partition file
//
// *****************************************************************************
void readPartitionFile(const char* partitionFile)
{
    unsigned int nodeid, degree, partitionid;
    
    // Read in partition file to find the max node
    fin = fopen(partitionFile, "r");
    if (fin != NULL)
    {
        while (fscanf(fin, "%u %u %u", &nodeid, &degree, &partitionid) == 3)
        {
            if (nodeid > maxNode)
            {
                maxNode = nodeid;
            }
        }
        fclose(fin);
    }
    else
    {
        cerr << "Error! Can't open partition input file '" << partitionFile << "'" << endl;
        exit(EXIT_FAILURE);
    }

    // Create arrays of max node size for partitions and neighbors
    partitionArr = new unsigned int*[maxNode];
    neighborArr = new unsigned int*[maxNode];
    
    for (i = 0; i < maxNode; i++)
    {
        partitionArr[i] = new unsigned int[2];
    }

    // Read in partition file again to assign degree and partition id into the
    // partition array. Create an array for each node's neighbors based on the
    // degree of the node.
    fin = fopen(partitionFile, "r");
    if (fin != NULL)
    {
        while (fscanf(fin, "%u %u %u", &nodeid, &degree, &partitionid) == 3)
        {
            partitionArr[nodeid-1][DEGREE] = degree;
            partitionArr[nodeid-1][PARTITION] = partitionid;
            neighborArr[nodeid-1] = new unsigned int[degree];
        }
    }
    else
    {
        cerr << "Error! Can't open partition input file '" << partitionFile << "'" << endl;
        exit(EXIT_FAILURE);
    }
}

// *****************************************************************************
// Method: simulateRandomWalk
//
// Purpose:
//   Simulate the act of running several random walks. This will loop through
//   each node and calculate the node's credit based off of their neighbor's
//   credits divided by their neighbor's degree. Each node's credit will get
//   updated to the node's line in the output file during the looping phase.
//   Each process will write out the nodes (and credit) in its partition only.
//
// *****************************************************************************
void simulateRandomWalk()
{
    int oldCreditIdx = 0, newCreditIdx = 1;
    double credit;
    unsigned int idx = 0;
    size_t fileLineSize = 4096, slen = 0;
    time_t startTimer;
    MPI_Status status;
    
    // Create temporary output file name based on process rank
    ostringstream oss;
    oss << rank;
    string tempFileName = string("tempOutput.") + oss.str();
    
    // Create send and receive buffers for round credit MPI communication
    double *send_buffer = new double[maxNode];
    double *recv_buffer = new double[maxNode];
    
    // Simulate the random walk for the given number of rounds
    for (int round = 0; round < numRounds; round++)
    {
        // Start each round at the same time across all processes
        MPI_Barrier(MPI_COMM_WORLD);
        time(&startTimer);
        
        // Assign the send buffer to the local credits calculated last round
        for (i = 0; i < maxNode; i++)
        {
            if (partitionArr[i][PARTITION] == (unsigned)rank)
            {
                send_buffer[i] = creditsArr[oldCreditIdx][i];
            }
        }
        
        // Loop through all partition to partition combinations. Use MPI to
        // send and receive the previously calculated credits of last round so
        // that each partition can have the entire list of credits locally and
        // avoid sending other MPI communication to get any external credits
        // at a later time.
        for (int sendRank = 0; sendRank < number_of_tasks; sendRank++)
        {
            for (int recvRank = 0; recvRank < number_of_tasks; recvRank++)
            {
                if (sendRank != recvRank)
                {
                    if (sendRank == rank)
                    {
                        // Send local credit to a different partition
                        MPI_Send(send_buffer, maxNode, MPI_DOUBLE, recvRank, 5, MPI_COMM_WORLD);
                    }
                    else if (recvRank == rank)
                    {
                        // Receive credit from a different partition
                        MPI_Recv(recv_buffer, maxNode, MPI_DOUBLE, sendRank, 5, MPI_COMM_WORLD, &status);
                        for (k = 0; k < maxNode; k++)
                        {
                            // If the partition id of the node matches the
                            // send rank, that credit in the buffer needs to
                            // get copied to the credit array for last round.
                            if (partitionArr[k][PARTITION] == (unsigned)sendRank)
                            {
                                creditsArr[oldCreditIdx][k] = recv_buffer[k];
                            }
                        }

                    }
                }
            }
        }
        
        // Let all processes catch up from the MPI communication
        MPI_Barrier(MPI_COMM_WORLD);
        
        //** NOTE **
        //** We now have all credits across all processes and can simulate **
        
        char fileLine [fileLineSize];
        
        // Open output file to read previous credits line
        fin = fopen(outputFileName.c_str(), "r");
        if (fin == NULL)
        {
            cerr << "Error! Can't open output file '" << outputFileName << "'" << endl;
            exit(EXIT_FAILURE);
        }
        
        // Open temp output file to write new credits
        fout = fopen(tempFileName.c_str(), "w");
        if (fout == NULL)
        {
            cerr << "Error! Can't open temp output file '" << tempFileName << "'" << endl;
            exit(EXIT_FAILURE);
        }
        
        i = 0;
        
        // Spin through each file line, append the new credit to the node's
        // line and write back out to temp file
        while (fgets(fileLine, fileLineSize, fin) != NULL)
        {
            if (i >= maxNode)
            {
                cerr << "Error! Counter is past size of nodes" << endl;
                cerr << "i=" << i << " max node size=" << maxNode << endl;
                exit(EXIT_FAILURE);
            }
            
            // Increment i until the current node is in the process's partition
            // and the node has neighbors
            while ((partitionArr[i][PARTITION] != (unsigned)rank or partitionArr[i][DEGREE] <= 0) && i < maxNode)
            {
                i++;
            }
            
            // Credit for a node is equal to each neighbor's credit divided by
            // each neighbor's degree size
            credit = 0;
            for (j = 0; j < partitionArr[i][DEGREE]; j++)
            {
                idx = neighborArr[i][j];
                credit += creditsArr[oldCreditIdx][idx] / partitionArr[idx][DEGREE];
            }
            creditsArr[newCreditIdx][i] = credit;
            
            // Strip off line break from end of line
            slen = strlen(fileLine);
            if (fileLine[slen-1] == '\n')
            {
                fileLine[slen-1] = 0;
            }
            
            // Write out previous file line then tab then node's credit for
            // current round
            fprintf(fout, "%s\t%f\n", fileLine, credit);
            
            i++;
        }
        fclose(fin);
        fclose(fout);
        
        // Double the size for the file line if the current line length is
        // within 100 characters of the end of the file line
        if (slen > fileLineSize - 100)
        {
            fileLineSize *= 2;
        }
        
        // Remove last round's output file
        remove(outputFileName.c_str());
        
        // Rename this round's temporary file to given output file name
        rename(tempFileName.c_str(), outputFileName.c_str());
        
        // Toggle credit indices
        oldCreditIdx = newCreditIdx;
        newCreditIdx = (newCreditIdx + 1) % 2;
        
        // Print out time for round
        cout << " ---- time for round " << (round + 1) << ", partition " << rank << " = " << getDiffTime(startTimer) << "sec" << endl;

        // First barrier for MPI processes to get to this point
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Second barrier for std output to get written out
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Print out the overall time for all partitions for that round of
        // processing
        if (rank == 0)
        {
            cout << "total time for round " << (round + 1) << ": " << getDiffTime(startTimer) << "sec" << endl;
        }
    }
    
    // Make sure all processes delete their buffers at the same time
    MPI_Barrier(MPI_COMM_WORLD);
    
    delete[] send_buffer;
    delete[] recv_buffer;
}

// *****************************************************************************
// Method: validateInputArguments
//
// Purpose:
//   Validates the input arguments against necessary requirements
//
// Arguments:
//   argc - number of command line arguments
//   argv - char array containing each argument
//
// *****************************************************************************
void validateInputArguments(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "Error! Correct syntax is: prog Nodes2partition EdgeView numRounds" << endl;
        exit(EXIT_FAILURE);
    }
    
    numRounds = atoi(argv[3]);
    if (numRounds <= 0)
    {
        cerr << "Error! Number of rounds must be > 0, not " << argv[3] << endl;
        exit(EXIT_FAILURE);
    }
}

// *****************************************************************************
// Method: main
//
// Purpose:
//   Main function called when program is invoked
//
// Arguments:
//   argc - number of command line arguments
//   argv - char array containing each argument
//
// *****************************************************************************
int main(int argc, char* argv[])
{
    // Start up MPI
    MPI_Init(&argc, &argv);

    // extract current processes rank into the rank variable
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // extract the total number of processes into the number_of_tasks variable
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_tasks);

    validateInputArguments(argc, argv);

    // Calculate time it takes to read in input files
    time_t readFileTimer;
    time(&readFileTimer);
    
    readPartitionFile(argv[1]);

    readEdgeFile(argv[2]);

    cout << "time to read input files, partition " << rank << " = " << getDiffTime(readFileTimer) << "sec" << endl;
    
    createOutputFile();
    
    simulateRandomWalk();
    
    cleanUp();

    // Clean up MPI
    MPI_Finalize();

    return 0;
}
