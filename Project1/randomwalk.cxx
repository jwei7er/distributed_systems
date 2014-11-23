/*==============================================================================

 Author: Jordan Weiler
 Date:   April 24, 2014

 This program simulates the behavior of many random walks on a large graph.
 It takes an input file consisting of positive integer edges of an undirected
 graph, an output file name, and the number of times to simulate. The output
 file contains the nodes, the number of degrees for each node, and each node's
 credit for each round in the simulation.

 Each node starts with a credit of 1.0. For each subsequent round of
 simulation, the node's credit is determined by the credit of each neighbor
 divided by each neighbor's degree.

 For example, if node i has four neighbors with credit values of 0.1, 0.23,
 0.35, and 0.64 and node degrees of 5, 27, 19, and 8 respectively, the credit
 for node i after the round would be:

 c(t+1, i) = (0.1/5) + (0.23/27) + (0.35/19) + (0.64/8)
           = 0.126940

==============================================================================*/

#include <iostream>
#include <vector>
#include <time.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <cstdio>

using namespace std;

// Global variables
unsigned int maxNode = 0, totalNodes = 0;
int loops;
double credit;
double** creditsArr;
FILE *fin, *fout;
size_t i, j;

// *****************************************************************************
// Function: printTimer
//
// Purpose:
//   Print the timer and message for a section of code
//
// Arguments:
//   start - start timer
//   msg - message to send to cout with the time of execution
//
// Returns:
//   void
// *****************************************************************************
void printTimer(time_t start, string msg)
{
    time_t stop;
    time(&stop);
    cout << msg << " " << difftime(stop, start) << "sec" << endl;
}

// *****************************************************************************
// Function: readInputFile
//
// Purpose:
//   Read in the input file and creates the vector of neighbors for each node
//
// Arguments:
//   inputFile - name for input file
//
// Returns: 
//   neighbors - vector of vectors containing all neighbors for each node
// *****************************************************************************
vector< vector<int> > readInputFile(const char* inputFile)
{
    time_t startTimer;
    time(&startTimer);
    
    // Open input file and get the max node
    fin = fopen(inputFile, "r");
    if (fin != NULL)
    {
        unsigned int a, b;
        while (fscanf(fin, "%u %u", &a, &b) == 2)
        {
            if (a > maxNode) maxNode = a;
            if (b > maxNode) maxNode = b;
        }
        fclose(fin);
    }
    else
    {
        cerr << "Error! Can't open input file '" << inputFile << "'" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Create neighbors vector of vectors to the max node size
    vector< vector<int> > neighbors(maxNode);
    
    // Open input file again and add each edge to the neighbors vector
    fin = fopen(inputFile, "r");
    if (fin != NULL)
    {
        unsigned int a, b;
        while (fscanf(fin, "%u %u", &a, &b) == 2)
        {
            neighbors[a-1].push_back(b-1);
            neighbors[b-1].push_back(a-1);
        }
        fclose(fin);
    }
    else
    {
        cerr << "Error! Can't open input file '" << inputFile << "'" << endl;
        exit(EXIT_FAILURE);
    }
    
    printTimer(startTimer, "time to read input file =");
    
    return neighbors;
}

// *****************************************************************************
// Function: createOutputFile
//
// Purpose:
//   Create the output file shell with the node names and size of their
//   neighbors. This file will get read back into the program and each loop
//   will append the node's credit to the end of each node's line.
//
// Arguments:
//   neighbors - vector of vectors containing all neighbors for each node
//   outputFile - name for output file
//
// Returns: 
//   void
// *****************************************************************************
void createOutputFile(vector< vector<int> > neighbors, const char* outputFile)
{
//    time_t startTimer;
//    time(&startTimer);
    
    // Create a 2d array for old and new credits
    creditsArr = new double*[maxNode];
    for (i = 0; i < maxNode; i++)
    {
        creditsArr[i] = new double[2];
    }
    
    // Open output file and write out each node, then tab, then neighbor size
    int n = 0;
    fout = fopen(outputFile, "w");
    if (fout != NULL)
    {
        for (i = 0; i < maxNode; i++)
        {
            // Only write node out if it has neighbors (i.e. found in input file)
            // Initialize credit to 1
            if (neighbors[i].size() > 0)
            {
                creditsArr[i][0] = 1.0;
                totalNodes += 1;
                
                char buffer [50];
                n = sprintf (buffer, "%lu\t%lu\n", (i+1), neighbors[i].size());
                fwrite (buffer, sizeof(char), n, fout);
            }
            else
            {
                creditsArr[i][0] = 0;
                creditsArr[i][1] = 0;
            }
        }
        fclose(fout);
    }
    else
    {
        cerr << "Error! Can't open output file '" << outputFile << "'" << endl;
        exit(EXIT_FAILURE);
    }
    
//    printTimer(startTimer, "time to create output file =");
}


// *****************************************************************************
// Function: simulateRandomWalk
//
// Purpose:
//   Simulate the act of running several random walks. This will loop through
//   each node and calculate the node's credit based off of their neighbor's
//   credits versus their neighbor's degree. Each node's credit will get
//   appended to the end of the output file during the looping phase.
//
// Arguments:
//   neighbors - vector of vectors containing all neighbors for each node
//   outputFile - name for output file
//
// Returns: 
//   void
// *****************************************************************************
void simulateRandomWalk(vector< vector<int> > neighbors, const char* outputFile)
{
    int oldCreditIdx = 0, newCreditIdx = 1;
    double totalCredits = 0;
    unsigned int idx = 0;
    const string tempFileName = "tempOutput";
    size_t fileLineSize = 4096, slen = 0;
    time_t startTimer;
  
    // Number of times to loop through calculating credits for each node
    for (int round = 0; round < loops; round++)
    {
        time(&startTimer);
        
        char fileLine [fileLineSize];
        
        // Open output file to read previous credits
        fin = fopen(outputFile, "r");
        if (fin == NULL)
        {
            cerr << "Error! Can't open output file '" << outputFile << "'" << endl;
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
        totalCredits = 0;
        // Spin through each file line, append the new credit to the node's line, write out to temp file
        while (fgets(fileLine, fileLineSize, fin) != NULL)
        {
            if (i >= maxNode)
            {
                cerr << "Error! Counter is past size of nodes" << endl;
                cerr << "i=" << i << " max node size=" << maxNode << endl;
                exit(EXIT_FAILURE);
            }
            
            // Increment i until a node has neighbors
            while (neighbors[i].size() <= 0 && i < maxNode)
                i++;
            
            // Credit for a node is equal to each neighbor's credit / each neighbor's degree size
            credit = 0;
            for (j = 0; j < neighbors[i].size(); j++)
            {
                idx = neighbors[i][j];
                credit += creditsArr[idx][oldCreditIdx] / neighbors[idx].size();
            }
            creditsArr[i][newCreditIdx] = credit;
            
            // Strip off line break from end of line
            slen = strlen(fileLine);
            if (fileLine[slen-1] == '\n')
                fileLine[slen-1] = 0;
            
            // Write out previous file line then tab then node's credit for current loop
            fprintf(fout, "%s\t%f\n", fileLine, credit);
            //            //fprintf(fout, "%s", oss.str().c_str());
            
            //n = sprintf(fileLine, "%s\t%f\n", fileLine, credit);
            //fwrite(fileLine, sizeof(char), n, fout);
            
            i++;
            totalCredits += credit;
        }
        fclose(fin);
        fclose(fout);
        
        // Double the size for the file line if the current line len is within 100 characters
        if (slen > fileLineSize - 100)
        {
            fileLineSize *= 2;
        }
        
        // Remove last round's output file
        remove(outputFile);
        
        // Rename this round's temporary file to given output file name
        rename(tempFileName.c_str(), outputFile);
        
        // Toggle credit indices
        oldCreditIdx = newCreditIdx;
        newCreditIdx = (newCreditIdx + 1) % 2;
        
        // Print out time for round
        ostringstream oss;
        oss << "round " << (round + 1) << " =";
        printTimer(startTimer,oss.str());
        
//        // Display the credit total offset from the correct value
//        cout << "Credits: " << (totalNodes - totalCredits) << endl;
    }
}

// *****************************************************************************
// Function: validateInputArguments
//
// Purpose:
//   Validates the input arguments against necessary requirements.
//
// Arguments:
//   argc - number of command line arguments
//   argv - char array containing each argument
//
// Returns: 
//   void
// *****************************************************************************
void validateInputArguments(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "Error! Correct syntax is: prog infile outfile rounds" << endl;
        exit(EXIT_FAILURE);
    }
    
    loops = atoi(argv[3]);
    if (loops <= 0)
    {
        cerr << "Error! Number of loops must be positive integer > 0" << endl;
        exit(EXIT_FAILURE);
    }
}

// *****************************************************************************
// Function: cleanUp
//
// Purpose: Deallocates memory blocks
//
// Arguments:
//   none
//
// Returns: 
//   void
// *****************************************************************************
void cleanUp()
{
    for (i = 0; i < maxNode; i++)
    {
        if (creditsArr[i] != NULL) delete[] creditsArr[i];
    }
    if (creditsArr != NULL) delete[] creditsArr;
}


// *****************************************************************************
// Function: main
//
// Purpose:
//   Main function called when program is invoked.
//
// Arguments:
//   argc - number of command line arguments
//   argv - char array containing each argument
//
// Returns:
//   void
// *****************************************************************************
int main(int argc, char* argv[]) {
    validateInputArguments(argc, argv);

    vector< vector<int> > neighbors = readInputFile(argv[1]);

    createOutputFile(neighbors, argv[2]);

    simulateRandomWalk(neighbors, argv[2]);

    cleanUp();
    
    return 0;
}

