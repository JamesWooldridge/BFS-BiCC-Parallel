#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include <time.h>
#include "omp.h"
#include <assert.h>

const int MAX_BUFFER_LENGTH = 16384;

/**
 *	Reads the graph from the given file, and fills the arrays of source and dest vertices
 *	Each edge in the graph is represented by (srcVerts[i], destVerts[i])
 */
Graph * parseGraphFile(char* filename) {
	int n, m;
	int *srcVerts, *destVerts;
	FILE *file = fopen(filename, "r");
	Graph *graph;
	if(file != NULL) {
		char line[1000];
		int lineNum = 0, vertIndex = 0;
		while(fgets(line, sizeof line, file) != NULL) {
			// Strip the trailing new line character
			size_t ln = strlen(line) - 1;
			if(line[ln] == '\n') {
				line[ln] = '\0';
			}

			if(lineNum == 0) {
				// First line, get n and m values
				char *tok = strtok(line, " ");
				n = strtol(tok, NULL, 10);

				tok = strtok(NULL, " ");
				m = strtol(tok, NULL, 10);

				// Initialise srcVerts and destVerts to number of edges
				srcVerts = (int *) malloc(sizeof(int) * m);
				destVerts = (int *) malloc(sizeof(int) * m);
			} else {
				// Edge definition, put into src and dest arrays
				char *tok = strtok(line, " ");
				int src = strtol(tok, NULL, 10);
				srcVerts[vertIndex] = src;

				tok = strtok(NULL, " ");
				int dest = strtol(tok, NULL, 10);
				destVerts[vertIndex] = dest;
				vertIndex++;
			}
			lineNum++;
		}

		// Generate the graph with the edges from the file
		graph = initGraph(n, m, srcVerts, destVerts);
		free(srcVerts);
		free(destVerts);
	}

	fclose(file);
	return graph;
}

int cas(int* p, int oldval, int newval) {
    int ret = 0;
    #pragma omp critical (CAS)
    {
        int v = *p;
        if(v == oldval) {
            *p = newval;
            ret = 1;
        }
    }
    #pragma omp flush (p)
    return ret;
}

/**
 *  BFS to calculate P and L arrays of the given graph
 */
void bfs(Graph *graph, int v, int *maxDepth, int P[], int L[]) {
    // List of vertices to visit - appended to during the search
    int *toVisit = malloc(sizeof(int) * graph->numVertices);
    toVisit[0] = v;

    P[v] = v;
    L[v] = 0;

    // start and end correspond to the start and end of a level of vertices in the toVisit list
    int start = 0; int end = 1;

    #pragma omp parallel shared(start, end)
    {
        // Buffer of vertices for each thread, containing vertices to add to toVisit
        int vertexBuffer[MAX_BUFFER_LENGTH];
        
        // Initialise P and L values to -1 for our visit check later
        #pragma omp for
        for(int curr = 0; curr < graph->numVertices; curr++) {
            if(curr != v) {
                P[curr] = -1;
                L[curr] = -1;
            }
        }

        // When start == end, we've hit the end of the toVisit list, so everything is visited
        while(start != end) {
            int oldEnd = end;
            int bufferSize = 0;

            #pragma omp barrier

            #pragma omp for
            // Go through each vertex at this level
            for(int currV = start; currV < oldEnd; currV++) {
                v = toVisit[currV];

                // Get the offset into the edgeArray for v
                int *adjVertices = adjacentVertices(graph, v);
                int adjEnd = outDegree(graph, v);
                for(int currU = 0; currU < adjEnd; currU++) {
                    int u = adjVertices[currU];

                    // Set the P value using CAS in case another thread got in here too
                    if(cas(&P[u], -1, v)) {
                        L[u] = L[v] + 1;
                        if(bufferSize < MAX_BUFFER_LENGTH) {
                            // CAS successful, room left in buffer - add u to threads buffer
                            vertexBuffer[bufferSize++] = u;
                        } else {
                            // CAS successful, no room in buffer - move buffer into visit list
                            
                            // Move shared end index so we have room to add vertices from our buffer
                            int visitOffset = __sync_fetch_and_add(&end, MAX_BUFFER_LENGTH);
                            assert(visitOffset + MAX_BUFFER_LENGTH <= graph->numVertices);
                            for(int i = 0; i < MAX_BUFFER_LENGTH; i++) {
                                // Put each vertex in our buffer into our reserved space in the visit list
                                toVisit[visitOffset + i] = vertexBuffer[i];
                            }

                            // Buffer is empty, add this neighbour for next time
                            vertexBuffer[0] = u;
                            bufferSize = 1;
                        }
                    }
                }
            }

            // Check if there's things in the buffer to move into visit list
            if(bufferSize) {
                // Move shared end index so we have room to add vertices from our buffer
                int visitOffset = __sync_fetch_and_add(&end, bufferSize);
                assert(visitOffset + bufferSize <= graph->numVertices);
                for(int i = 0; i < bufferSize; i++) {
                    // Put each vertex in our buffer into our reserved space in the visit list
                    toVisit[visitOffset + i] = vertexBuffer[i];
                }
            }

            #pragma omp single 
            {
                // Move on to next level!
                start = oldEnd;
            }

            #pragma omp barrier
        }
    }

    free(toVisit);
}

/**
 *	Discovers the biconnected components in the given graph
 */
void bicc(Graph *graph) {
	// Get memory for the needed arrays
	int numVertices = graph->numVertices;
	int arrSize = sizeof(int) * numVertices;
	int *P = (int *) malloc(arrSize);
	int *L = (int *) malloc(arrSize);
	int *Art = (int *) malloc(arrSize);
	int *Low = (int *) malloc(arrSize);
	int *Par = (int *) malloc(arrSize);

	int maxDepth = 0;

	/* === 2: for all v ∈ V do === */
	#pragma omp parallel for
	for(int i = 0; i < numVertices; i++) {
		/* === 3: Art(v) ← false === */
		Art[i] = 0;
		/* === 5: Low(v) ← v === */
		Low[i] = i;
		/* === 6: Par(v) ← v === */
		Par[i] = i;
		L[i] = -1;
        P[i] = -1;
	}

	/* === 7: Select a root vertex r === */
	int root = graph->maxDegreeVertex;

	/* === 8: P, L, LQ ← BFS(G, r) === */
	bfs(graph, root, &maxDepth, P, L);

	// printf("i\tP\tL\n");
	// for(int i = 0; i < numVertices; i++) {
	// 	printf("%d\t%d\t%d\n", i, P[i], L[i]);
	// }
}

int main(int argc, char *argv[]) {
	// Check that all params are given
	if(argc != 3) {								
		fprintf(stderr, "Usage: ./bic graphFile outputFile\n");
		return 1;
	}

	// Parse the input file and generate graph
	char* filename = argv[1];
	Graph *graph = parseGraphFile(filename);

	// Find biconnected components
	bicc(graph);

	destroyGraph(graph);
}
