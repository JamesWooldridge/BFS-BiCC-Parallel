#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"

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

/**
 *	BFS to calculate P and L arrays of the given graph
 */
void bfs(Graph *graph, int v, int *maxDepth, int P[], int L[]) {
	// printf("starting from %d\n", v);
	// List of vertices to visit
	int *vlist = malloc(sizeof(int) * graph->numVertices);
	vlist[0] = v;

	P[v] = v;
	L[v] = 0;
	int start = 0; int end = 1;
	while(start != end) {
		// printf("start = %d, end = %d\n", start, end);
		int oldEnd = end;
		for(int curr = start; curr < oldEnd; curr++) {
			v = vlist[curr];
			// printf("\tv = %d\n", v);

			int *adjVertices = adjacentVertices(graph, v);
			int adjEnd = outDegree(graph, v);
			for(int o = 0; o < adjEnd; o++) {
				int u = adjVertices[o];
				if(P[u] == -1) {
					P[u] = v;
					L[u] = L[v] + 1;

					// printf("\t\tP[%d] = %d, L[%d] = %d\n", u, P[u], u, L[u]);

					vlist[end++] = u;
					// printf("\t\t\tadding %u to list at %d\n", u, end-1);
				}
			}
		}

		start = oldEnd;
	}
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

	printf("\ni\tP\tL\n");
	for(int i = 0; i < numVertices; i++) {
		printf("%d\t%d\t%d\n", i, P[i], L[i]);
	}

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