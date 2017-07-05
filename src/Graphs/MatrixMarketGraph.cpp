//
// Created by william on 30/01/17.
//

#include <zconf.h>
#include "../../inc/Graphs/MatrixMarketGraph.h"
#include "../../inc/Graphs/mmio.h"
#include "../../src/Graphs/mmio.c"
#include <QDebug>

MatrixMarketGraph::MatrixMarketGraph(char *filePath) {
    read(filePath);
}

MatrixMarketGraph::~MatrixMarketGraph() {
//    fprintf(stderr, "Deleting MatrixGraph\n");
}

// Taken from the example @ http://math.nist.gov/MatrixMarket/mmio-c.html
void MatrixMarketGraph::read(char *filePath) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int rows, cols, edgs;
    int i, *I, *J;

    if ((f = fopen(filePath, "r")) == NULL) {
//        fprintf(stderr, "File not Opened");
        return;
    }

    if (mm_read_banner(f, &matcode) != 0) {
//        printf("\nCould not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    //  if (!mm_is_coordinate(matcode) || mm_is_integer(matcode)) {
    if (!mm_is_coordinate(matcode) || !mm_is_pattern(matcode) || !mm_is_symmetric(matcode)) {
//        printf("Sorry, this application only supports graphs that are:");
//        printf("Matrix Market type: [%s][%s][%s]\n", MM_COORDINATE_STR, MM_PATTERN_STR, MM_SYMM_STR);
//        printf("and not: [%s]\n", MM_INT_STR);
        exit(1);
    }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &rows, &cols, &edgs);

    if (ret_code != 0)
        exit(1);

    /* reseve memory for matrices */
    I = (int *) malloc(edgs * sizeof(int));
    J = (int *) malloc(edgs * sizeof(int));

    for (i = 0; i < edgs; i++) {
        fscanf(f, "%d %d\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    fclose(f);

//    fprintf(stdout, "\n%s\n", filePath);
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, rows, cols, edgs);
    //  for (i = 0; i < edgs; i++)
    //    fprintf(stderr, "%d %d\n", I[i], J[i]);

    numVertices = (unsigned long) rows;

    //Initialise Adj Matrix
//    for (int i = 0; i < numVertices; ++i) {
//        vector<int> row;
//        adjacencyMatrix.push_back(row);
//        for (int j = 0; j < numVertices; ++j) {
//            adjacencyMatrix[i].push_back(0);
//        }
//    }

    for (int j = 0; j < numVertices; ++j) {
        vertices.push_back(new Vertex(0, 0, 0));
        vertices[j]->setColour(0, 0, 0);
    }

    initialiseAdjacencyMatrix();

    //Attach points to each other
    for (int k = 0; k < edgs; ++k) {
        if(I[k] == J[k]){
            //fprintf(stderr, "%d %d\n", I[k], J[k]);
        }
        else {
//            vertices[I[k]]->attachPoint(vertices[J[k]]);
            edges.push_back(new Edge(vertices[I[k]], vertices[J[k]]));
            adjacencyMatrix[I[k]][J[k]] = 1;
            adjacencyMatrix[J[k]][I[k]] = 1;
            set.insert(I[k]);
            set.insert(J[k]);
        }
    }

    std::set<int>::iterator it;
    for(it = set.begin(); it != set.end(); it++){
        qDebug() << *it;
    }

    for (int i = 0; i < numVertices; ++i) {
        for (int j = i; j < numVertices; ++j) {
            if (j == i) continue;
            if (adjacencyMatrix[i][j] == 1) {
                vertices[i]->degree++;
                vertices[j]->degree++;
            }
        }
    }

    edgeList.clear();
    getEdgeListFromAdjacencyMatrix();

    numEdges = edgeList.size();
    setRandomColour();
    //  for (int i = 0; i < edgeList.size(); ++i) {
    //    fprintf(stderr, "%d,%d\n", edgeList[i][0], edgeList[i][1]);
    //  }
}

int *MatrixMarketGraph::split(string str)
{
    return NULL;
}


