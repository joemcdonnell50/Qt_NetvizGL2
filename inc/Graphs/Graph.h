//
// Created by werl on 25/10/16.
//

#include <stdio.h>
#include <string.h>
#include <cstdio>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <list>
#include <stack>
#include <sys/time.h>
#include "../Vertex.h"
#include "../Edge.h"

#ifndef NETVIZGL_GRAPH_H
#define NETVIZGL_GRAPH_H

using namespace std;

class Graph {
public:
    Graph();
    virtual ~Graph();

    unsigned long numVertices;
    unsigned long numEdges;

    vector<Vertex *> vertices;
    vector<Edge *> edges;
    vector<vector<int>> adjacencyMatrix;
    vector<int *> rawDataFromFile;
    vector<int *> edgeList;
    vector<int> adjList [];           //////////////////
    set<int> set;

    void update();
    void draw();
    static unsigned int hash3(unsigned int h1, unsigned int h2, unsigned int h3);

    vector<vector<int> > getAdjacencyMatrix() const;

protected:
    void createGraphData();
    void createVertices();
    void createEdges();
    void initialiseAdjacencyMatrix();
    void getEdgeListFromAdjacencyMatrix();
    void intialiseAdjacencyList();     //////////////////
    void getAdjacencyMatrixFromEdgeList();
    void setRandomColour();

//private:
    void read(char *filePath);
    int *split(string str);
};

#endif //NETVIZGL_GRAPHREADER_H
