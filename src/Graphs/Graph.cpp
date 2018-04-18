//
// Created by werl on 25/10/16.
//

#include "../../inc/Graphs/Graph.h"
#include <stack>
#include <QDebug>

Graph::Graph() {
    numVertices = 0;
}

Graph::~Graph() {
    vertices.clear();
    adjacencyMatrix.clear();
    edgeList.clear();
    edges.clear();
    //adjList.clear(); /////////////////////
    set.clear();
}

void Graph::initialiseAdjacencyMatrix()
{
    for (int i = 0; i < numVertices; ++i) {
        vector<int> row;
        adjacencyMatrix.push_back(row);
        for (int j = 0; j < numVertices; ++j) {
            adjacencyMatrix[i].push_back(0);
        }
    }
}

void Graph::getEdgeListFromAdjacencyMatrix()
{
    int *temp = new int[2];
    for (int i = 0; i < numVertices; ++i) {
        for (int j = i; j < numVertices; ++j) {
            if (adjacencyMatrix[i][j] == 1 && i != j) {
                temp[0] = i;
                temp[1] = j;
                edgeList.push_back(temp);
                temp = new int[2];
            }
        }
    }
}

void Graph::intialiseAdjacencyList()       //////////////////////
{
    vector<int> adjList [numVertices];
    for (int i = 0; i < numVertices; ++i){
        //if (graph->edgeList[i][0] == graph->edgeList[i+1][0] )
        int v = edgeList[i][0];
        int u = edgeList[i][1];
                //cout << v << " " << u << endl;
                adjList[v].push_back(u);
                adjList[u].push_back(v);
                //cout << "made it " << endl;
    }
}


void Graph::getAdjacencyMatrixFromEdgeList()
{
    for (int k = 0; k < edgeList.size(); ++k) {
        adjacencyMatrix[edgeList[k][0]][edgeList[k][1]] = 1;
        adjacencyMatrix[edgeList[k][1]][edgeList[k][0]] = 1;
    }
}

vector<vector<int> > Graph::getAdjacencyMatrix() const
{
    return adjacencyMatrix;
}

void Graph::createGraphData()
{
    set.clear();
    for (int i = 0; i < edgeList.size(); ++i) {
        set.insert(edgeList[i][0]);
        set.insert(edgeList[i][1]);
    }
    numVertices = set.size();

    adjacencyMatrix.clear();
    initialiseAdjacencyMatrix();
    getAdjacencyMatrixFromEdgeList();
    edgeList.clear();
    getEdgeListFromAdjacencyMatrix();
    intialiseAdjacencyList();          ////////////////
    vertices.clear();
    createVertices();
    edges.clear();
    createEdges();
    numEdges = edges.size();
}

void Graph::createVertices()
{
    for (int j = 0; j < set.size(); ++j) {
        vertices.push_back(new Vertex(0, 0, 0));
        vertices[j]->setColour(0, 0, 0);
    }
}

void Graph::createEdges()
{
    for (int k = 0; k < edgeList.size(); ++k) {
        edges.push_back(new Edge(vertices[edgeList[k][0]], vertices[edgeList[k][1]]));
        vertices[edgeList[k][0]]->degree++;
        vertices[edgeList[k][1]]->degree++;
    }
}

void Graph::update() {
    for (int i = 0; i < vertices.size(); ++i) {
        vertices[i]->update();
    }
    for (int i = 0; i < edges.size(); ++i) {
        edges[i]->update();
    }
}

void Graph::draw() {
    for (int i = 0; i < vertices.size(); ++i) {
        vertices[i]->draw();
    }
    for (int i = 0; i < vertices.size(); ++i) {
        vertices[i]->drawText();
    }
    for (int i = 0; i < edges.size(); ++i) {
        edges[i]->draw();
    }
    for (int i = 0; i < edges.size(); ++i) {
        edges[i]->drawText();
    }
}

unsigned int Graph::hash3(unsigned int h1, unsigned int h2, unsigned int h3) {
    return ((h1 * 2654435789U) + h2) * 2654435789U + h3;
}

void Graph::setRandomColour(){
    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < vertices.size(); ++j) {
        vertices[j]->setColour(((double) rand() / (RAND_MAX)),
                               ((double) rand() / (RAND_MAX)),
                               ((double) rand() / (RAND_MAX)));
    }
}

void Graph::read(char *filePath) {
    string inString;
    ifstream inFile;

    inFile.open(filePath);
    if (inFile.is_open()) {
    }
    else {
        exit(0);
    }
    while (!inFile.eof()) {
        getline(inFile, inString);
        if (inString.size() > 1){
            rawDataFromFile.push_back(split(inString));
        }
    }
    inFile.close();
}

int *Graph::split(string str) {
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end);

    int *ret = new int[tokens.size()];
    for (int i = 0; i < tokens.size(); ++i){
        ret[i] = atoi(tokens[i].c_str());
    }
    return ret;
}

