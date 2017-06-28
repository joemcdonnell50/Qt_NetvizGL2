//
// Created by werl on 11/11/16.
//

#include <cstdlib>
#include <fstream>
#include "../../inc/Graphs/EdgeGraph.h"
#include <QDebug>
#include <QSet>

EdgeGraph::EdgeGraph(char *filePath) : Graph(filePath) {
    read(filePath);
}


/*
//EdgeGraph::EdgeGraph(char *filePath, vector<int *> newEdgeList)
//    : Graph(filePath) {
//    edgeList = newEdgeList;

//    for (int i = 0; i < edgeList.size(); ++i) {
//        if (edgeList[i][0] > numVertices) {
//            numVertices = (unsigned long) edgeList[i][0];
//        }
//        if (edgeList[i][1] > numVertices) {
//            numVertices = (unsigned long) edgeList[i][1];
//        }
//    }
//    numVertices++;

//    createGraphData();

//    numEdges = edgeList.size();

//    //  for (int i = 0; i < edgeList.size(); ++i) {
//    //    fprintf(stderr, "%d,%d\n", edgeList[i][0], edgeList[i][1]);
//    //  }
//}
*/

/*
//EdgeGraph::EdgeGraph(vector<int *> newEdgeList) : EdgeGraph((char *) "./tempGraph", newEdgeList){
//}
*/

EdgeGraph::~EdgeGraph() {
    //    fprintf(stderr, "Deleting EdgeGraph\n");
}

void EdgeGraph::read(char *filePath) {
    string inString;
    ifstream inFile;
    inFile.open(filePath);
    if (inFile.is_open()) {
        //        fprintf(stdout, "\nOpened: %s \n", filePath);
    }
    else {
        //        fprintf(stderr, "\nFailed to open %s \n", filePath);
        exit(0);
    }

    while (!inFile.eof()) {
        getline(inFile, inString);
        if (inString.size() > 1){
            edgeList.push_back(split(inString));
        }
    }
    inFile.close();

    //    qDebug() << (set.find("3") != set.end());
    numVertices = set.size();

    createGraphData();

    numEdges = edges.size();

    //    std::set<string>::iterator itt;
    //    int ii = 0;
    //    for(itt = set.begin(); itt != set.end(); itt++){
    //        vertices[ii]->setName(*itt);
    //        const char *s = vertices[ii]->getName().c_str();
    //        qDebug() << s;
    //        ii++;
    //    }

    /*
    //  for (int i = 0; i < edgeList.size(); ++i) {
    //    fprintf(stderr, "%d,%d\n", edgeList[i][0], edgeList[i][1]);
    //  }
    */
}

int *EdgeGraph::split(string str) {
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;

    std::vector<std::string> tokens(beg, end);

    vector<string>::iterator it;
    for(it=tokens.begin();it!=tokens.end();it++){
        if((set.find(*it) == set.end())){
            set.insert(*it);
        }
    }
    int *ret = new int[tokens.size()];
    for (int i = 0; i < tokens.size(); ++i)
        ret[i] = atoi(tokens[i].c_str());

    return ret;
}

void EdgeGraph::createGraphData()
{
    vertices.clear();
    for (int j = 0; j < numVertices; ++j) {
        vertices.push_back(new Vertex(0, 0, 0));
        vertices[j]->setColour(0, 0, 0);
    }
    adjacencyMatrix.clear();
    initialiseAdjacencyMatrix();
    for (int k = 0; k < edgeList.size(); ++k) {
        if(edgeList[k][0] != edgeList[k][1])
            edges.push_back(new Edge(vertices[edgeList[k][0]], vertices[edgeList[k][1]]));
        vertices[edgeList[k][0]]->degree++;
        vertices[edgeList[k][1]]->degree++;
        adjacencyMatrix[edgeList[k][0]][edgeList[k][1]] = 1;
        adjacencyMatrix[edgeList[k][1]][edgeList[k][0]] = 1;
    }
}

//string *EdgeGraph::getEdgeList(string str){
//    std::istringstream buf(str);
//    std::istream_iterator<std::string> beg(buf), end;

//    std::vector<std::string> tokens(beg, end);

//    string *sss = new string[tokens.size()];
//    for (int i = 0; i < tokens.size(); ++i){
//        sss[i] = tokens[i];
//    }
//    return sss;
//}
