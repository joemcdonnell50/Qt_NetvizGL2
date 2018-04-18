#ifndef NETVIZGL_MULTILEVELGEM_H
#define NETVIZGL_MULTILEVELGEM_H

#include "Algorithm.h"
#include <sys/time.h>
#include <zconf.h>

class MultiLevelGEM : public Algorithm {
public:
    MultiLevelGEM(Graph *g);
    void apply() override;
    void initialPlacement() override;
    void placement();
    double Repulsion(double dist, int j);
    double Attraction(double dist, int j);


    double edgeIndex = 0;
    double energy1 = 0;
    double energy = 0;
    bool doPlace;
    double area;
    double k;
    double Width;
    double Length;
    double t;
    double C;
    int counter;

    vector<int> connectedNodes;
    vector<int> adjList [];
    vector<int*> newEdgeList;
    vector<int> seenVertices;
    Graph *temp;

};

#endif //NETVIZGL_MULTILEVELGEM_H
