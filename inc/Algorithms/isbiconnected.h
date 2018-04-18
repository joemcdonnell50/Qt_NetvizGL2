#ifndef ISBICONNECTED_H
#define ISBICONNECTED_H

#include "Algorithm.h"


class IsBiconnected : public Algorithm {

public:
    IsBiconnected(Graph *g);
    void apply() override;
    bool isBC();
    bool isBCUtil(int u, bool visited[], int disc[],int low[],int parent[]);
    void intialiseAdjList();
    list<int> *adj;
};

#endif // ISBICONNECTED_H
