//
// Created by werl on 01/11/16.
//

#ifndef NETVIZGL_SIMPLEFORCEDIRECTED_H
#define NETVIZGL_SIMPLEFORCEDIRECTED_H

#include "../Graphs/Graph.h"
#include "Algorithm.h"

class SimpleForceDirected : public Algorithm {
//private:
//    Graph *graph;
public:
    SimpleForceDirected(Graph *g);

    void apply() override;
    void initialPlacement() override;
    void calculateForces();
    void updateNodePosition();
};

#endif //NETVIZGL_SIMPLEFORCEDIRECTED_H
