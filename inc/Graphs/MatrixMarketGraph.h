//
// Created by william on 30/01/17.
//

#ifndef NETVIZGL_MATRIXMARKETGRAPH_H
#define NETVIZGL_MATRIXMARKETGRAPH_H

#include "Graph.h"

class MatrixMarketGraph : public Graph {
public:
    MatrixMarketGraph(char *filePath);

//    void draw() override;
//    void update() override;

    ~MatrixMarketGraph() override;

private:
    void read(char *filePath) override;
};

#endif //NETVIZGL_MATRIXMARKETGRAPH_H
