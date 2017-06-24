//
// Created by werl on 11/11/16.
//

#ifndef NETVIZGL_ADJACENCYGRAPH_H
#define NETVIZGL_ADJACENCYGRAPH_H

#include "Graph.h"

class AdjacencyGraph : public Graph {
public:
    AdjacencyGraph(char *filePath);

//    virtual void draw();
//    virtual void update();

    virtual ~AdjacencyGraph();

private:
    virtual void read(char *filePath);
    virtual int *split(string str);
};

#endif //NETVIZGL_ADJACENCYGRAPH_H
