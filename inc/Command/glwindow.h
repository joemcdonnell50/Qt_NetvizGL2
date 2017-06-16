#ifndef GLWINDOW_H
#define GLWINDOW_H

#include "../Graphs/Graph.h"

class GLWindow{

public:
    virtual Graph *getGraph() = 0;
    virtual char* getPath() = 0;
    virtual void setGraph(Graph *g) = 0;
    virtual void setPath(char *p) = 0;
};


#endif // GLWINDOW_H
