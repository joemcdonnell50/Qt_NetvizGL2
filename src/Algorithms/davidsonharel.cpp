#include <sys/time.h>
#include <zconf.h>
#include "../../inc/Algorithms/DavidsonHarel.h"

DavidsonHarel::DavidsonHarel(Graph *g) : Algorithm(g) {
    W = 40;
    L = 30;
    area = W * L;
    t = graph->numVertices;
    k = sqrt(area / (double)t);
    m_temperature = 1000;
    m_shrinkingFactor = 0.8;
    m_diskRadius = 100.0;
    adj = new list<int>[graph->numVertices];

    m_numberOfIterations = 30 * t;
    initialPlacement();
}

void DavidsonHarel::apply() {
    //calculate the repulsive forces for every node
    Vertex *v;
    Vertex *u;

        for (int i = 0; i < graph->numVertices; ++i) {
            v = graph->vertices[i];
            v->forceX = 0;
            v->forceY = 0;

            for (int j = 0; j < graph->numVertices; ++j) {
                if (i == j) continue;

                u = graph->vertices[j];
                double xDist = (v->posX - u->posX);
                double yDist = (v->posY - u->posY);
                double dist = sqrt((xDist * xDist) + (yDist * yDist));

                if (dist < 0.00000000002) dist = 0.00000000002;

                double div = (dist + 1.0)*(dist + 1.0);

                double repulsion = 1.0 / div;  // /

                v->forceX += xDist / dist * repulsion;
                v->forceY += yDist / dist * repulsion;
            }
        }

    //    double lengthSum = 0.0;
    //    m_multiplier = 2.0;
    //    m_prefEdgeLength = 0.0;
    //    for (int i = 0; i < graph->numVertices; ++i){

    //    }
        double prefEdgeLength = 100;

        for (int i = 0; i < graph->numEdges; ++i) {

            v = graph->edges[i]->base;
            u = graph->edges[i]->connect;

            double xDist = (v->posX - u->posX);
            double yDist = (v->posY - u->posY);
            double dist = sqrt((xDist * xDist) + (yDist * yDist));
            //cout << dist << endl;
            if (dist < 0.00000000002) dist = 1.0;

            double attraction = dist - prefEdgeLength;
            //attraction *= attraction;

            v->forceX -= xDist / dist * attraction;
            v->forceY -= yDist / dist * attraction;

            u->forceX += xDist / dist * attraction;
            u->forceY += yDist / dist * attraction;
        }

        for (int i = 0; i < graph->numVertices; ++i) {
            v = graph->vertices[i];
            v->posX += v->forceX * 0.0015;
            v->posY += v->forceY * 0.0015;

        }
        t *= .9;

}



void DavidsonHarel::initialPlacement() {
    //char *digit = new char[64];

    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < graph->numVertices; ++j) {
        //sprintf(digit, "%d", j);
        //graph->vertices[j]->setText(digit);
        graph->vertices[j]->posX = ((double) rand()) / RAND_MAX * (W) - W / 2;
        graph->vertices[j]->posY = ((double) rand()) / RAND_MAX * (L) - L / 2;
        graph->vertices[j]->posZ = 0;

    }
}
