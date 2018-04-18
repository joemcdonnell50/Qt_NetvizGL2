//
// Created by werl on 10/02/17.
//

#include <sys/time.h>
#include <zconf.h>
#include "../../inc/Algorithms/FruchtermanReingold.h"

FruchtermanReingold::FruchtermanReingold(Graph *g) : Algorithm(g) {
    W = 40;
    L = 30;
    area = W * L;
    t = graph->numVertices;
    k = sqrt(area / (double)t);
    initialPlacement();
}

void FruchtermanReingold::apply() {
    calculateRepulsiveForces();
    calculateAttractiveForces();
    updateNodePosition();
    t *= .9;
}

void FruchtermanReingold::calculateRepulsiveForces(){
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

            double repulsion = k * k / dist;
            v->forceX += xDist / dist * repulsion;
            v->forceY += yDist / dist * repulsion;
        }
    }
}

void FruchtermanReingold::calculateAttractiveForces(){
    Vertex *v;
    Vertex *u;
    for (int i = 0; i < graph->numEdges; ++i) {
        v = graph->edges[i]->base;
        u = graph->edges[i]->connect;

        double xDist = (v->posX - u->posX);
        double yDist = (v->posY - u->posY);
        double dist = sqrt((xDist * xDist) + (yDist * yDist));

        if (dist < 0.00000000002) dist = 0.00000000002;

        double attraction = dist * dist / k;

        v->forceX -= xDist / dist * attraction;
        v->forceY -= yDist / dist * attraction;

        u->forceX += xDist / dist * attraction;
        u->forceY += yDist / dist * attraction;
    }

}

void FruchtermanReingold::updateNodePosition(){
    Vertex *v;
    for (int i = 0; i < graph->numVertices; ++i) {
        v = graph->vertices[i];
        v->posX += v->forceX * 0.0015;
        v->posY += v->forceY * 0.0015;
    }

}

void FruchtermanReingold::initialPlacement() {
    char *digit = new char[64];

    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < graph->numVertices; ++j) {
        //sprintf(digit, "%d", j);
        graph->vertices[j]->setText(digit);
        graph->vertices[j]->posX = ((double) rand()) / RAND_MAX * (W) - W / 2;
        graph->vertices[j]->posY = ((double) rand()) / RAND_MAX * (L) - L / 2;
        graph->vertices[j]->posZ = 0;

    }
}
