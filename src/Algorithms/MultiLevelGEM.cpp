#include <sys/time.h>
#include <zconf.h>
#include <QQueue>
#include "../../inc/Algorithms/MultiLevelGEM.h"
#include "../../inc/Graphs/EdgeGraph.h"

bool complete = false;

MultiLevelGEM::MultiLevelGEM(Graph *g) : Algorithm(g) {
    Width = 40;
    Length = 30;
    area = Width * Length;
    t = graph->numVertices;
    k = sqrt(area / (double)t);
    counter = 30000;           //number of iterations
    initialPlacement();
    doPlace = true;
    complete = false;

}

void MultiLevelGEM::apply() {
    //multilevel GEM
    if (edgeIndex < graph->numEdges) {

        if (doPlace) {
            placement();
            //      doPlace = false;
        }

        energy = 999999999;
        while (energy > (10 + (seenVertices.size() * 0.1))) {
            energy = 0;
            Vertex *v;
            Vertex *u;


            for (int i = 0; i < seenVertices.size(); ++i) {
                v = graph->vertices[seenVertices[i]];
                v->forceX = 0;
                v->forceY = 0;

                for (int j = 0; j < seenVertices.size(); ++j) {
                    if (i == j) continue;

                    u = graph->vertices[seenVertices[j]];
                    double xDist = (v->posX - u->posX);
                    double yDist = (v->posY - u->posY);
                    double dist = sqrt((xDist * xDist) + (yDist * yDist));

                    if (dist < 0.00000000002) dist = 0.00000000002;

                    double repulsion = k * k / dist;
                    v->forceX += xDist / dist * repulsion;
                    v->forceY += yDist / dist * repulsion;
                }
            }

            for (int i = 0; i < edgeIndex; ++i) {

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

            for (int i = 0; i < seenVertices.size(); ++i) {
                v = graph->vertices[seenVertices[i]];
                v->posX += v->forceX * 0.005;
                v->posY += v->forceY * 0.005;

                if ((v->forceX + v->forceY) > energy)
                    energy = (v->forceX + v->forceY);
            }

        }
    } else {
        Vertex *v;
        Vertex *u;

    while (!complete && counter -->0){
        double desiredSqu = 5.0 * 5.0;

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
                double distSq = dist * dist;

                if (distSq < 0.00000000002) distSq = 0.00000000002;

                //double repulsion = k * k / dist;
                v->forceX += xDist * desiredSqu / distSq;
                v->forceY += yDist * desiredSqu / distSq;
            }
        }

        for (int i = 0; i < graph->numEdges; ++i) {

            v = graph->edges[i]->base;
            u = graph->edges[i]->connect;

            double xDist = (v->posX - u->posX);
            double yDist = (v->posY - u->posY);
            double dist = sqrt((xDist * xDist) + (yDist * yDist));
            double distSq = dist * dist;


            if (distSq < 0.00000000002) distSq = 0.00000000002;


            double attraction = dist * dist / k;

            v->forceX -= xDist * distSq / (desiredSqu * v->degree);
            v->forceY -= yDist * distSq / (desiredSqu * v->degree);

            u->forceX += xDist * distSq / (desiredSqu * u->degree);
            u->forceY += yDist * distSq / (desiredSqu * u->degree);

        }

        for (int i = 0; i < graph->numVertices; ++i) { //updates positions of vertex v
            v = graph->vertices[i];



            v->posX += v->forceX * 0.0015;
            v->posY += v->forceY * 0.0015;

            //if ((v->forceX + v->forceY) > energy)
                //energy = (v->forceX + v->forceY);

            //if ((v->forceX + v->forceY) > m_globalTemperature)
                    //m_globalTemperature = (v->forceX + v->forceY);

            //cout << "after calcs 2" << endl;
        }
        counter--;
        if (counter == 0) {complete = true;}
    }

  }
}



void MultiLevelGEM::placement() {

    double radius = 2;

    int v = graph->edgeList[edgeIndex][0];
    int a = graph->edgeList[edgeIndex][1];


    //calculate J
    double connectedEdges = 1;
    while (edgeIndex + connectedEdges < (graph->numEdges) &&
           graph->edgeList[edgeIndex + connectedEdges][0] == v) {
        connectedEdges++;
    }

    //seenVertices.clear();
    vector<int> connectedNodes;
    for (int l = 0; l < graph->numEdges; ++l) {
        if (graph->edgeList[l][0] == v) {
            connectedNodes.push_back(graph->edgeList[l][1]);
            //seenVertices.push_back(graph->edgeList[l][1]);
        }
    }

    for (int l = 0; l < graph->numEdges; ++l) {
        if (graph->edgeList[l][1] == v) {
            connectedNodes.push_back(graph->edgeList[l][0]);
            //seenVertices.push_back(graph->edgeList[l][0]);
        }
    }

    for (int i = 0; i < connectedEdges; ++i) {
        a = connectedNodes[i];

        graph->vertices[v]->posZ = 0;
        graph->vertices[a]->posZ = 0;

        //cerr << "placement v:" << v << " u:" << a << endl;

        graph->vertices[a]->posX = graph->vertices[v]->posX +
                cos((2 * M_PI * edgeIndex) / connectedNodes.size()) * radius;

        graph->vertices[a]->posY = graph->vertices[v]->posY +
                sin((2 * M_PI * edgeIndex) / connectedNodes.size()) * radius;

        edgeIndex++;
        //usleep(1000000);
    }

    bool isSeen = false;
    for (int i = 0; i < connectedNodes.size(); ++i) {
        for (int j = 0; j < seenVertices.size(); ++j) {
            if (connectedNodes[i] == seenVertices[j])
                isSeen = true;
        }
        if (!isSeen)
            seenVertices.push_back(connectedNodes[i]);
        isSeen = false;
    }

    //  GLWindow::Ins()->algorithm->graph = graph;
    //  GLWindow::Ins()->graph = graph;
    // << endl;
}

double MultiLevelGEM::Repulsion(double dist, int j){
    Vertex *u;
    u = graph->vertices[j];

    double repulse = -C * u->degree * k * k / dist;
    return repulse;

}

double MultiLevelGEM::Attraction(double dist, int j){
    Vertex *v;
    v = graph->vertices[j];
    double repulse = Repulsion(dist, j);
    double attract = (dist - k / v->degree) - repulse;
    return attract;

}



void MultiLevelGEM::initialPlacement() {

    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < graph->numVertices; ++j) {
        //sprintf(digit, "%d", j);
        //graph->vertices[j]->setText(digit);
        graph->vertices[j]->posX = ((double) rand()) / RAND_MAX * (Width) - Width / 2;
        graph->vertices[j]->posY = ((double) rand()) / RAND_MAX * (Length) - Length / 2;
        graph->vertices[j]->posZ = 0;

    }
}


