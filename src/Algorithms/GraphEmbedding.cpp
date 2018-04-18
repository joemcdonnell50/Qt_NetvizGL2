#include <zconf.h>
#include <sys/time.h>
#include <iostream>
#include "../../inc/Algorithms/GraphEmbedding.h"

bool finished = false;

GraphEmbedding::GraphEmbedding(Graph *g) : Algorithm(g) {
    W = 40;
    L = 30;
    area = W * L;
    t = graph->numVertices;
    //k = sqrt(area / (double)t);
    m_numberOfRounds = 40000;                 // The maximal number of rounds per node.
    m_minimalTemperature = 0.005;            // The minimal temperature. //should 0.005
    m_initialTemperature = 12.0;            // The initial temperature.
    m_gravitationalConstant = 1.0 / 16.0;  // The gravitational constant.
    m_desiredLength = 5.0;          // The desired edge length.
    m_maximalDisturbance = 0;	         // The maximal disturbance.
    m_rotationAngle = (M_PI/3.0);       // The opening angle for rotations.
    m_oscillationAngle  = (M_PI/2.0);  // The opening angle for oscillations.
    m_rotationSensitivity = 0.01;     // The rotation sensitivity.
    m_oscillationSensitivity = 0.3;  // The oscillation sensitivity.
    m_minDistCC = 20;              // The minimal distance between connected components.
    m_globalTemperature = m_initialTemperature;
    finished = false;
    initialPlacement();

}

void GraphEmbedding::apply() {

    Vertex *v;
    int n = graph->numVertices;

    m_globalTemperature = m_initialTemperature;
    m_barycenterX = 0;
    m_barycenterY = 0;
    for(int i = 0; i < graph->numVertices; ++i) {
        v = graph->vertices[i];
        m_barycenterX += v->degree * v->posX;
        m_barycenterY += v->degree * v->posY;
    }

//    for(int i = 0; i < graph->numVertices; ++i) {
//        v = graph->vertices[i];
//        v->forceX = (m_barycenterX / n - v->posX) * m_gravitationalConstant;
//        v->forceY = (m_barycenterY / n - v->posY) * m_gravitationalConstant;
//    }
//    m_cos = cos(m_oscillationAngle / 2.0);
//    m_sin = sin(M_PI / 2 + m_rotationAngle / 2.0);





    int temperature = m_numberOfRounds;


    while (!finished && temperature -->0){    //(m_globalTemperature > m_minimalTemperature){ // && counter-- > 0
        Vertex *v;
        Vertex *u;

        double desiredSqu = m_desiredLength * m_desiredLength;

        for (int i = 0; i < graph->numVertices; ++i) {
            v = graph->vertices[i];
            v->forceX = 0;
            v->forceY = 0;
            //cout << v->forceX << endl;

            for (int j = 0; j < graph->numVertices; ++j) {
                if (i == j) continue;

                u = graph->vertices[j];
                double xDist = (v->posX - u->posX);
                double yDist = (v->posY - u->posY);
                double dist = sqrt((xDist * xDist) + (yDist * yDist));
                double distSq = dist * dist;

                if (distSq < 0.00000000002) distSq = 0.00000000002;

                double repulsion = desiredSqu / distSq;
                v->forceX += xDist * repulsion;
                v->forceY += yDist * repulsion;
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

            double nodeWeight = v->degree / 2 + 1.0;

            double attraction = distSq / (desiredSqu * nodeWeight);

            v->forceX -= xDist * attraction;
            v->forceY -= yDist * attraction;

            u->forceX += xDist * attraction;
            u->forceY += yDist * attraction;

        }

        for (int i = 0; i < graph->numVertices; ++i) { //updates positions of vertex v
            v = graph->vertices[i];
            v->posX += v->forceX * 0.0015;
            v->posY += v->forceY * 0.0015;
        }
        temperature--;      //m_globalTemperature += (v->forceX + v->forceY) / n;
        //cout << counter << endl;
        if (temperature == 0) {finished = true;}
    }


}

//idea for using oiscillations for placing node

//void updateNode(Vertex *v){

//    int n = graph->numVertices;
//    double impulseLength;

//    impulseLength = sqrt((m_newImpulseX * m_newImpulseX) + (m_newImpulseY * m_newImpulseY));
//    if(impulseLength > 0.0) {

//        // scale impulse by node temperature
//        m_newImpulseX *= v->temperature / impulseLength;
//        m_newImpulseY *= v->temperature / impulseLength;

//        // move node
//        v->posX += m_newImpulseX;
//        v->posY += m_newImpulseY;

//        // adjust barycenter
//        m_barycenterX += v->degree * v->m_newImpulseX;
//        m_barycenterY += v->degree * v->m_newImpulseY;

//        impulseLength = sqrt((m_newImpulseX * m_newImpulseX) + (m_newImpulseY * m_newImpulseY))
//                * sqrt((v->forceX * v->forceX) + (v->forceY * v->forceY));
//        if(impulseLength > 0.0) {

//            m_globalTemperature -= v->temperature / n;

//            // compute sine and cosine of angle between old and new impulse
//            double sinBeta, cosBeta;
//            sinBeta = (m_newImpulseX * v->forceX
//                       - m_newImpulseY * v->forceY)
//                    / impulseLength;
//            cosBeta = (m_newImpulseX * v->forceX
//                       + m_newImpulseY * v->forceY)
//                    / impulseLength;

//            // check for rotation
//            if(sinBeta > m_sin)
//                v->skewGauge += m_rotationSensitivity;

//            // check for oscillation
//            if(sqrt(cosBeta * cosBeta) > m_cos)
//                v->temperature *=
//                        (1 + cosBeta * m_oscillationSensitivity);

//            // cool down according to skew gauge
//            v->temperature *= (1.0 - sqrt(v->skewGauge * v->skewGauge);
//            if(v->temperature > m_initialTemperature)
//                v->temperature = m_initialTemperature;

//            // adjust global temperature
//            m_globalTemperature += v->temperature / n;
//        }

//        // save impulse
//        v->forceX = m_newImpulseX;
//        v->forceY = m_newImpulseY;
//    }



//}

double GraphEmbedding::getWeight(int i) {
    vector<vector<int>> adjMat;
    double degree = 0.0;
    adjMat = graph->getAdjacencyMatrix();
    for (int j = 0; j < graph->numVertices; ++j) {

        if (adjMat[i][j]==1)		//adj_mat is the adjacency matrix
          degree++;
    }
    degree = degree / 2.5 + 1.0;

    return degree;


}


void GraphEmbedding::initialPlacement() {
    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < graph->numVertices; ++j) {
        graph->vertices[j]->posX = ((double) rand()) / RAND_MAX * (W) - W / 2;
        graph->vertices[j]->posY = ((double) rand()) / RAND_MAX * (L) - L / 2;
    }

}
