#ifndef GRAPHEMBEDDING_H
#define GRAPHEMBEDDING_H

#include "Algorithm.h"
#include <sys/time.h>
#include <zconf.h>

class GraphEmbedding : public Algorithm {
 public:
    GraphEmbedding(Graph *g);
    void apply() override;
    void initialPlacement() override;
    double area;
    double k;
    double W;
    double L;
    double t;
    void placement();
    double getWeight(int i);
    double edgeIndex = 0;
    double energy = 0;
    vector<int*> newEdgeList;
    vector<int> seenVertices;

    // algorithm parameters

    int m_numberOfRounds;
    double m_minimalTemperature;
    double m_initialTemperature;
    double m_gravitationalConstant;
    double m_desiredLength;
    double m_maximalDisturbance;
    double m_rotationAngle;
    double m_oscillationAngle;
    double m_rotationSensitivity;
    double m_oscillationSensitivity;
    double m_minDistCC;

    // node data used by the algorithm

    vector<double> m_impulseX;
    vector<double> m_impulseY;
    vector<double> m_localTemperature;
    vector<double> m_skewGauge;

    // other data used by the algorithm

    double m_barycenterX;
    double m_barycenterY;
    double m_newImpulseX;
    double m_newImpulseY;
    double m_globalTemperature;
    double m_cos;
    double m_sin;
};

#endif // GRAPHEMBEDDING_H
