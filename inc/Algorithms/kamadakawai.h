#ifndef KAMADAKAWAI_H
#define KAMADAKAWAI_H

#include "Algorithm.h"
#include <sys/time.h>
#include <zconf.h>

class KamadaKawai : public Algorithm {
public:
    KamadaKawai(Graph *g);
    void apply() override;
    void initialPlacement() override;

    int BFS(int u, int v);
    double max_x;
    double min_x;
    double max_y;
    double min_y;
    double epsilon;
    double area;
    double k;
    double C;
    double highT;
    double Width;
    double Length;
    double w;
    double h;
    double t;
    double L0;
    double maxiter;
    double kConst;
    int dij_bfs;

    vector<vector<int>> dij;
    vector<vector<double>> kij;
    vector<vector<double>> lij;
    vector<double> D1;
    vector<double> D2;
    vector<double> minx;
    vector<double> maxx;
    vector<double> miny;
    vector<double> maxy;

};

#endif // KAMADAKAWAI_H
