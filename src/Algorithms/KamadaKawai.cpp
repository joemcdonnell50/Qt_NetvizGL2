#include <sys/time.h>
#include <zconf.h>
#include <QQueue>
#include "../../inc/Algorithms/KamadaKawai.h"

KamadaKawai::KamadaKawai(Graph *g) : Algorithm(g) {
    Width = 40;
    Length = 30;
    area = Width * Length;
    t = graph->numVertices;
    highT = graph->numVertices;
    k = sqrt(area / (double)t);
    w = sqrt((double)t);
    h = w;
    kConst = graph->numVertices;
    L0 = sqrt(t);
    maxiter = t * 10; //should be 10 * noNodes
    max_x = 0.0;
    max_y = 0.0;
    min_x = 0.0;
    min_y = 0.0;
    epsilon = 0.01;
    vector<vector<int>> dij(highT, vector<int>(highT));
    vector<vector<double>> kij(highT, vector<double>(highT));
    vector<vector<double>> lij(highT, vector<double>(highT));
    vector<double> D1(highT);
    vector<double> D2(highT);
    vector<double> minx(highT);
    vector<double> maxx(highT);
    vector<double> miny(highT);
    vector<double> maxy(highT);
    initialPlacement();
}

void KamadaKawai::apply() {
    Vertex *v;
    Vertex *u;
    double L;
    int dij_value = 0;

    double max_dij;


    //intialise limits
    for (int i = 0; i < graph->numVertices; ++i) {
        minx[i] = -w/2;
        maxx[i] =  w/2;
        miny[i] = -h/2;
        maxy[i] =  h/2;
    }

    //Find and intialise dij, the shortest paths between v->u
    for (int i = 0; i < graph->numVertices; ++i) {
        v = graph->vertices[i];

        for (int j = 0; j < graph->numVertices; ++j) {
            if (i == j) continue;

            u = graph->vertices[j];

            dij_value = BFS(i, j);
            //cout << "dij is " << dij_value << "on this loop " << i << endl;
            dij[i][j] = dij_value;
            //cout << "dij " << dij[i][j] << " " << i << " " << j << endl;
        }
    }

    max_dij = 0;
    //Find longest shortest path
    for (int i = 0; i < graph->numVertices; ++i) {
        for (int j = i + 1; j < graph->numVertices; ++j) {
            if (dij[i][j] > max_dij){
                max_dij = dij[i][j];
            }
        }
    }
    for (int i = 0; i < graph->numVertices; ++i) {
        for (int j = 0; j < graph->numVertices; ++j) {
            if (dij[i][j] > max_dij){
                dij[i][j] = max_dij;
            }
        }
    }

    L = L0 / max_dij;  //L0 is sqrt(noVertices)
    //compute strength of spring kij and ideal length of spring lij and fill vectors
    for (int i = 0; i < graph->numVertices; ++i) {
        for (int j = 0; j < graph->numVertices; ++j) {
            if (i == j) continue;
            double tmp = dij[i][j] * dij[i][j];
            kij[i][j] = kConst / tmp;
            lij[i][j] = L * dij[i][j]; //L is desirable length of an edge
        }
    }

    //initialise delta for x and y of node v
    for (int i = 0; i < graph->numVertices; ++i) {
        v = graph->vertices[i];
        double delta_one = 0.0, delta_two = 0.0;

        for (int j = 0; j < graph->numVertices; ++j) {
            if (i == j) continue;

            u = graph->vertices[j];
            double xDist = (v->posX - u->posX);
            double yDist = (v->posY - u->posY);
            double dist = sqrt((xDist * xDist) + (yDist * yDist));
            delta_one += kij[i][j] * (xDist - lij[i][j] * xDist / dist); //delta for x
            delta_two += kij[i][j] * (yDist - lij[i][j] * yDist / dist); //delta for y
        }
        D1[i] = delta_one;
        D2[i] = delta_two;

    }

    for(int i = 0; i < maxiter; ++i){ //Max iterations is 10 * no of nodes

        double myD1, myD2, A, B, C;
        double max_delta, delta_x, delta_y;
        double old_x, old_y, new_x, new_y;

        myD1 = 0.0, myD2 = 0.0, A=0.0, B=0.0, C=0.0;

        //Select max delta
        int m = 0; max_delta = -1;
        for (int i = 0; i < graph->numVertices; ++i) {

            double delta = sqrt((D1[i] * D1[i]) + (D2[i] * D2[i]));

            if (delta > max_delta) {
                m = i;
                max_delta = delta;
            }
        }

        //cout << max_delta << endl;



        if (max_delta < epsilon) { break; }

        //cout << "max_delta is greater than epsilon" << endl;

        //cout << m << " " << i << " times" << endl;

        v = graph->vertices[m];

        old_x = v->posX;
        old_y = v->posY;


        for (int i = 0; i < graph->numVertices; ++i){
            double dx, dy, dist, den;
            if (i == m) continue;
            v = graph->vertices[i];

            dx = old_x - v->posX; //old_x is current x-coordinate of selected node
            dy = old_y - v->posY; //old_y is current y-coordinate of selected node
            dist = sqrt(dx * dx + dy * dy);
            den = dist * (dx * dx + dy * dy);
            A += kij[m][i] * (1 - lij[m][i] * dy * dy / den);
            B += kij[m][i] * (lij[m][i] * dx * dy / den);
            C += kij[m][i] * (1 - lij[m][i] * dx * dx / den);
        }
        myD1 = D1[m]; //partial derivative of x
        myD2 = D2[m]; //partial derivative of y

        //Solve linear equations per (11) and (12) in paper
        delta_y = (B * myD1 - myD2 * A) / (C * A - B * B);
        delta_x = - (myD1 + B * delta_y) / A;

        new_x = old_x + delta_x;
        new_y = old_y + delta_y;

        //update node position
        v->posX = new_x;
        v->posY = new_y;

        //v = graph->vertices[m];

        //Limits
        if (new_x < minx[m]) { new_x = minx[m]; }
        if (new_x > maxx[m]) { new_x = maxx[m]; }
        if (new_y < miny[m]) { new_y = miny[m]; }
        if (new_y > maxy[m]) { new_y = maxy[m]; }

        //cout << minx[m] << endl;
        //cout << maxx[m] << endl;

        //                if (graph->vertices[v]->posX > max_x)
        //                  max_x = graph->vertices[v]->posX;
        //                if (graph->vertices[v]->posX < min_x)
        //                  min_x = graph->vertices[v]->posX;
        //                if (graph->vertices[v]->posY > max_y)
        //                  max_y = graph->vertices[v]->posY;
        //                if (graph->vertices[v]->posY < min_y)
        //                  min_y = graph->vertices[v]->posY;


        //Update partial derivatives for moved node
        D1[m] = D2[m] = 0.0;
        for (int i = 0; i < graph->numVertices; ++i){

            double old_dx, old_dy, old_mi, new_dx, new_dy, new_mi_dist, old_mi_dist;
            if (i == m) continue;
            v = graph->vertices[i];

            old_dx = old_x - v->posX;
            old_dy = old_y - v->posY;
            old_mi_dist = sqrt(old_dx * old_dx + old_dy * old_dy);
            new_dx = new_x - v->posX;
            new_dy = new_y - v->posY;
            new_mi_dist = sqrt(new_dx * new_dx + new_dy * new_dy);

            D1[i] -= kij[m][i] * (-old_dx + lij[m][i]) * (old_dx / old_mi_dist);

            D2[i] -= kij[m][i] * (-old_dy + lij[m][i]) * (old_dy / old_mi_dist);

            D1[i] += kij[m][i] * (-new_dx + lij[m][i]) * (new_dx / new_mi_dist);

            D2[i] += kij[m][i] * (-new_dy + lij[m][i]) * (new_dy / new_mi_dist);

            D1[m] += kij[m][i] * (new_dx - lij[m][i]) * (new_dx / new_mi_dist);

            D2[m] += kij[m][i] * (new_dy - lij[m][i]) * (new_dy / new_mi_dist);

        }
            v->posX += new_x;
            v->posY += new_y;




    }


    dij.clear();
    kij.clear();
    lij.clear();
    D1.clear();
    D2.clear();


}

int KamadaKawai::BFS(int u, int v) {
    //initialise adjacency list
    vector<int> adjList [graph->numVertices];
    for (int i = 0; i < graph->numVertices; ++i){        
        int v = graph->edgeList[i][0];
        int u = graph->edgeList[i][1];                
                adjList[v].push_back(u);
                adjList[u].push_back(v);
    }   
    int n = graph->numVertices;  // visited[n] for keeping track of visited node
    vector<bool> visited(n, 0);

    vector<int> distance(n, 1); // Initialise distances as 0
    QQueue <int> Q;  // queue to do BFS.
    distance[u] = 0;
    Q.enqueue(u);
    visited[u] = true;

    while (!Q.empty())
    {
        int x = Q.front();
        Q.dequeue();
        for (int i = 0; i < adjList[x].size(); i++)
        {
            if (visited[adjList[x][i]])
                continue;         
            distance[adjList[x][i]] = distance[x] + 1; // update distance for i
            Q.enqueue(adjList[x][i]);
            visited[adjList[x][i]] = 1;
        }
    }
    return distance[v];
//    vector<int*>::iterator iter;
//    for (iter = graph->edgeList.begin(); iter != graph->edgeList.end(); ++iter)
//        std::cout << **iter << endl;

}


void KamadaKawai::initialPlacement() {
    char *digit = new char[64];
    struct timeval time;
    gettimeofday(&time, NULL);
    srand(Graph::hash3(time.tv_sec, time.tv_usec, getpid()));
    for (int j = 0; j < graph->numVertices; ++j) {
        //sprintf(digit, "%d", j);
        graph->vertices[j]->setText(digit);
        graph->vertices[j]->posX = ((double) rand()) / RAND_MAX * (Width) - Width / 2;
        graph->vertices[j]->posY = ((double) rand()) / RAND_MAX * (Length) - Length / 2;
        graph->vertices[j]->posZ = 0;

    }

}
