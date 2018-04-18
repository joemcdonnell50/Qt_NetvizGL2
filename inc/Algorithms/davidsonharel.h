#ifndef DAVIDSONHAREL_H
#define DAVIDSONHAREL_H

#include "Algorithm.h"
#include <sys/time.h>
#include <zconf.h>

class DavidsonHarel : public Algorithm {
public:
    DavidsonHarel(Graph *g);
    void apply() override;
    void initialPlacement() override;
    list<int> *adj;
    double area;
    double k;
    double W;
    double L;
    double t;
    void placement();
    int m_temperature;          //The temperature during the annealing process.
    double m_shrinkingFactor;   //The factor for radius.
    double m_diskRadius;        //The radius of the disk around the old position of a vertex where the new position will be.
    double m_energy;            //The current energy of the system.
    int m_numberOfIterations;   //The number of iterations per temperature step.
};
#endif // DAVIDSONHAREL_H
