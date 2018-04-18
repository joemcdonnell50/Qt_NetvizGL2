#include "algorithmfactory.h"

AlgorithmFactory::AlgorithmFactory()
{

}

Algorithm *AlgorithmFactory::getAlgorithm(char a, Graph *g)
{
    if(a == '1')
        return new SimpleForceDirected(g);
    if(a == '2')
        return new FruchtermanReingold(g);
    if(a == '3')
        return new MultiForce(g);
    if(a == '4')
        return new MultiLevelGEM(g);
    if(a == '5')
        return new DavidsonHarel(g);
    if(a == '6')
        return new GraphEmbedding(g);
    return NULL;
}
