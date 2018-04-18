#include "testthread.h"
#include "inc/Algorithms/SimpleForceDirected.h"
#include "inc/Algorithms/FruchtermanReingold.h"
#include "inc/Algorithms/MultiForce.h"
#include <QDebug>

TestThread::TestThread(QObject *p) : QThread(p)
{
    isPaused = false;
    algorithm = 0;
}

TestThread::~TestThread()
{
    delete algorithm;
    algorithm = 0;
}

void TestThread::run()
{
    bool finished = false;
    //int count = 0;
    while(1){
        sync.lock();
        if(isPaused)
            pauseCond.wait(&sync); // in this place, your thread will stop to execute until someone calls resume
        sync.unlock();
        if(algorithm!=NULL){  //&& !finished
            algorithm->apply();
            //count++;
            //if (count == 10)
                //finished = true;
       }
        //        usleep(5000);
    }
}

//Algorithm *TestThread::getAlgorithm(char a, Graph *g)
//{
//    if(a == '1')
//        return new SimpleForceDirected(g);
//    if(a == '2')
//        return new FruchtermanReingold(g);
//    if(a == '3')
//        return new MultiForce(g);
//    return NULL;
//}

void TestThread::resume()
{
    sync.lock();
    isPaused = false;
    sync.unlock();
    pauseCond.wakeAll();
}

void TestThread::pause()
{
    sync.lock();
    isPaused = true;
    sync.unlock();
}

void TestThread::setAlgorithm(Algorithm *a)
{
    algorithm = a;
}
