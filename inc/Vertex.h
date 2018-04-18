//
// Created by werl on 23/09/16.
//

#ifndef NETVIZGL_SPHERE_H
#define NETVIZGL_SPHERE_H

//#include "edge.h"
#include <cmath>
#include <vector>
#include <mutex>
#include <GL/glu.h>


using namespace std;
class Vertex {
private:
    GLdouble colourR;
    GLdouble colourG;
    GLdouble colourB;
    GLdouble *vertices;
    GLdouble *colours;
    GLuint *indices;
    GLint indIndex;
    GLdouble *pos;
    std::mutex mtx;



public:
    static const unsigned int rings = 12, sectors = 12;
    static constexpr double radius = 0.01;

    double temperature = 0.0;
    int degree = 0;
    int level;
    int vertexNumber = 0;
    string name;
    GLdouble posX, posY, posZ;
    GLdouble force;
    GLdouble forceX, forceY, forceZ;
    GLdouble velocityX, velocityY, velocityZ;
    bool selected;
    //    vector<Vertex *> attachedPoints;
    //    vector<Edge *> edges;

    Vertex(GLdouble offsetx, GLdouble offsety, GLdouble offsetz);
    virtual ~Vertex();

    void update();
    void draw();
    void drawText();
    void setColour(GLdouble r, GLdouble g, GLdouble b);
    GLdouble *getColour(GLdouble *colours);
    bool isPointerOver(double x, double y, int width, int height);
    double getDepth();
    void *getScreenPosition(GLdouble *pos);
    void setPositionFromScreen(double x, double y);
    //    void attachPoint(Vertex *p);

    //  FTPixmapFont *font;
    char *text;
    void setText(const char *t);

    string getName() const;
    void setName(const string &value);
};

#endif //NETVIZGL_SPHERE_H
