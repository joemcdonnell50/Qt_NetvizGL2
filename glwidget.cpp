#include "glwidget.h"
#include <QDebug>
#include <GL/glu.h>
#include <QTimer>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include <QCloseEvent>
#include <QPixmap>
#include <QtWidgets>
#include <QFile>
#include <QTextStream>
#include "inc/Graphs/EdgeGraph.h"
#include "inc/Graphs/AdjacencyGraph.h"
#include "inc/Graphs/MatrixMarketGraph.h"
#include "inc/Algorithms/Algorithm.h"
#include "inc/Algorithms/FruchtermanReingold.h"
#include "inc/Algorithms/SimpleForceDirected.h"
#include "inc/Algorithms/MultiForce.h"
#include "inc/Command/LoadGraph.h"
#include "inc/Command/SelectVertex.h"
#include "inc/Command/DragVertex.h"
#include "inc/Command/DeleteVertex.h"
#include "inc/Centrality/DegreeCentrality.h"
#include "inc/Centrality/DistanceCentrality.h"
#include "inc/Centrality/Betweenness.h"
#include "testthread.h"
#include "algorithmfactory.h"


GLWidget::GLWidget(QWidget *parent):QOpenGLWidget(parent)
{
    setFocus();

    translateX = 0;
    translateY = 0;
    translateZ = 1.8;
    isMouseLeftDown = false;
    isMouseMiddleDown = false;
    isMouseRightDown = false;
    isKeyCtrlDown = false;

    c = new LoadGraph(this);
    sv = new SelectVertex(this);
    dragv = new DragVertex(this);

    t = NULL;
    graph = NULL;
    selectedNode = NULL;

    // update timer
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(update()));
    timer->start(20);

}

GLWidget::~GLWidget()
{
    delete timer;
    delete t;
    delete graph;
}

void GLWidget::initializeGL()
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glClearColor(1, 1 ,1 ,1);
}

void GLWidget::paintGL()
{
    /*
     * ---draw here---
     */

    glViewport(0, 0, width(), height());

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluPerspective(45, (double)width()/(double)height(), 0.01, 10);
    glTranslatef(translateX, translateY, -translateZ);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

/*
//    glEnable(GL_LINE_SMOOTH);

//    glPointSize(30/(translateZ+0.84));
//    glLineWidth(1);

//    glVertexPointer(3, GL_FLOAT, 0, vertex);
//    glDrawArrays(GL_POINTS, 0, 1);
//    glDrawArrays(GL_LINE_LOOP, 0, 3);
*/

    /*
     * draw graph
     */
    if(graph!=NULL){
        graph->update();
        graph->draw();
    }

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);

    /*
     * ---end of drawing---
     */
}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    mouseX = event->pos().x();
    mouseY = event->pos().y();
    if(event->button() & Qt::LeftButton){
        isMouseLeftDown = true;
    }
    if(event->button() & Qt::MiddleButton){
        isMouseMiddleDown = true;
    }
    //Error on right-click is here
    if(event->button() & Qt::RightButton){
        glLoadIdentity();
        gluPerspective(45, (double)width()/(double)height(), 0.01, 10);
        glTranslatef(translateX, translateY, -translateZ);
        glViewport(0, 0, width(), height());
        if(graph != NULL){
            sv->execute();
        }

        isMouseRightDown = true;
    }
    update();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if(event->button() & Qt::LeftButton){
        isMouseLeftDown = false;
    }
    if(event->button() & Qt::MiddleButton){
        isMouseMiddleDown = false;
    }
    if(event->button() & Qt::RightButton){
        isMouseRightDown = false;
    }

    mouseX = event->pos().x();
    mouseY = event->pos().y();
    if(t != NULL){
        t->resume();
    }
    update();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    mouseDiffX = mouseX - event->pos().x();
    mouseDiffY = mouseY - event->pos().y();
/*
//    if (wind->mouseMIDDLE) {
//        wind->yaw += (xpos - wind->mouseX) / 8;
//        wind->pitch += (ypos - wind->mouseY) / 8;
//        wind->mouseX = xpos;
//        wind->mouseY = ypos;
//    }
*/
    if (isMouseLeftDown && !isKeyCtrlDown) {
        translateX += (((double)event->pos().x() - mouseX) / (double)height()) * 0.82 * translateZ;
        translateY += ((mouseY - (double)event->pos().y()) / (double)height()) * 0.82 * translateZ;
    }
    if (isMouseLeftDown && isKeyCtrlDown) {
        t->pause();
        dragv->execute();
    }

    mouseX = event->pos().x();
    mouseY = event->pos().y();
    update();
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    translateZ += ((double)event->delta() / 20 /120) * translateZ;
    if (translateZ < 0.16) {
        translateZ = 0.16;
    }
    if(translateZ > 10){
        translateZ = 10;
    }
    update();
}

void GLWidget:: keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_Delete){
        if(selectedNode){
            deletev = new DeleteVertex(this, selectedVertexNumber);
            deletev->execute();
            selectedNode = 0;
            selectedVertexNumber = -1;
            changeAlgorithm(algorithm);
        }
    }
    if(event->key() == Qt::Key_Escape){
        translateX = 0;
        translateY = 0;
    }
    if (event->key() == Qt::Key_Control){
        isKeyCtrlDown = true;
    }
    if(event->key() == Qt::Key_F11){
        saveScreenshot();
    }
    update();
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_A){

    }
    if (event->key() == Qt::Key_Control){
        isKeyCtrlDown = false;
    }
    update();
}

void GLWidget::changeAlgorithm(char a)
{
    terminateThread();
    algorithm = a;
    if (graph != NULL){
        t = new TestThread(this);
        t->setAlgorithm(AlgorithmFactory::getAlgorithm(a, graph));
//        t->getAlgorithm(a, graph);
        t->start();
    }
    setFocus();
}

void GLWidget::terminateThread()
{
    if(t!=NULL){
        t->terminate();
        t->wait();
        delete t;
        t = NULL;
    }
}

void GLWidget::setSelectedVertexColour(int r, int g, int b)
{
    GLdouble red = (GLdouble)r / 255;
    GLdouble green = (GLdouble)g / 255;
    GLdouble blue = (GLdouble)b / 255;
    if(selectedNode)
        selectedNode->setColour(red, green, blue);
    setFocus();
}

int GLWidget::getWidth()
{
    return width();
}

int GLWidget::getHeight()
{
    return height();
}

void GLWidget::loadGraph(char *p)
{
    terminateThread();
    path = p;
    c->execute();
    changeAlgorithm('1');
    setFocus();
}

void GLWidget::degreeC()
{
    DegreeCentrality cc;
    if (graph != NULL)
        cc.calcApply(graph);
    setFocus();
}

void GLWidget::distanceC()
{
    DistanceCentrality cc;
    if (graph != NULL)
        cc.calcApply(graph);
    setFocus();
}

void GLWidget::betweennessC()
{
    Betweenness cc;
    if (graph != NULL)
        cc.calcApply(graph);
    setFocus();
}

void GLWidget::saveScreenshot()
{
    originalPixmap = this->grab();
    const QString format = "png";
    QString initialPath = QStandardPaths::writableLocation(QStandardPaths::PicturesLocation);
    if (initialPath.isEmpty())
        initialPath = QDir::currentPath();
    initialPath += tr("/untitled.") + format;

    QFileDialog fileDialog(this, tr("Save As"), initialPath);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setDirectory(initialPath);
    QStringList mimeTypes;
    foreach (const QByteArray &bf, QImageWriter::supportedMimeTypes())
        mimeTypes.append(QLatin1String(bf));
    fileDialog.setMimeTypeFilters(mimeTypes);
    fileDialog.selectMimeTypeFilter("image/" + format);
    fileDialog.setDefaultSuffix(format);
    if (fileDialog.exec() != QDialog::Accepted)
        return;
    const QString fileName = fileDialog.selectedFiles().first();
    if (!originalPixmap.save(fileName)) {
        QMessageBox::warning(this, tr("Save Error"), tr("The image could not be saved to \"%1\".")
                             .arg(QDir::toNativeSeparators(fileName)));
    }
    setFocus();
}

void GLWidget::saveFile()
{
    const QString format = "txt";
    QString initialPath = QStandardPaths::writableLocation(QStandardPaths::PicturesLocation);
    if (initialPath.isEmpty())
        initialPath = QDir::currentPath();
    initialPath += tr("/untitled.") + format;

    QFileDialog fileDialog(this, tr("Save As"), initialPath);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setDirectory(initialPath);
    fileDialog.setDefaultSuffix(format);
    if (fileDialog.exec() != QDialog::Accepted)
        return;
    const QString fileName = fileDialog.selectedFiles().first();
    QFile file(fileName);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    for(int i = 0; i < graph->edgeList.size(); i++){
        QString str = QString::fromStdString(std::to_string(graph->edgeList[i][0]) + " " + std::to_string(graph->edgeList[i][1]) + "\n");
        out << str;
    }
    file.flush();
    file.close();
    setFocus();
}
