#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDir>
#include <QMessageBox>
#include <QDebug>
#include <QKeyEvent>
#include <QSet>
#include <QColor>
#include <QColorDialog>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    qDebug() << "close";
    ui->widget->terminateThread();
}

void MainWindow::on_pushButton_Open_clicked()
{
    QString file = QFileDialog::getOpenFileName(this,
                                                "choose a file",
                                                "C:/Users/New User/Documents/GitHub/Qt_NetvizGL2/Graphs",
                                                "All (*.*);; Text (*.txt);; Matrix (*.mtx)");
    if(!file.isEmpty()){
        char *p = new char[file.toStdString().size()];
        strcpy(p, file.toStdString().c_str());
        ui->widget->loadGraph(p);
        ui->radioButton_Simple->setChecked(true);
    }
}

void MainWindow::on_pushButton_color_clicked()
{
    QColor color;
    color = QColorDialog::getColor(Qt::white, this);
    qDebug() << "RGB ( " << color.red() << ", " << color.green() << ", " << color.blue() << " )";
    ui->widget->setSelectedVertexColour(color.red(), color.green(), color.blue());
}

void MainWindow::on_radioButton_Simple_clicked()
{
    ui->widget->changeAlgorithm('1');
}

void MainWindow::on_radioButton_FR_clicked()
{
    ui->widget->changeAlgorithm('2');
}

void MainWindow::on_radioButton_Multiforce_clicked()
{
    ui->widget->changeAlgorithm('3');
}

void MainWindow::on_radioButton_MultiLevelGEM_clicked()
{
    ui->widget->changeAlgorithm('4');
}

void MainWindow::on_radioButton_DavidsonHarel_clicked()
{
    ui->widget->changeAlgorithm('5');
}

void MainWindow::on_radioButton_GEM_clicked()
{
    ui->widget->changeAlgorithm('6');
}

void MainWindow::on_pushButton_degree_clicked()
{
    ui->widget->degreeC();
}

void MainWindow::on_pushButton_distance_clicked()
{
    ui->widget->distanceC();
}

void MainWindow::on_pushButton_between_clicked()
{
    ui->widget->betweennessC();
}

void MainWindow::on_pushButton_save_clicked()
{
    ui->widget->saveFile();
}
