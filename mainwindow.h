#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void closeEvent(QCloseEvent *event) override;

private slots:
    void on_pushButton_Open_clicked();

    void on_pushButton_color_clicked();

    void on_radioButton_Simple_clicked();

    void on_radioButton_FR_clicked();

    void on_radioButton_Multiforce_clicked();

    void on_radioButton_MultiLevelGEM_clicked();

    void on_radioButton_DavidsonHarel_clicked();

    void on_radioButton_GEM_clicked();

    void on_pushButton_degree_clicked();

    void on_pushButton_distance_clicked();

    void on_pushButton_between_clicked();

    void on_pushButton_save_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
