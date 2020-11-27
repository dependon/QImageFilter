#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_cv1Btn_clicked();

    void on_openBtn_clicked();

    void on_fanseBtn_clicked();

    void on_fushBtn_clicked();

    void on_pushButton_clicked();

    void on_gugaiBtn_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;

    QImage *m_img{nullptr};
};

#endif // MAINWINDOW_H
