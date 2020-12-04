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

    void on_duibiSlider_valueChanged(int value);

    void on_lightSlider_valueChanged(int value);

    void on_GrayScaleBtn_clicked();

    void on_coolBtn_clicked();

    void on_warnBtn_clicked();

    void on_skinBtn_clicked();

private:
    Ui::MainWindow *ui;

    QImage *m_img{nullptr};
    QImage *m_imgCopy{nullptr};
};

#endif // MAINWINDOW_H
