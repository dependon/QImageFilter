#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
class QIntValidator;
class QLineEdit;
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

    void on_gugaiBtn_clicked();

    void on_pushButton_2_clicked();

    void on_duibiSlider_valueChanged(int value);

    void on_lightSlider_valueChanged(int value);

    void on_GrayScaleBtn_clicked();

    void on_coolBtn_clicked();

    void on_warnBtn_clicked();

    void on_skinBtn_clicked();

    void on_Horizontal_clicked();

    void on_vertical_clicked();

    void on_lapace_clicked();

    void on_sobel_clicked();

    void on_stauration_clicked();

    void on_transparency_valueChanged(int value);

    void on_TransparentBtn_clicked();

    void on_saveBtn_clicked();

    void on_changeColorBtn_clicked();

    void on_pushButton_3_clicked();
    void on_GaussianBlurBtn_clicked();

    void on_maiskBtn_clicked();

    void on_byjcBtn_clicked();

    void on_lktqBtn_clicked();

    void on_hisBtn_clicked();

    void on_greayHisBtn_clicked();

    void on_averageBtn_clicked();

    void on_GaussianFilterBtn_clicked();

    void on_zzBtn_clicked();

    void on_BilateralBtn_clicked();

    void on_edgeDetectionBtn_clicked();

    void on_laplacianEdgeDetectionBtn_clicked();

private:
    Ui::MainWindow *ui;

    QImage *m_img{nullptr};
    QImage m_imgCopy{nullptr};
};

#endif // MAINWINDOW_H
