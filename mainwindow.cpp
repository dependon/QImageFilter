#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDebug>
#include "api.h"
#include <QDateTime>
#include <QStandardPaths>
#include <QDebug>
#include <QRegExp>
#include <QRegExpValidator>
#include <QIntValidator>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle(tr("QImageFilter"));

    QRegularExpression regExp("^([0-9]|[1-9][0-9]|[1][0-9]{2}|2[0-4][0-9]|25[0-5])$");
    QValidator *validator = new QRegularExpressionValidator(regExp, this);

    ui->R1->setValidator(validator);
    ui->R2->setValidator(validator);
    ui->G1->setValidator(validator);
    ui->G2->setValidator(validator);
    ui->B1->setValidator(validator);
    ui->B2->setValidator(validator);
    ui->GaussianBlurEdit->setValidator(validator);

    QString path = "F:/lmh/github_code/QImageFilter/test/1.jpeg";
    if (m_img) {
        delete m_img;
        m_img = nullptr;
    }
    m_img = new QImage(path);

    ui->label->setPixmap(QPixmap::fromImage(*m_img).scaled(800, 600));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_cv1Btn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::oldImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_openBtn_clicked()
{
    QString path = QFileDialog::getOpenFileName();
    if (m_img) {
        delete m_img;
        m_img = nullptr;
    }
    m_img = new QImage(path);

    ui->label->setPixmap(QPixmap::fromImage(*m_img).scaled(800, 600));
    update();
}

void MainWindow::on_fanseBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::InverseColorImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

#include <QPainter>
void MainWindow::on_gugaiBtn_clicked()
{
    if (m_img) {
        QString path = QFileDialog::getOpenFileName();

        QImage m_imgB(path);

        QImage tmpFrame = m_imgB.scaled(QSize(m_img->width(), m_img->height()));
        QPainter painter;
        painter.begin(m_img);
        painter.drawImage(0, 0, tmpFrame);
        painter.end();
    }
    ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800, 600));
    update();
}

void MainWindow::on_pushButton_2_clicked()
{
    if (m_img) {
        QImage *img = new QImage(*m_img);
        m_imgCopy = QImageAPI::QImageD_RunBEEPSHorizontalVertical(*img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }


}

void MainWindow::on_duibiSlider_valueChanged(int value)
{
    if (m_img) {

        m_imgCopy = QImageAPI::lightContrastImage(*m_img, ui->duibiSlider->value(), ui->lightSlider->value());
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }

}

void MainWindow::on_lightSlider_valueChanged(int value)
{

    if (m_img) {

        m_imgCopy = QImageAPI::lightContrastImage(*m_img, ui->duibiSlider->value(), ui->lightSlider->value());
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_GrayScaleBtn_clicked()
{
    if (m_img) {

        m_imgCopy = QImageAPI::GrayScaleImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_coolBtn_clicked()
{
    if (m_img) {

        m_imgCopy = QImageAPI::coolImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_warnBtn_clicked()
{
    if (m_img) {

        m_imgCopy = QImageAPI::warnImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_skinBtn_clicked()
{
    if (m_img) {
        m_imgCopy = *m_img;
        m_imgCopy = QImageAPI::skinImage(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_Horizontal_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::Horizontal(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_vertical_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::Vertical(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_lapace_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::LaplaceSharpen(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_sobel_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::SobelEdge(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_stauration_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::StaurationImg(*m_img, 50);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_transparency_valueChanged(int value)
{
    if (m_img) {
        m_imgCopy = QImageAPI::transparencyImg(value, *m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_TransparentBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::Transparent2Png(*m_img, QColor(ui->R1->text().toInt(),ui->G1->text().toInt(),ui->B1->text().toInt()));
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_saveBtn_clicked()
{
//    QString strDir = QFileDialog::getExistingDirectory(
//                    this
//                    ,tr("Open Directory")
//                    ,"/home"
//                    ,QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    QString desktop = QStandardPaths::writableLocation(QStandardPaths::HomeLocation) + "/Desktop";
    QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), desktop, tr(".png")); //选择路径
    if (!filename.contains(".png")) {
        filename = filename + ".png";
    }
    m_imgCopy.save(filename);
}

void MainWindow::on_changeColorBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::changeColor2Png(*m_img, QColor(ui->R1->text().toInt(),ui->G1->text().toInt(),ui->B1->text().toInt())
                                               ,QColor(ui->R2->text().toInt(),ui->G2->text().toInt(),ui->B2->text().toInt()));
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_pushButton_3_clicked()
{
    if (m_img) {
//        m_imgCopy = QImageAPI::changeColor2Png(*m_img, QColor(220,83,39),QColor(249,188,61));
        m_imgCopy = QImageAPI::smooth(*m_img, 0.4);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}


void MainWindow::on_GaussianBlurBtn_clicked()
{
    int x = ui->GaussianBlurEdit->text().toInt();
    if (m_img) {
        m_imgCopy = QImageAPI::applyGaussianBlur(*m_img, x);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_maiskBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::applyMosaic(*m_img,ui->maiskEdit->text().toInt());
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_byjcBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::detectEdges(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_lktqBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::ContourExtraction(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_hisBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::equalizeHistogram(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_greayHisBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::drawHistogram(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_averageBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::applyMeanFilter(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_GaussianFilterBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::applyGaussianFilter(*m_img,ui->GaussianFilterEdit->text().toInt());
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_zzBtn_clicked()
{
    if (m_img) {
        m_imgCopy = QImageAPI::applyMedianFilter(*m_img);
        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}

void MainWindow::on_BilateralBtn_clicked()
{
    if (m_img) {

       // 应用双边滤波
        double sigmaS = ui->sigmaSEdit->text().toDouble(); // 空间域标准差
        double sigmaR = ui->sigmaREdit->text().toDouble(); // 灰度值域标准差
        m_imgCopy = QImageAPI::applyBilateralFilter(*m_img,sigmaS,sigmaR);

        ui->labelcl->setPixmap(QPixmap::fromImage(m_imgCopy).scaled(800, 600));
        update();
    }
}
