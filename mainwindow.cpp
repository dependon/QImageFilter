#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDebug>
#include "api.h"
#include <QDateTime>

#include <QDebug>
void skinImage(QImage *img, QImage *imgCopy)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            int R = qRed(line[x]);
            int G = qGreen(line[x]);
            int B = qBlue(line[x]);
            //            int Y=0.257 *qRed(line[x])+0.564 *qGreen(line[x])+0.098 *qBlue(line[x])+16;
            //            int Cb=-0.148 *qRed(line[x])-0.291 *qGreen(line[x])+0.439 *qBlue(line[x])+128;
            //            int Cr=0.439 *qRed(line[x])-0.368 *qGreen(line[x])-0.071 *qBlue(line[x])+128;

            int Y = ((R << 6) + (R << 1) + (G << 7) + (G << 4) + (B << 4) + (B << 3) + 3840) >> 8;
            int Cb = (-((R << 5) + (R << 2) + (R << 1)) - ((G << 6) + (G << 3) + (G << 2)) + ((B << 6) + (R << 5) + (R << 4)) + 32768) >> 8;
            int Cr = (((R << 6) + (R << 5) + (R << 4)) - ((G << 6) + (G << 4) + (G << 3) + (G << 2) + (G << 1)) - ((B << 4) + (B << 1)) + 32768) >> 8;

            //            qDebug()<<Cb<<Cb1;
            //            qDebug()<<Cr<<Cr1;
            if ((Cb > 77 && Cb < 127) && (Cr > 133 && Cr < 173)) {
                //imgCopy->setPixel(x,y, qRgb(0, 0, 0));
            } else
                imgCopy->setPixel(x, y, qRgb(255, 255, 255));
        }

    }
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle(tr("QImageFilter"));
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

void MainWindow::on_fushBtn_clicked()
{
    return;
    if (m_img) {
        for (int height = 0; height < m_img->height(); height++) {
            for (int width = 0; width < m_img->width(); width++) {


                QColor frontColor = m_img->pixel(width, height);
                QColor afterColor;
                float frontred = frontColor.red();
                float frontgreen = frontColor.green();
                float frontblue = frontColor.blue();
                if ((frontred / frontgreen) < 0.7) {
                    m_img->setPixel(width, height, qRgb(0, 0, 0));
                    qDebug() << 1111;
                } else {
//                    qDebug()<<frontred <<frontblue;
                }

            }
        }
        ui->labelcl->setPixmap(QPixmap::fromImage(*m_img));
        update();
    }
}

void MainWindow::on_pushButton_clicked()
{
    return;
    int index = 0;
    while (index++ < 1)
        if (m_img) {
            for (int height = 0; height < m_img->height(); height++) {
                for (int width = 0; width < m_img->width(); width++) {
                    QColor frontColor = m_img->pixel(width, height);
                    QColor afterColor;
                    float frontred = frontColor.red();
                    float frontgreen = frontColor.green();
                    float frontblue = frontColor.blue();
                    float afterred = 0.7 * frontred + 0.769 * frontgreen + 0.189 * frontblue;
                    float aftergreen = 0.7 * frontred + 0.686 * frontgreen + 0.168 * frontblue;
                    float afterblue = 0.5 * frontred + 0.534 * frontgreen + 0.131 * frontblue;

                    if (afterred > 255) {
                        afterred = 0;
                    }
                    if (aftergreen > 255) {
                        aftergreen = 0;
                    }
                    if (afterblue > 255) {
                        afterblue = 255;
                    }
                    m_img->setPixel(width, height, qRgb(afterred, aftergreen, afterblue));
                }
            }
            ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800, 600));
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
        skinImage(m_img, &m_imgCopy);
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
