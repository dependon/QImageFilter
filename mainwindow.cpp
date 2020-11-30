#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDebug>
#include "api.h"
#include <QDateTime>
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

void MainWindow::on_cv1Btn_clicked()
{
    if(m_img)
    {
        for(int height=0;height<m_img->height();height++)
        {
            for(int width=0;width<m_img->width();width++)
            {
                QColor frontColor=m_img->pixel(width,height);
                QColor afterColor;
                float frontred=frontColor.red();
                float frontgreen=frontColor.green();
                float frontblue=frontColor.blue();
                float afterred=0.393 *frontred+0.769 *frontgreen+0.189 *frontblue;
                float aftergreen=0.349 *frontred+0.686 *frontgreen+0.168 *frontblue;
                float afterblue=0.272 *frontred+0.534 *frontgreen+0.131 *frontblue;

                if(afterred>255)
                {
                    afterred=255;
                }
                if(aftergreen>255)
                {
                    aftergreen=255;
                }
                if(afterblue>255)
                {
                    afterblue=255;
                }
                m_img->setPixel(width,height,qRgb(afterred,aftergreen,afterblue));
            }
        }
        ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800,600));
        update();
    }
}

void MainWindow::on_openBtn_clicked()
{
    QString path=QFileDialog::getOpenFileName();
    if(m_img)
    {
        delete m_img;
        m_img=nullptr;
    }
    m_img=new QImage(path);

    ui->label->setPixmap(QPixmap::fromImage(*m_img).scaled(800,600));
    update();
}

void MainWindow::on_fanseBtn_clicked()
{

    if(m_img)
    {
        for(int height=0;height<m_img->height();height++)
        {
            for(int width=0;width<m_img->width();width++)
            {
                QColor frontColor=m_img->pixel(width,height);
                QColor afterColor;
                float afterred=255-frontColor.red();
                float aftergreen=255-frontColor.green();
                float afterblue=255-frontColor.blue();
                if(afterred>255)
                {
                    afterred=255;
                }
                if(aftergreen>255)
                {
                    aftergreen=255;
                }
                if(afterblue>255)
                {
                    afterblue=255;
                }
                m_img->setPixel(width,height,qRgb(afterred,aftergreen,afterblue));
            }
        }
        ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800,600));
        update();
    }
}

void MainWindow::on_fushBtn_clicked()
{
    if(m_img)
    {
        for(int height=0;height<m_img->height();height++)
        {
            for(int width=0;width<m_img->width();width++)
            {

//                QRgb pixel = pimg.pixel(k,i);

//                QColor rgb(pixel);                // 获取像素的rgb值

//                if(rgb.red()/(float)rgb.blue()<0.7)
//                        pimg.setPixel(k,i, QColor(0,0,0) );   // 黑色
//                else
//                        pimg.setPixel(k,i, rgb);
                QColor frontColor=m_img->pixel(width,height);
                QColor afterColor;
                float frontred=frontColor.red();
                float frontgreen=frontColor.green();
                float frontblue=frontColor.blue();
                if((frontred/frontgreen)<0.7)
                {
                    m_img->setPixel(width,height,qRgb(0,0,0));
                    qDebug()<<1111;
                }
                else {
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
    int index=0;
    while(index++ <1)
    if(m_img)
    {
        for(int height=0;height<m_img->height();height++)
        {
            for(int width=0;width<m_img->width();width++)
            {
                QColor frontColor=m_img->pixel(width,height);
                QColor afterColor;
                float frontred=frontColor.red();
                float frontgreen=frontColor.green();
                float frontblue=frontColor.blue();
//                float afterred=0.393 *frontred+0.769 *frontgreen+0.189 *frontblue;
//                float aftergreen=0.349 *frontred+0.686 *frontgreen+0.168 *frontblue;
//                float afterblue=0.272 *frontred+0.534 *frontgreen+0.131 *frontblue;
                float afterred=0.7 *frontred+0.769 *frontgreen+0.189 *frontblue;
                float aftergreen=0.7 *frontred+0.686 *frontgreen+0.168 *frontblue;
                float afterblue=0.5 *frontred+0.534 *frontgreen+0.131 *frontblue;
//                if(frontred>100 ||frontgreen>100 ||frontblue>100)
//                {
//                    afterred=255;
//                    aftergreen=255;
//                    afterblue=255;
//                }
//                else {
//                    afterred=0;
//                    aftergreen=0;
//                    afterblue=0;
//                }
                if(afterred>255)
                {
                    afterred=0;
                }
                if(aftergreen>255)
                {
                    aftergreen=0;
                }
                if(afterblue>255)
                {
                    afterblue=255;
                }
                m_img->setPixel(width,height,qRgb(afterred,aftergreen,afterblue));
            }
        }
        ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800,600));
        update();
    }
}
#include <QPainter>
void MainWindow::on_gugaiBtn_clicked()
{
    if(m_img)
    {
        QString path=QFileDialog::getOpenFileName();

        QImage m_imgB(path);

         QImage tmpFrame = m_imgB.scaled(QSize(m_img->width(), m_img->height()));
         QPainter painter;
         painter.begin(m_img);
         painter.drawImage(0, 0, tmpFrame);
         painter.end();
    }
    ui->labelcl->setPixmap(QPixmap::fromImage(*m_img).scaled(800,600));
    update();
}

void MainWindow::on_pushButton_2_clicked()
{
    if(m_img)
    {
        QImage *img=new QImage(*m_img);
        QImageD_RunBEEPSHorizontalVertical(m_img,img);
        qDebug()<<QDateTime::currentSecsSinceEpoch();
        ui->labelcl->setPixmap(QPixmap::fromImage(*img).scaled(800,600));
        update();
    }


}
