#ifndef API_H
#define API_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <QVector>
#include <QImage>

int Bound(int range_left,int data,int range_right)
{
    int index=data;
    if(data>range_right)
    {
        index=range_right;
    }
    else if (data<range_left) {
        index=range_left;
    }
    return index;
}

//垂直方向递归
void RunBEEPSVerticalHorizontal(double *data,int width,int height,double spatialDecay,double *exp_table,double *g_table)
{
    int length0=height*width;
    double* g= new double[length0];
    int m = 0;
    for (int k2 = 0;k2<height;++k2)
    {
        int n = k2;
        for (int k1 = 0;k1<width;++k1)
        {
            g[n]=data[m++];
            n += height;
        }
    }
    double*p = new double[length0];
    double*r = new double[length0];
    memcpy(p, g, sizeof(double) * length0);
    memcpy(r, g, sizeof(double) * length0);
    for (int k1 = 0;k1<width; ++k1)
    {
        int startIndex=k1 * height;
        double mu = 0.0;
        for (int k=startIndex+1,K =startIndex+height;k<K;++k)
        {
            int div0=fabs(p[k]-p[k-1]);
            mu =exp_table[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);//文献中的公式1，这里做了一下修改，效果影响不大
        }
        for (int k =startIndex+height-2;startIndex <= k;--k)
        {
            int div0=fabs(r[k]-r[k+1]);
            mu =exp_table[div0];
            r[k] = r[k+1] * mu + r[k] * (1.0-mu) ;//文献公式3
        }
    }
    double rho0=1.0/(2-spatialDecay);
    for (int k = 0;k <length0;++k)
    {
        r[k]= (r[k]+p[k])*rho0-g_table[(int)g[k]];
    }
    m = 0;
    for (int k1=0;k1<width;++k1)
    {
        int n = k1;
        for (int k2 =0;k2<height;++k2)
        {
            data[n] = r[m++];
            n += width;
        }
    }

    memcpy(p,data, sizeof(double) * length0);
    memcpy(r,data, sizeof(double) * length0);
    for (int k2 = 0; k2<height;++k2)
    {
        int startIndex=k2 * width;
        double mu = 0.0;
        for (int k=startIndex+1, K=startIndex+width;k<K;++k)
        {
            int div0=fabs(p[k]-p[k-1]);
            mu =exp_table[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);
        }
        for (int k=startIndex+width-2;startIndex<=k;--k)
        {
            int div0=fabs(r[k]-r[k+1]);
            mu =exp_table[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu) ;
        }
    }

    double init_gain_mu=spatialDecay/(2-spatialDecay);
    for (int k = 0; k <length0; k++)
    {
        data[k]=(p[k]+r[k])*rho0-data[k]*init_gain_mu;//文献中的公式5
    }
    delete[] p;
    delete[] r;
    delete[] g;
}

//水平方向递归
void RunBEEPSHorizontalVertical(double *data,int width,int height,double spatialDecay,double *exptable,double *g_table)
{
    int length=width*height;
    double* g = new double[length];
    double* p = new double[length];
    double* r = new double[length];
    memcpy(p,data, sizeof(double) * length);
    memcpy(r,data, sizeof(double) * length);
    double rho0=1.0/(2-spatialDecay);
    for (int k2 = 0;k2 < height;++k2)
    {
        int startIndex=k2 * width;
        for (int k=startIndex+1,K=startIndex+width;k<K;++k)
        {
            int div0=fabs(p[k]-p[k-1]);
            double mu =exptable[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);//文献公式1

        }

        for (int k =startIndex + width - 2;startIndex <= k;--k)
        {
            int div0=fabs(r[k]-r[k+1]);
            double mu =exptable[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu);//文献公式3
        }
        for (int k =startIndex,K=startIndex+width;k<K;k++)
        {
            r[k]=(r[k]+p[k])*rho0- g_table[(int)data[k]];
        }
    }

    int m = 0;
    for (int k2=0;k2<height;k2++)
    {
        int n = k2;
        for (int k1=0;k1<width;k1++)
        {
            g[n] = r[m++];
            n += height;
        }
    }

    memcpy(p, g, sizeof(double) * height * width);
    memcpy(r, g, sizeof(double) * height * width);
    for (int k1=0;k1<width;++k1)
    {
        int startIndex=k1 * height;
        double mu = 0.0;
        for (int k =startIndex+1,K =startIndex+height;k<K;++k)
        {
            int div0=fabs(p[k]-p[k-1]);
            mu =exptable[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);
        }
        for (int k=startIndex+height-2;startIndex<=k;--k)
        {
            int div0=fabs(r[k]-r[k+1]);
            mu =exptable[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu);
        }
    }

    double init_gain_mu=spatialDecay/(2-spatialDecay);
    for (int k = 0;k <length;++k)
    {
        r[k]= (r[k]+p[k])*rho0- g[k]*init_gain_mu;
    }
    m = 0;
    for (int k1=0;k1<width;++k1)
    {
        int n = k1;
        for (int k2=0;k2<height;++k2)
        {
            data[n]=r[m++];
            n += width;
        }
    }

    delete p;
    delete r;
    delete g;
}

//qimage磨皮
void QImageD_RunBEEPSHorizontalVertical(QImage *img,QImage *imgCopy,double spatialDecay=0.02,double photometricStandardDeviation=10)
{
    if(!img||!imgCopy){
        return ;
    }

    double c=-0.5/(photometricStandardDeviation * photometricStandardDeviation); //-1/2 *光度标准偏差的平方
    double mu=spatialDecay/(2-spatialDecay);

    double *exptable=new double[256];;
    double *g_table=new double[256];;
    for (int i=0;i<=255;i++)
    {
        exptable[i]=(1-spatialDecay)* exp(c*i*i);
        g_table[i]=mu*i;
    }
    int width=img->width();
    int height=img->height();
        int length=width*height;
    double* data2Red= new double[length];
    double* data2Green= new double[length];
    double* data2Blue= new double[length];

    int i=0;

    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            QRgb rgb=imgCopy->pixel(x,y);
            data2Red[i]=qRed(rgb);
            data2Green[i]=qGreen(rgb);
            data2Blue[i]=qBlue(rgb);
            i++;
        }
    }


    double* gRed = new double[length];
    double* pRed = new double[length];
    double* rRed = new double[length];

    double* gGreen = new double[length];
    double* pGreen = new double[length];
    double* rGreen = new double[length];

    double* gBlue = new double[length];
    double* pBlue = new double[length];
    double* rBlue = new double[length];
    memcpy(pRed,data2Red, sizeof(double) * length);
    memcpy(rRed,data2Red, sizeof(double) * length);

    memcpy(pGreen,data2Green, sizeof(double) * length);
    memcpy(rGreen,data2Green, sizeof(double) * length);

    memcpy(pBlue,data2Blue, sizeof(double) * length);
    memcpy(rBlue,data2Blue, sizeof(double) * length);


    double rho0=1.0/(2-spatialDecay);
    for (int k2 = 0;k2 < height;++k2)
    {
        int startIndex=k2 * width;
        double mu=0.0;
        for (int k=startIndex+1,K=startIndex+width;k<K;++k)
        {
            int div0Red=fabs(pRed[k]-pRed[k-1]);
            mu =exptable[div0Red];
            pRed[k] = pRed[k - 1] * mu + pRed[k] * (1.0 - mu);//公式1

            int div0Green=fabs(pGreen[k]-pGreen[k-1]);
            mu =exptable[div0Green];
            pGreen[k] = pGreen[k - 1] * mu + pGreen[k] * (1.0 - mu);//公式1

            int div0Blue=fabs(pBlue[k]-pBlue[k-1]);
            mu =exptable[div0Blue];
            pBlue[k] = pBlue[k - 1] * mu + pBlue[k] * (1.0 - mu);//公式1

        }

        for (int k =startIndex + width - 2;startIndex <= k;--k)
        {
            int div0Red=fabs(rRed[k]-rRed[k+1]);
            double mu =exptable[div0Red];
            rRed[k] = rRed[k + 1] * mu + rRed[k] * (1.0 - mu);//公式3

            int div0Green=fabs(rGreen[k]-rGreen[k+1]);
            mu =exptable[div0Green];
            rGreen[k] = rGreen[k + 1] * mu + rGreen[k] * (1.0 - mu);//公式3

            int div0Blue=fabs(rBlue[k]-rBlue[k+1]);
            mu =exptable[div0Blue];
            rBlue[k] = rBlue[k + 1] * mu + rBlue[k] * (1.0 - mu);//公式3
        }
        for (int k =startIndex,K=startIndex+width;k<K;k++)
        {
            rRed[k]=(rRed[k]+pRed[k])*rho0- g_table[(int)data2Red[k]];
            rGreen[k]=(rGreen[k]+pGreen[k])*rho0- g_table[(int)data2Green[k]];
            rBlue[k]=(rBlue[k]+pBlue[k])*rho0- g_table[(int)data2Blue[k]];
        }
    }

    int m = 0;
    for (int k2=0;k2<height;k2++)
    {
        int n = k2;
        for (int k1=0;k1<width;k1++)
        {
            gRed[n] = rRed[m];
            gGreen[n] = rGreen[m];
            gBlue[n] = rBlue[m];
            m++;
            n += height;
        }
    }

    memcpy(pRed, gRed, sizeof(double) * height * width);
    memcpy(rRed, gRed, sizeof(double) * height * width);

    memcpy(pGreen, gGreen, sizeof(double) * height * width);
    memcpy(rGreen, gGreen, sizeof(double) * height * width);

    memcpy(pBlue, gBlue, sizeof(double) * height * width);
    memcpy(rBlue, gBlue, sizeof(double) * height * width);

    for (int k1=0;k1<width;++k1)
    {
        int startIndex=k1 * height;
        double mu = 0.0;
        for (int k =startIndex+1,K =startIndex+height;k<K;++k)
        {
            int div0Red=fabs(pRed[k]-pRed[k-1]);
            mu =exptable[div0Red];
            pRed[k] = pRed[k - 1] * mu + pRed[k] * (1.0 - mu);

            int div0Green=fabs(pGreen[k]-pGreen[k-1]);
            mu =exptable[div0Green];
            pGreen[k] = pGreen[k - 1] * mu + pGreen[k] * (1.0 - mu);

            int div0Blue=fabs(pBlue[k]-pBlue[k-1]);
            mu =exptable[div0Blue];
            pBlue[k] = pBlue[k - 1] * mu + pBlue[k] * (1.0 - mu);
        }
        for (int k=startIndex+height-2;startIndex<=k;--k)
        {
            int div0Red=fabs(rRed[k]-rRed[k+1]);
            mu =exptable[div0Red];
            rRed[k] = rRed[k + 1] * mu + rRed[k] * (1.0 - mu);

            int div0Green=fabs(rGreen[k]-rGreen[k+1]);
            mu =exptable[div0Green];
            rGreen[k] = rGreen[k + 1] * mu + rGreen[k] * (1.0 - mu);

            int div0Blue=fabs(rBlue[k]-rBlue[k+1]);
            mu =exptable[div0Blue];
            rBlue[k] = rBlue[k + 1] * mu + rBlue[k] * (1.0 - mu);
        }
    }

    double init_gain_mu=spatialDecay/(2-spatialDecay);
    for (int k = 0;k <length;++k)
    {
        rRed[k]= (rRed[k]+pRed[k])*rho0- gRed[k]*init_gain_mu;

        rGreen[k]= (rGreen[k]+pGreen[k])*rho0- gGreen[k]*init_gain_mu;

        rBlue[k]= (rBlue[k]+pBlue[k])*rho0- gBlue[k]*init_gain_mu;

    }

    m = 0;
    for (int k1=0;k1<width;++k1)
    {
        int n = k1;
        for (int k2=0;k2<height;++k2)
        {

            data2Red[n]=rRed[m];
            data2Green[n]=rGreen[m];
            data2Blue[n]=rBlue[m];
            imgCopy->setPixel(k1,k2,qRgb(data2Red[n],data2Green[n],data2Blue[n]));
            m++;
            n += width;
        }
    }


    delete pRed;
    delete rRed;
    delete gRed;

    delete pGreen;
    delete rGreen;
    delete gGreen;

    delete pBlue;
    delete rBlue;
    delete gBlue;
}
void warnImage(QImage *img,QImage *imgCopy,int index=30)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    QColor frontColor;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){
            frontColor=line[x];
            float r=frontColor.red()+index;
            float g=frontColor.green()+index;
            float b=frontColor.blue();
            r=Bound(0,r,255);
            g=Bound(0,g,255);
            imgCopy->setPixel(x,y, qRgb(r, g, b));
        }

    }
}

void coolImage(QImage *img,QImage *imgCopy,int index=30)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    QColor frontColor;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){
            frontColor=line[x];
            float r=frontColor.red();
            float g=frontColor.green();
            float b=frontColor.blue()+index;
            b=Bound(0,b,255);
            imgCopy->setPixel(x,y, qRgb(r, g, b));
        }

    }
}

void GrayScaleImage(QImage *img,QImage *imgCopy)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){
            int average = (qRed(line[x]) + qGreen(line[x]) + qRed(line[x]))/3;
            imgCopy->setPixel(x,y, qRgb(average, average, average));
        }

    }
}

void lightContrastImage(QImage *img,QImage *imgCopy,int light=100,int Contrast=150)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){
            float r=light*0.01*qRed(line[x])-150+Contrast;
            float g=light*0.01*qGreen(line[x])-150+Contrast;
            float b=light*0.01*qBlue(line[x])-150+Contrast;
            r=Bound(0,r,255);
            g=Bound(0,g,255);
            b=Bound(0,b,255);
            imgCopy->setPixel(x,y, qRgb(r, g, b));
        }

    }
}

void InverseColorImage(QImage *img,QImage *imgCopy)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){

            imgCopy->setPixel(x,y, qRgb(255-qRed(line[x]), 255-qGreen(line[x]), 255-qBlue(line[x])));
        }

    }
}
void oldImage(QImage *img,QImage *imgCopy)
{
    if(!img||!imgCopy){
        return ;
    }

    QRgb * line;
    for(int y = 0; y<imgCopy->height(); y++){
        line = (QRgb *)imgCopy->scanLine(y);
        for(int x = 0; x<imgCopy->width(); x++){
            float r=0.393 *qRed(line[x])+0.769 *qGreen(line[x])+0.189 *qBlue(line[x]);
            float g=0.349 *qRed(line[x])+0.686 *qGreen(line[x])+0.168 *qBlue(line[x]);
            float b=0.272 *qRed(line[x])+0.534 *qGreen(line[x])+0.131 *qBlue(line[x]);
            r=Bound(0,r,255);
            g=Bound(0,g,255);
            b=Bound(0,b,255);
            imgCopy->setPixel(x,y, qRgb(r, g, b));
        }

    }
}
#endif // API_H
