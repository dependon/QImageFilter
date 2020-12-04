#ifndef API_H
#define API_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <QVector>
#include <QImage>
#include <QDebug>
#include <QThread>
#include "application.h"
#include <QThread>

class api:public QObject
{
    Q_OBJECT
public:
    static int Bound(int range_left,int data,int range_right);
    static void RGBToYUV(int Red, int Green, int Blue, int* Y,int* U,int* V);;
    static void YUVToRGB(int Y, int U, int V, int* Red, int* Green, int* Blue);;
    //垂直方向递归原始函数
    static void RunBEEPSVerticalHorizontal(double *data,int width,int height,double spatialDecay,double *exp_table,double *g_table);
    //水平方向递归原始函数
    static void RunBEEPSHorizontalVertical(double *data,int width,int height,double spatialDecay,double *exptable,double *g_table);
    //qimage磨皮
    static void QImageD_RunBEEPSHorizontalVertical(QImage *img,QImage *imgCopy,double spatialDecay=0.02,double photometricStandardDeviation=10);
    static void warnImage(QImage *img,QImage *imgCopy,int index=30);
    static void coolImage(QImage *img,QImage *imgCopy,int index=30);
    static void GrayScaleImage(QImage *img,QImage *imgCopy);
    static void lightContrastImage(QImage *img,QImage *imgCopy,int light=100,int Contrast=150);
    static void InverseColorImage(QImage *img,QImage *imgCopy);
    static void oldImage(QImage *img,QImage *imgCopy);

};

#endif // API_H
