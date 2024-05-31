/*
 * Copyright (C) 2020 ~ 2021 LiuMingHang.
 *
 * Author:     LiuMingHang <liuminghang0821@gmail.com>
 *
 * Maintainer: LiuMingHang <liuminghang0821@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef IIMAGEAPI_H
#define IIMAGEAPI_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <QVector>
#include <QImage>
#include <QColor>

class QImageAPI
{

public:

    QImageAPI();
    //Take RGB maximum
    static int RgbMax(int red, int green, int blue);
    //Take RGB minimum
    static int RgbMin(int red, int green, int blue);
    //Boundary judgment
    static int Bound(int range_left, int data, int range_right);
    //qimage Skin grinding
    static QImage QImageD_RunBEEPSHorizontalVertical(const QImage &img, double spatialDecay = 0.02, double photometricStandardDeviation = 10);
    //Warm color filter
    static QImage warnImage(const QImage &img, int index = 30);
    //Cool color filter
    static QImage coolImage(const QImage &img,  int index = 30);
    //Grayscale filter
    static QImage GrayScaleImage(const QImage &img);
    //Brightness and saturation
    static QImage lightContrastImage(const QImage &img, int light = 100, int Contrast = 150);
    //Anti color filter
    static QImage InverseColorImage(const QImage &img);
    //Old photo filter
    static QImage oldImage(const QImage &img);

    //laplacian sharpening
    static QImage LaplaceSharpen(const QImage &img);

    //Sobel Edge Detector
    static QImage SobelEdge(const QImage &img);

    //Greyscale
    static QImage GreyScale(const QImage &img);

    //Contour acquisition 轮廓提取
    static QImage ContourExtraction(const QImage &img);

    //Flip horizontally
    static QImage Horizontal(const QImage &img);

    //Flip vertical
    static QImage Vertical(const QImage &origin);

    //Binarization
    static QImage Binaryzation(const QImage &img);

    //Metal wire drawing effect
    static QImage Metal(const QImage &img);

    //Adjust image brightness
    static QImage Brightness(int delta, const QImage &img);

    //Transparency
    static QImage transparencyImg(int delta, const QImage &img);

    //Saturation (- 100 - 100)
    static QImage StaurationImg(const QImage &origin, int saturation);

    //指定颜色透明
    static QImage Transparent2Png(const QImage &bmp,QColor color);

    //指定颜色透明
    static QImage changeColor2Png(const QImage &bmp,QColor oldcolor,QColor newcolor);

    //磨皮算法
    static QImage smooth(QImage inputImage, float sigma);

    //高斯模糊(有点慢)
    static QImage applyGaussianBlur(const QImage &oldimage, int radius);

    //马赛克
    static QImage applyMosaic(const QImage &oldImage, int blockSize);

    //边缘检测 Prewitt
    static QImage applyPrewitt(const QImage &image);

    //边缘检测detect
    static QImage detectEdges(const QImage &inputImage);

    //皮肤识别
    static QImage skinImage(const QImage &img);

    //直方图均衡化
    static QImage equalizeHistogram(const QImage &inputImage);

    //得到图像灰度直方图
    static QImage drawHistogram(const QImage &image);

    //均值滤波
    static QImage applyMeanFilter(const QImage &inputImage);

    // 对图像进行高斯滤波
    static QImage applyGaussianFilter(const QImage &inputImage, double sigma);

    // 对图像进行中值滤波
    static QImage applyMedianFilter(const QImage &inputImage);

    //test
    static QImage applyBilateralFilter(const QImage &inputImage, double sigmaS, double sigmaR);
};


#endif // IIMAGEAPI_H
