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
#include "api.h"

#include <QPainter>
#include <QDebug>
#include <QDateTime>
#include <QtMath>
#include <QRgb>
#include <QColor>

typedef struct  {
    uint8_t R;
    uint8_t G;
    uint8_t B;
    uint8_t L;
} RGBL;

typedef struct {
    float H;
    float S;
    float V;
} HSV;

// Prewitt算子
const int prewitt_x[3][3] = {
    {-1, 0, 1},
    {-1, 0, 1},
    {-1, 0, 1}
};

const int prewitt_y[3][3] = {
    {-1, -1, -1},
    {0, 0, 0},
    {1, 1, 1}
};

static void RGB_TO_HSV(const RGBL *input, HSV *output) // convert RGB value to HSV value
{
    float r, g, b, minRGB, maxRGB, deltaRGB;

    r = input->R / 255.0f;
    g = input->G / 255.0f;
    b = input->B / 255.0f;
    minRGB = QImageAPI::RgbMin(r, g, b);
    maxRGB = QImageAPI::RgbMax(r, g, b);
    deltaRGB = maxRGB - minRGB;

    output->V = maxRGB;
    if (maxRGB != 0.0f)
        output->S = deltaRGB / maxRGB;
    else
        output->S = 0.0f;
    if (output->S <= 0.0f) {
        output->H = 0.0f;
    } else {
        if (r == maxRGB) {
            output->H = (g - b) / deltaRGB;
        } else {
            if (g == maxRGB) {
                output->H = 2.0f + (b - r) / deltaRGB;
            } else {
                if (b == maxRGB) {
                    output->H = 4.0f + (r - g) / deltaRGB;
                }
            }
        }
        output->H = output->H * 60.0f;
        if (output->H < 0.0f) {
            output->H += 360;
        }
        output->H /= 360;
    }

}

static void HSV_TO_RGB(HSV *input, RGBL *output) //convert HSV value to RGB value
{
    float R, G, B;
    int k;
    float aa, bb, cc, f;
    if (input->S <= 0.0f)
        R = G = B = input->V;
    else {
        if (input->H == 1.0f)
            input->H = 0.0f;
        input->H *= 6.0f;
        k = (int)floor(input->H);
        f = input->H - k;
        aa = input->V * (1.0f - input->S);
        bb = input->V * (1.0f - input->S * f);
        cc = input->V * (1.0f - (input->S * (1.0f - f)));
        switch (k) {
        case 0:
            R = input->V;
            G = cc;
            B = aa;
            break;
        case 1:
            R = bb;
            G = input->V;
            B = aa;
            break;
        case 2:
            R = aa;
            G = input->V;
            B = cc;
            break;
        case 3:
            R = aa;
            G = bb;
            B = input->V;
            break;
        case 4:
            R = cc;
            G = aa;
            B = input->V;
            break;
        case 5:
            R = input->V;
            G = aa;
            B = bb;
            break;
        }
    }
    output->R = (unsigned char)(R * 255);
    output->G = (unsigned char)(G * 255);
    output->B = (unsigned char)(B * 255);
}

void adjustBrightness(RGBL &rgb_v, int step)
{
    HSV hsv_v;
    RGB_TO_HSV(&rgb_v, &hsv_v);
    rgb_v.L = hsv_v.V;
    rgb_v.L += step;
    if (rgb_v.L <= 0) {
        rgb_v.L = 1;
    } else if (rgb_v.L >= 100) {
        rgb_v.L = 100;
    }

    hsv_v.V = rgb_v.L / 100.0;
    HSV_TO_RGB(&hsv_v, &rgb_v);
}


QImageAPI::QImageAPI()
{

}

int QImageAPI::RgbMax(int red, int green, int blue)
{
    int ret = 0;
    if (red >= green) {
        ret = red;
    } else {
        ret = green;
    }
    if (ret < blue) {
        ret = blue;
    }
    return ret;
}

int QImageAPI::RgbMin(int red, int green, int blue)
{
    int ret = 0;
    if (red >= green) {
        ret = green;
    } else {
        ret = red;
    }
    if (ret > blue) {
        ret = blue;
    }
    return ret;
}

int QImageAPI::Bound(int range_left, int data, int range_right)
{
    int index = data;
    if (data > range_right) {
        index = range_right;
    } else if (data < range_left) {
        index = range_left;
    }
    return index;
}

QImage QImageAPI::QImageD_RunBEEPSHorizontalVertical(const QImage &img, double spatialDecay, double photometricStandardDeviation)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    QImage imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);

    double c = -0.5 / (photometricStandardDeviation * photometricStandardDeviation);
    double mu = spatialDecay / (2 - spatialDecay);

    double *exptable = new double[256];
    double *g_table = new double[256];
    for (int i = 0; i <= 255; i++) {
        exptable[i] = (1 - spatialDecay) * exp(c * i * i);
        g_table[i] = mu * i;
    }
    int width = img.width();
    int height = img.height();
    int length = width * height;
    double *data2Red = new double[length];
    double *data2Green = new double[length];
    double *data2Blue = new double[length];


    int size = imgCopy.width() * imgCopy.height();
    uint8_t *rgb = imgCopy.bits();
    for (int i = 0; i < size; i++) {
        data2Red[i] = rgb[i * 3];
        data2Green[i] = rgb[i * 3 + 1];
        data2Blue[i] = rgb[i * 3 + 2];
    }


    double *gRed = new double[length];
    double *pRed = new double[length];
    double *rRed = new double[length];

    double *gGreen = new double[length];
    double *pGreen = new double[length];
    double *rGreen = new double[length];

    double *gBlue = new double[length];
    double *pBlue = new double[length];
    double *rBlue = new double[length];
    memcpy(pRed, data2Red, sizeof(double) * length);
    memcpy(rRed, data2Red, sizeof(double) * length);

    memcpy(pGreen, data2Green, sizeof(double) * length);
    memcpy(rGreen, data2Green, sizeof(double) * length);

    memcpy(pBlue, data2Blue, sizeof(double) * length);
    memcpy(rBlue, data2Blue, sizeof(double) * length);


    double rho0 = 1.0 / (2 - spatialDecay);
    for (int k2 = 0; k2 < height; ++k2) {
        int startIndex = k2 * width;
        double mu = 0.0;
        for (int k = startIndex + 1, K = startIndex + width; k < K; ++k) {
            int div0Red = fabs(pRed[k] - pRed[k - 1]);
            mu = exptable[div0Red];
            pRed[k] = pRed[k - 1] * mu + pRed[k] * (1.0 - mu);

            int div0Green = fabs(pGreen[k] - pGreen[k - 1]);
            mu = exptable[div0Green];
            pGreen[k] = pGreen[k - 1] * mu + pGreen[k] * (1.0 - mu);

            int div0Blue = fabs(pBlue[k] - pBlue[k - 1]);
            mu = exptable[div0Blue];
            pBlue[k] = pBlue[k - 1] * mu + pBlue[k] * (1.0 - mu);

        }

        for (int k = startIndex + width - 2; startIndex <= k; --k) {
            int div0Red = fabs(rRed[k] - rRed[k + 1]);
            double mu = exptable[div0Red];
            rRed[k] = rRed[k + 1] * mu + rRed[k] * (1.0 - mu);

            int div0Green = fabs(rGreen[k] - rGreen[k + 1]);
            mu = exptable[div0Green];
            rGreen[k] = rGreen[k + 1] * mu + rGreen[k] * (1.0 - mu);

            int div0Blue = fabs(rBlue[k] - rBlue[k + 1]);
            mu = exptable[div0Blue];
            rBlue[k] = rBlue[k + 1] * mu + rBlue[k] * (1.0 - mu);
        }
        for (int k = startIndex, K = startIndex + width; k < K; k++) {
            rRed[k] = (rRed[k] + pRed[k]) * rho0 - g_table[(int)data2Red[k]];
            rGreen[k] = (rGreen[k] + pGreen[k]) * rho0 - g_table[(int)data2Green[k]];
            rBlue[k] = (rBlue[k] + pBlue[k]) * rho0 - g_table[(int)data2Blue[k]];
        }
    }

    int m = 0;
    for (int k2 = 0; k2 < height; k2++) {
        int n = k2;
        for (int k1 = 0; k1 < width; k1++) {
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

    for (int k1 = 0; k1 < width; ++k1) {
        int startIndex = k1 * height;
        double mu = 0.0;
        for (int k = startIndex + 1, K = startIndex + height; k < K; ++k) {
            int div0Red = fabs(pRed[k] - pRed[k - 1]);
            mu = exptable[div0Red];
            pRed[k] = pRed[k - 1] * mu + pRed[k] * (1.0 - mu);

            int div0Green = fabs(pGreen[k] - pGreen[k - 1]);
            mu = exptable[div0Green];
            pGreen[k] = pGreen[k - 1] * mu + pGreen[k] * (1.0 - mu);

            int div0Blue = fabs(pBlue[k] - pBlue[k - 1]);
            mu = exptable[div0Blue];
            pBlue[k] = pBlue[k - 1] * mu + pBlue[k] * (1.0 - mu);
        }
        for (int k = startIndex + height - 2; startIndex <= k; --k) {
            int div0Red = fabs(rRed[k] - rRed[k + 1]);
            mu = exptable[div0Red];
            rRed[k] = rRed[k + 1] * mu + rRed[k] * (1.0 - mu);

            int div0Green = fabs(rGreen[k] - rGreen[k + 1]);
            mu = exptable[div0Green];
            rGreen[k] = rGreen[k + 1] * mu + rGreen[k] * (1.0 - mu);

            int div0Blue = fabs(rBlue[k] - rBlue[k + 1]);
            mu = exptable[div0Blue];
            rBlue[k] = rBlue[k + 1] * mu + rBlue[k] * (1.0 - mu);
        }
    }

    double init_gain_mu = spatialDecay / (2 - spatialDecay);
    for (int k = 0; k < length; ++k) {
        rRed[k] = (rRed[k] + pRed[k]) * rho0 - gRed[k] * init_gain_mu;

        rGreen[k] = (rGreen[k] + pGreen[k]) * rho0 - gGreen[k] * init_gain_mu;

        rBlue[k] = (rBlue[k] + pBlue[k]) * rho0 - gBlue[k] * init_gain_mu;

    }

    m = 0;
    int nRowBytes = (width * 24 + 31) / 32 * 4;
    int  lineNum_24 = 0;
    for (int k1 = 0; k1 < width; ++k1) {
        int n = k1;
        for (int k2 = 0; k2 < height; ++k2) {

            data2Red[n] = rRed[m];
            data2Green[n] = rGreen[m];
            data2Blue[n] = rBlue[m];
            lineNum_24 = k2 * nRowBytes;
            rgb[lineNum_24 + k1 * 3] = data2Red[n];
            rgb[lineNum_24 + k1 * 3 + 1] = data2Green[n];
            rgb[lineNum_24 + k1 * 3 + 2] = data2Blue[n];
            m++;
            n += width;
        }
    }
    delete []data2Red;
    data2Red = nullptr;
    delete []data2Green ;
    data2Green = nullptr;
    delete []data2Blue;
    data2Blue = nullptr;

    delete []pRed;
    pRed = nullptr;
    delete []rRed;
    rRed = nullptr;
    delete []gRed;
    gRed = nullptr;

    delete []pGreen;
    pGreen = nullptr;
    delete []rGreen;
    rGreen = nullptr;
    delete []gGreen;
    gGreen = nullptr;

    delete []pBlue;
    pBlue = nullptr;
    delete []rBlue;
    rBlue = nullptr;
    delete []gBlue;
    gBlue = nullptr;

    delete []exptable;
    exptable = nullptr;
    delete []g_table;
    g_table = nullptr;

    qDebug() << "磨皮结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::warnImage(const QImage &img, int index)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    }
    QColor frontColor;
    int size = img.width() * img.height();

    for (int i = 0; i < size ; i++) {
        int r = rgb[i * 3] + index;
        int g = rgb[i * 3 + 1] + index;
        int b = rgb[i * 3 + 2] ;

        rgb[i * 3] = r > 255 ? 255 : r;
        rgb[i * 3 + 1] = g > 255 ? 255 : g;
        rgb[i * 3 + 2] = b > 255 ? 255 : b;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::coolImage(const QImage &img,  int index)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    }
    QColor frontColor;
    int size = img.width() * img.height();

    for (int i = 0; i < size ; i++) {
        int r = rgb[i * 3] ;
        int g = rgb[i * 3 + 1] ;
        int b = rgb[i * 3 + 2] + index;

        rgb[i * 3] = r > 255 ? 255 : r;
        rgb[i * 3 + 1] = g > 255 ? 255 : g;
        rgb[i * 3 + 2] = b > 255 ? 255 : b;
    }

    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::GrayScaleImage(const QImage &img)
{
#if 1
    QImage grayImage = img.convertToFormat(QImage::Format_Grayscale8);
    return grayImage;
#else
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    }
    QColor frontColor;
    int size = img.width() * img.height();

    for (int i = 0; i < size ; i++) {
        int average = (rgb[i * 3] + rgb[i * 3 + 1] + rgb[i * 3 + 2]) / 3;
        rgb[i * 3] = average > 255 ? 255 : average;
        rgb[i * 3 + 1] = average > 255 ? 255 : average;
        rgb[i * 3 + 2] = average > 255 ? 255 : average;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
#endif
}

QImage QImageAPI::lightContrastImage(const QImage &img,  int light, int Contrast)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    }
    int r;
    int g;
    int b;
    int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        r = light * 0.01 * rgb[i * 3] - 150 + Contrast;
        g = light * 0.01 * rgb[i * 3 + 1] - 150 + Contrast;
        b = light * 0.01 * rgb[i * 3 + 2]  - 150 + Contrast;
        r = Bound(0, r, 255);
        g = Bound(0, g, 255);
        b = Bound(0, b, 255);
        rgb[i * 3] = r;
        rgb[i * 3 + 1] = g;
        rgb[i * 3 + 2] = b;
    }

    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::InverseColorImage(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    } int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        rgb[i * 3] = 255 - rgb[i * 3] ;
        rgb[i * 3 + 1] = 255 - rgb[i * 3 + 1]  ;
        rgb[i * 3 + 2] = 255 - rgb[i * 3 + 2]  ;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::oldImage(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    if (nullptr == rgb) {
        return QImage();
    } int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        float r = 0.393 * rgb[i * 3] + 0.769 * rgb[i * 3 + 1] + 0.189 * rgb[i * 3 + 2];
        float g = 0.349 * rgb[i * 3] + 0.686 * rgb[i * 3 + 1] + 0.168 * rgb[i * 3 + 2];
        float b = 0.272 * rgb[i * 3] + 0.534 * rgb[i * 3 + 1] + 0.131 * rgb[i * 3 + 2];
        r = Bound(0, r, 255);
        g = Bound(0, g, 255);
        b = Bound(0, b, 255);
        rgb[i * 3] = r;
        rgb[i * 3 + 1] = g ;
        rgb[i * 3 + 2] = b  ;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}


QImage QImageAPI::LaplaceSharpen(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;
    int width = img.width();
    int height = img.height();
    int window[3][3] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(width, height, QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    QImage imgCopyrgbImg = QImage(img).convertToFormat(QImage::Format_RGB888);
    uint8_t *rgbImg = imgCopyrgbImg.bits();
    uint8_t *rgb = imgCopy.bits();

    int nRowBytes = (width * 24 + 31) / 32 * 4;
    int  lineNum_24 = 0;
    for (int x = 1; x < img.width(); x++) {
        for (int y = 1; y < img.height(); y++) {
            int sumR = 0;
            int sumG = 0;
            int sumB = 0;



            for (int m = x - 1; m <= x + 1; m++)
                for (int n = y - 1; n <= y + 1; n++) {
                    if (m >= 0 && m < width && n < height) {
                        lineNum_24 = n * nRowBytes;
                        sumR += rgbImg[lineNum_24 + m * 3] * window[n - y + 1][m - x + 1];
                        sumG += rgbImg[lineNum_24 + m * 3 + 1] * window[n - y + 1][m - x + 1];
                        sumB += rgbImg[lineNum_24 + m * 3 + 2] * window[n - y + 1][m - x + 1];
                    }
                }


            int old_r = rgbImg[lineNum_24 + x * 3];
            sumR += old_r;
            sumR = qBound(0, sumR, 255);

            int old_g = rgbImg[lineNum_24 + x * 3 + 1];
            sumG += old_g;
            sumG = qBound(0, sumG, 255);

            int old_b = rgbImg[lineNum_24 + x * 3 + 2];
            sumB += old_b;
            sumB = qBound(0, sumB, 255);
            lineNum_24 = y * nRowBytes;
            rgb[lineNum_24 + x * 3] = sumR;
            rgb[lineNum_24 + x * 3 + 1] = sumG;
            rgb[lineNum_24 + x * 3 + 2] = sumB;
        }
    }

    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::SobelEdge(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    double *Gx = new double[9];
    double *Gy = new double[9];

    /* Sobel */
    Gx[0] = 1.0; Gx[1] = 0.0; Gx[2] = -1.0;
    Gx[3] = 2.0; Gx[4] = 0.0; Gx[5] = -2.0;
    Gx[6] = 1.0; Gx[7] = 0.0; Gx[8] = -1.0;

    Gy[0] = -1.0; Gy[1] = -2.0; Gy[2] = - 1.0;
    Gy[3] = 0.0; Gy[4] = 0.0; Gy[5] = 0.0;
    Gy[6] = 1.0; Gy[7] = 2.0; Gy[8] = 1.0;

    QRgb pixel;
    QImage grayImage = GreyScale(img);
    int height = grayImage.height();
    int width = grayImage.width();
    QImage imgCopy = QImage(width, height, QImage::Format_RGB888);

    uint8_t *rgbImg = grayImage.bits();
    uint8_t *rgb = imgCopy.bits();

    int nRowBytes = (width * 24 + 31) / 32 * 4;
    int  lineNum_24 = 0;

    float *sobel_norm = new float[width * height];
    float max = 0.0;
    QColor my_color;

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            double value_gx = 0.0;
            double value_gy = 0.0;

            for (int k = 0; k < 3; k++) {
                for (int p = 0; p < 3; p++) {
                    if ((x + 1 + 1 - k < width) && (y + 1 + 1 - p < height)) {
                        pixel = grayImage.pixel(x + 1 + 1 - k, y + 1 + 1 - p);
                        lineNum_24 = (y + 1 + 1 - p) * nRowBytes;
                        value_gx += Gx[p * 3 + k] * rgbImg[lineNum_24 + (x + 1 + 1 - k) * 3];
                        value_gy += Gy[p * 3 + k] * rgbImg[lineNum_24 + (x + 1 + 1 - k) * 3];
                    }
                }
                sobel_norm[x + y * width] = abs(value_gx) + abs(value_gy);

                max = sobel_norm[x + y * width] > max ? sobel_norm[x + y * width] : max;
            }
        }
    }

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            my_color.setHsv(0, 0, 255 - int(255.0 * sobel_norm[i + j * width] / max));

            lineNum_24 = j * nRowBytes;
            rgb[lineNum_24 + i * 3] = my_color.red();
            rgb[lineNum_24 + i * 3 + 1] = my_color.green();
            rgb[lineNum_24 + i * 3 + 2] = my_color.blue();
        }
    }
    delete[] sobel_norm;

    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}


QImage QImageAPI::GreyScale(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;

    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        int average = (rgb[i * 3] * 299 + rgb[i * 3 + 1] * 587 + rgb[i * 3 + 1] * 114 + 500) / 1000;
        rgb[i * 3] = average;
        rgb[i * 3 + 1] = average;
        rgb[i * 3 + 2] = average;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;

}


QImage QImageAPI::Binaryzation(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();

    QImage imgCopy;

    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    int newGray = 0;
    int gray = 0;
    int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        gray = (rgb[i * 3] + rgb[i * 3 + 1] + rgb[i * 3 + 2]) / 3;
        if (gray > 128)
            newGray = 255;
        else
            newGray = 0;
        rgb[i * 3] = newGray;
        rgb[i * 3 + 1] = newGray;
        rgb[i * 3 + 2] = newGray;
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;

    return imgCopy;
}


QImage QImageAPI::ContourExtraction(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    int width = img.width();
    int height = img.height();
    int pixel[8];
    QImage binImg = Binaryzation(img);
    QImage newImg = QImage(width, height, QImage::Format_RGB888);
    newImg.fill(Qt::white);

    uint8_t *rgb = newImg.bits();
    uint8_t *binrgb = binImg.bits();
    int nRowBytes = (width * 24 + 31) / 32 * 4;
    int  lineNum_24 = 0;
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            memset(pixel, 0, 8);
            lineNum_24 = y * nRowBytes;
            if (binrgb[lineNum_24 + x * 3] == 0) {
                rgb[lineNum_24 + x * 3] = 0;
                rgb[lineNum_24 + x * 3 + 1] = 0;
                rgb[lineNum_24 + x * 3 + 2] = 0;
                pixel[0] = binrgb[(y - 1) * nRowBytes + (x - 1) * 3];
                pixel[1] = binrgb[(y) * nRowBytes + (x - 1) * 3];
                pixel[2] = binrgb[(y + 1) * nRowBytes + (x - 1) * 3];
                pixel[3] = binrgb[(y - 1) * nRowBytes + (x) * 3];
                pixel[4] = binrgb[(y + 1) * nRowBytes + (x) * 3];
                pixel[5] = binrgb[(y - 1) * nRowBytes + (x + 1) * 3];
                pixel[6] = binrgb[(y) * nRowBytes + (x + 1) * 3];
                pixel[7] = binrgb[(y + 1) * nRowBytes + (x + 1) * 3];

                if (pixel[0] + pixel[1] + pixel[2] + pixel[3] + pixel[4] + pixel[5] + pixel[6] + pixel[7] == 0) {
                    rgb[lineNum_24 + x * 3] = 255;
                    rgb[lineNum_24 + x * 3 + 1] = 255;
                    rgb[lineNum_24 + x * 3 + 2] = 255;
                }

            }
        }
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return newImg;
}


/*****************************************************************************
 *                                   Flip
 * **************************************************************************/
QImage QImageAPI::Horizontal(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    QImage copyImage(QSize(img.width(), img.height()), QImage::Format_ARGB32);
    copyImage = img.mirrored(true, false);
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return copyImage;

}


QImage QImageAPI::Metal(const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    QImage *baseImage = new QImage(img);
    QImage darkImage = QImageAPI::Brightness(-100, img);
    QImage greyImage = QImageAPI::GreyScale(darkImage);
    QPainter painter;

    QImage newImage = baseImage->scaled(QSize(img.width(), img.height()));

    painter.begin(&newImage);
    painter.setOpacity(0.5);
    painter.drawImage(0, 0, greyImage);
    painter.end();
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return newImage;
}


QImage QImageAPI::Brightness(int delta, const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();


    int r, g, b;
    QImage imgCopy;

    if (img.format() != QImage::Format_RGB888) {
        imgCopy = QImage(img).convertToFormat(QImage::Format_RGB888);
    } else {
        imgCopy = QImage(img);
    }
    uint8_t *rgb = imgCopy.bits();
    int size = img.width() * img.height();
    for (int i = 0; i < size ; i++) {
        r = rgb[i * 3] + delta;
        g = rgb[i * 3 + 1] + delta;
        b = rgb[i * 3 + 2] + delta;
        r = qBound(0, r, 255);
        g = qBound(0, g, 255);
        b = qBound(0, b, 255);
        rgb[i * 3] = r ;
        rgb[i * 3 + 1] =  g ;
        rgb[i * 3 + 2] =  b ;
    }

    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return imgCopy;
}

QImage QImageAPI::transparencyImg(int delta, const QImage &img)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    QImage newImage(img.width(), img.height(),
                    QImage::Format_ARGB32);
    QColor oldColor;
    int r, g, b;
    for (int x = 0; x < newImage.width(); x++) {
        for (int y = 0; y < newImage.height(); y++) {
            oldColor = QColor(img.pixel(x, y));

            r = oldColor.red() ;
            g = oldColor.green() ;
            b = oldColor.blue() ;

            newImage.setPixel(x, y, qRgba(r, g, b, delta));
        }
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return newImage;
}

QImage QImageAPI::StaurationImg(const QImage &origin, int saturation)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    int r, g, b, rgbMin, rgbMax;
    float k = saturation / 100.0f * 128;
    int alpha = 0;

    QImage newImage(origin);
    QColor tmpColor;

    for (int x = 0; x < newImage.width(); x++) {
        for (int y = 0; y < newImage.height(); y++) {
            tmpColor = QColor(origin.pixel(x, y));
            r = tmpColor.red();
            g = tmpColor.green();
            b = tmpColor.blue();

            rgbMin = RgbMin(r, g, b);
            rgbMax = RgbMax(r, g, b);

            int delta = (rgbMax - rgbMin);
            int value = (rgbMax + rgbMin);
            if (delta == 0) {
                continue;
            }
            int L = value >> 1;
            int S = L < 128 ? (delta << 7) / value : (delta << 7) / (510 - value);
            if (k >= 0) {
                alpha = k + S >= 128 ? S : 128 - k;
                alpha = 128 * 128 / alpha - 128;
            } else
                alpha = k;
            r = r + ((r - L) * alpha >> 7);
            g = g + ((g - L) * alpha >> 7);
            b = b + ((b - L) * alpha >> 7);
            r = Bound(0, r, 255);
            g = Bound(0, g, 255);
            b = Bound(0, b, 255);
            newImage.setPixel(x, y, qRgb(r, g, b));

        }
    }
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return newImage;

}


QImage QImageAPI::Vertical(const QImage &origin)
{
    qint64 startTime = QDateTime::currentMSecsSinceEpoch();
    QImage newImage(QSize(origin.width(), origin.height()), QImage::Format_ARGB32);
    newImage = origin.mirrored(false, true);
    qDebug() << "结束:" << QDateTime::currentMSecsSinceEpoch() - startTime;
    return newImage;
}

QImage QImageAPI::Transparent2Png(const QImage &bmp,QColor color)
{
    //BMP颜色格式转换成RGBA颜色格式
    QImage rebmp  = bmp.convertToFormat(QImage::Format_RGBA8888_Premultiplied,Qt::NoFormatConversion);
    int bmpWidth = bmp.width();
    int bmpHeight = bmp.height();
    //透明颜色
    QColor bmpBack= color;
    QColor bmpBackA(254,254,254,0);
    for(int i=0;i< bmpWidth;++i)
    {
        for(int j=0;j<bmpHeight;++j)
        {
            //如果身份证背景色等于 color,则设置为透明色 Color(254,254,254,0)
            if( bmp.pixelColor(i,j)== bmpBack)
            {
                rebmp.setPixelColor(i,j,bmpBackA);
            }
        }
    }
    return rebmp;
}

QImage QImageAPI::changeColor2Png(const QImage &bmp, QColor oldcolor, QColor newcolor)
{
    //BMP颜色格式转换成RGBA颜色格式
    QImage rebmp  = bmp.convertToFormat(QImage::Format_RGBA8888_Premultiplied,Qt::NoFormatConversion);
    int bmpWidth = bmp.width();
    int bmpHeight = bmp.height();
    //透明颜色
    QColor bmpBack= oldcolor;
    QColor bmpBackA(newcolor);
    for(int i=0;i< bmpWidth;++i)
    {
        for(int j=0;j<bmpHeight;++j)
        {
            if(bmp.pixelColor(i,j)==bmpBack)
            {
                rebmp.setPixelColor(i,j,bmpBackA);
            }
        }
    }
    return rebmp;
}

// 声明一个函数，用于对输入的图像进行磨皮处理
QImage QImageAPI::smooth(QImage inputImage, float sigma)
{
    int ksize = sigma * 5;
    if (ksize % 2 == 0) ksize += 1; // 去除偶数卷积核
    QImage outputImage = inputImage.copy();

    // 处理过程中使用的活动窗口滤波器
    QVector<float> kernel;
    const float PI = 3.141592653589793f;
    for (int i = -ksize / 2; i <= ksize / 2; i++)
    {
        kernel.append(1.0f / (2.0f * PI * sigma * sigma) * exp(-(i * i) / (2.0f * sigma * sigma)));
    }
    float sum_kernel = std::accumulate(kernel.constBegin(), kernel.constEnd(), 0.0f);

    int width = inputImage.width();
    int height = inputImage.height();
    int bytesPerPixel = inputImage.depth() / 8;

    for (int y = 0; y < height; ++y)
    {
        // 获取每一行的首地址
        uchar *inputScanLine = inputImage.scanLine(y);
        uchar *outputScanLine = outputImage.scanLine(y);

        // 处理每个像素点的 RGB 值
        for (int x = 0; x < width; ++x)
        {
            QVector<float> r, g, b;
            float sum_r = 0, sum_g = 0, sum_b = 0;
            // 遍历窗口内的像素点
            for (int i = -ksize / 2; i <= ksize / 2; i++)
            {
                int px = qBound(0, x + i, width - 1);
                uchar *pixel = inputScanLine + px * bytesPerPixel;
                r.append(pixel[2]);
                g.append(pixel[1]);
                b.append(pixel[0]);
            }
            // 按照像素点所在位置计算权重和，对 RGB 值进行同步处理
            for (int i = 0; i < r.size(); i++)
            {
                sum_r += kernel[i] * r[i];
                sum_g += kernel[i] * g[i];
                sum_b += kernel[i] * b[i];
            }
            sum_r = qBound(0.0f, sum_r / sum_kernel, 255.0f);
            sum_g = qBound(0.0f, sum_g / sum_kernel, 255.0f);
            sum_b = qBound(0.0f, sum_b / sum_kernel, 255.0f);

            // 将处理好的像素点写入图像中
            uchar *outputPixel = outputScanLine + x * bytesPerPixel;
            outputPixel[2] = static_cast<uchar>(sum_r);
            outputPixel[1] = static_cast<uchar>(sum_g);
            outputPixel[0] = static_cast<uchar>(sum_b);
        }
    }
    return outputImage;
}

// 高斯模糊函数
QImage QImageAPI::applyGaussianBlur(const QImage &oldimage, int radius)
{
    QImage image(oldimage);
    if (image.isNull() || radius <= 0)
        return QImage();

    QImage resultImage = image;
    const int size = radius * 2 + 1;
    const int sigma = radius / 2;
    const double sigmaSq = sigma * sigma;
    QVector<double> kernel(size);

    // 构建高斯核
    double sum = 0.0;
    for (int i = -radius; i <= radius; ++i)
    {
        double value = exp(-(i * i) / (2 * sigmaSq)) / (sqrt(2 * M_PI) * sigma);
        kernel[i + radius] = value;
        sum += value;
    }

    // 归一化
    for (int i = 0; i < size; ++i)
    {
        kernel[i] /= sum;
    }

    // 水平方向模糊
    for (int y = 0; y < image.height(); ++y)
    {
        for (int x = radius; x < image.width() - radius; ++x)
        {
            double red = 0, green = 0, blue = 0;
            for (int i = -radius; i <= radius; ++i)
            {
                QRgb pixel = image.pixel(x + i, y);
                red += qRed(pixel) * kernel[i + radius];
                green += qGreen(pixel) * kernel[i + radius];
                blue += qBlue(pixel) * kernel[i + radius];
            }
            resultImage.setPixel(x, y, qRgb(red, green, blue));
        }
    }

    // 垂直方向模糊
    for (int x = 0; x < image.width(); ++x)
    {
        for (int y = radius; y < image.height() - radius; ++y)
        {
            double red = 0, green = 0, blue = 0;
            for (int i = -radius; i <= radius; ++i)
            {
                QRgb pixel = resultImage.pixel(x, y + i);
                red += qRed(pixel) * kernel[i + radius];
                green += qGreen(pixel) * kernel[i + radius];
                blue += qBlue(pixel) * kernel[i + radius];
            }
            image.setPixel(x, y, qRgb(red, green, blue));
        }
    }
    return image;
}


QImage QImageAPI::applyMosaic(const QImage& oldImage, int blockSize) {
    if (oldImage.isNull() || blockSize <= 0) {
        return QImage(); // 返回空图片或处理错误
    }

    // 确保blockSize是偶数，并且不会使图像尺寸变得太小
    blockSize = (blockSize % 2 == 0) ? blockSize : blockSize + 1;
    if (oldImage.width() < blockSize || oldImage.height() < blockSize) {
        return oldImage; // 如果blockSize太大，直接返回原图
    }

    // 计算新图片的尺寸
    int newWidth = oldImage.width() / blockSize * blockSize;
    int newHeight = oldImage.height() / blockSize * blockSize;

    QImage newImage(newWidth, newHeight, oldImage.format());

    // 遍历每个块
    for (int y = 0; y < newHeight; y += blockSize) {
        for (int x = 0; x < newWidth; x += blockSize) {
            // 计算块的平均颜色
            QRgb averageColor = qRgb(0, 0, 0); // 初始化平均颜色为黑色
            int totalR = 0, totalG = 0, totalB = 0;
            int count = 0;

            for (int by = 0; by < blockSize && y + by < oldImage.height(); ++by) {
                for (int bx = 0; bx < blockSize && x + bx < oldImage.width(); ++bx) {
                    QRgb pixel = oldImage.pixel(x + bx, y + by);
                    totalR += qRed(pixel);
                    totalG += qGreen(pixel);
                    totalB += qBlue(pixel);
                    ++count;
                }
            }

            if (count > 0) { // 确保count不是0，避免除以0
                averageColor = qRgb(totalR / count, totalG / count, totalB / count);
            }

            // 用平均颜色填充整个块
            for (int by = 0; by < blockSize && y + by < newImage.height(); ++by) {
                for (int bx = 0; bx < blockSize && x + bx < newImage.width(); ++bx) {
                    newImage.setPixel(x + bx, y + by, averageColor);
                }
            }
        }
    }

    return newImage;
}

QImage QImageAPI:: applyPrewitt(const QImage &image)
{
    QImage result(image.size(), QImage::Format_Grayscale8);
    result.fill(Qt::black);

    // 确保图像是灰度图像
    QImage grayImage;
    if (image.format() != QImage::Format_Grayscale8) {
        grayImage = image.convertToFormat(QImage::Format_Grayscale8);
    } else {
        grayImage = image;
    }

    int width = grayImage.width();
    int height = grayImage.height();

    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int gx = 0, gy = 0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    QRgb pixel = grayImage.pixel(x + j, y + i);
                    int gray = qGray(pixel); // 从QRgb转换为灰度值
                    gx += gray * prewitt_x[i + 1][j + 1];
                    gy += gray * prewitt_y[i + 1][j + 1];
                }
            }

            // 计算梯度幅度（这里使用简单的平方和的平方根作为近似，但为简化可以使用绝对值之和）
            int gradientMagnitude = qAbs(gx) + qAbs(gy);

            // 设置阈值以决定边缘的亮度
            int threshold = 5; // 你可以调整这个阈值
            if (gradientMagnitude > threshold) {
                result.setPixel(x, y, 255); // 使用白色表示边缘
            }
        }
    }

    return result;
}

// Sobel算子
int sobelOperator(const QImage &image, int x, int y)
{
    int gx = 0, gy = 0;

    // Sobel算子
    int sobelX[3][3] = {{-1, 0, 1},
                        {-2, 0, 2},
                        {-1, 0, 1}};
    int sobelY[3][3] = {{-1, -2, -1},
                        {0, 0, 0},
                        {1, 2, 1}};

    // 遍历Sobel算子的3x3邻域
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            // 获取邻域内的像素值，超出边界的像素使用0代替
            int pixelX = qBound(0, x + i, image.width() - 1);
            int pixelY = qBound(0, y + j, image.height() - 1);
            QColor pixelColor(image.pixel(pixelX, pixelY));

            // 计算梯度值
            gx += sobelX[i + 1][j + 1] * pixelColor.red();
            gy += sobelY[i + 1][j + 1] * pixelColor.red();
        }
    }

    // 计算梯度的幅值
    int gradientMagnitude = qAbs(gx) + qAbs(gy);

    // 对梯度值进行归一化处理，确保在[0, 255]范围内
    gradientMagnitude = qBound(0, gradientMagnitude, 255);

    return gradientMagnitude;
}

// 边缘检测函数
QImage QImageAPI::detectEdges(const QImage &inputImage)
{
    QImage outputImage(inputImage.size(), inputImage.format());

    for (int y = 0; y < inputImage.height(); ++y) {
        for (int x = 0; x < inputImage.width(); ++x) {
            // 对每个像素应用Sobel算子
            int gradientMagnitude = sobelOperator(inputImage, x, y);
            // 将梯度值作为边缘强度，用灰度值表示
            outputImage.setPixelColor(x, y, QColor(gradientMagnitude, gradientMagnitude, gradientMagnitude));
        }
    }

    return outputImage;
}

QImage QImageAPI::skinImage(const QImage &img)
{
    QImage imgCopy = img;

    QRgb *line;
    for (int y = 0; y < imgCopy.height(); y++) {
        line = (QRgb *)imgCopy.scanLine(y);
        for (int x = 0; x < imgCopy.width(); x++) {
            int R = qRed(line[x]);
            int G = qGreen(line[x]);
            int B = qBlue(line[x]);

            int Y = ((R << 6) + (R << 1) + (G << 7) + (G << 4) + (B << 4) + (B << 3) + 3840) >> 8;
            int Cb = (-((R << 5) + (R << 2) + (R << 1)) - ((G << 6) + (G << 3) + (G << 2)) + ((B << 6) + (R << 5) + (R << 4)) + 32768) >> 8;
            int Cr = (((R << 6) + (R << 5) + (R << 4)) - ((G << 6) + (G << 4) + (G << 3) + (G << 2) + (G << 1)) - ((B << 4) + (B << 1)) + 32768) >> 8;

            if ((Cb > 77 && Cb < 127) && (Cr > 133 && Cr < 173)) {
                //imgCopy.setPixel(x,y, qRgb(0, 0, 0));
            } else
                imgCopy.setPixel(x, y, qRgb(255, 255, 255));
        }

    }
    return imgCopy;
}

// 直方图均衡化函数
QImage QImageAPI::equalizeHistogram(const QImage &inputImage)
{
    // 获取图像的大小
    int width = inputImage.width();
    int height = inputImage.height();

    // 计算像素总数
    int totalPixels = width * height;

    // 计算直方图
    int histogram[256] = {0};
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            QColor pixelColor(inputImage.pixel(x, y));
            int intensity = qRound(0.299 * pixelColor.red() + 0.587 * pixelColor.green() + 0.114 * pixelColor.blue());
            histogram[intensity]++;
        }
    }

    // 计算累积分布函数
    float cumulativeDistribution[256] = {0.0f};
    cumulativeDistribution[0] = static_cast<float>(histogram[0]) / totalPixels;
    for (int i = 1; i < 256; ++i) {
        cumulativeDistribution[i] = cumulativeDistribution[i - 1] + static_cast<float>(histogram[i]) / totalPixels;
    }

    // 对每个像素进行直方图均衡化
    QImage outputImage(inputImage.size(), inputImage.format());
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            QColor pixelColor(inputImage.pixel(x, y));
            int intensity = qRound(0.299 * pixelColor.red() + 0.587 * pixelColor.green() + 0.114 * pixelColor.blue());
            float newIntensity = 255.0f * cumulativeDistribution[intensity];
            outputImage.setPixel(x, y, qRgb(newIntensity, newIntensity, newIntensity));
        }
    }

    return outputImage;
}




// 绘制直方图
QImage QImageAPI::drawHistogram(const QImage &image)
{
    QVector<int> histogram(256, 0);

    // 计算图像的灰度直方图
    for (int y = 0; y < image.height(); ++y) {
        for (int x = 0; x < image.width(); ++x) {
            QColor color(image.pixel(x, y));
            int intensity = qGray(color.rgb());
            histogram[intensity]++;
        }
    }

    int maxValue = *std::max_element(histogram.constBegin(), histogram.constEnd());

    QImage histogramImage(256, maxValue, QImage::Format_RGB32);
    histogramImage.fill(Qt::white);

    QPainter painter(&histogramImage);
    painter.setPen(Qt::black);

    for (int i = 0; i < histogram.size(); ++i) {
        int height = histogram[i] * histogramImage.height() / maxValue;
        painter.drawLine(i, histogramImage.height(), i, histogramImage.height() - height);
    }

    return histogramImage;
}


// 定义一个简单的滤波核，这里使用3x3的均值滤波核
const int filterSize = 3;
const int filter[filterSize][filterSize] = {
    {1, 1, 1},
    {1, 1, 1},
    {1, 1, 1}
};

// 对图像进行均值滤波
QImage QImageAPI::applyMeanFilter(const QImage &inputImage)
{
    QImage outputImage = inputImage;

    for (int y = 1; y < inputImage.height() - 1; ++y) {
        for (int x = 1; x < inputImage.width() - 1; ++x) {
            int sumRed = 0, sumGreen = 0, sumBlue = 0;

            // 遍历滤波核
            for (int i = 0; i < filterSize; ++i) {
                for (int j = 0; j < filterSize; ++j) {
                    QColor color(inputImage.pixel(x + i - 1, y + j - 1));
                    sumRed += color.red() * filter[i][j];
                    sumGreen += color.green() * filter[i][j];
                    sumBlue += color.blue() * filter[i][j];
                }
            }

            // 计算平均值
            int avgRed = sumRed / (filterSize * filterSize);
            int avgGreen = sumGreen / (filterSize * filterSize);
            int avgBlue = sumBlue / (filterSize * filterSize);

            // 将新值设置为输出图像中的像素
            QColor newColor(avgRed, avgGreen, avgBlue);
            outputImage.setPixel(x, y, newColor.rgb());
        }
    }

    return outputImage;
}

// 高斯滤波核大小
const int GSfilterSize = 5;

// 生成高斯滤波核
double generateGaussian(int x, int y, double sigma)
{
    return exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
}

void generateGaussianFilter(double kernel[][GSfilterSize], double sigma)
{
    double sum = 0.0;
    int halfSize = GSfilterSize / 2;

    for (int x = -halfSize; x <= halfSize; ++x) {
        for (int y = -halfSize; y <= halfSize; ++y) {
            kernel[x + halfSize][y + halfSize] = generateGaussian(x, y, sigma);
            sum += kernel[x + halfSize][y + halfSize];
        }
    }

    for (int i = 0; i < GSfilterSize; ++i) {
        for (int j = 0; j < GSfilterSize; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

// 对图像进行高斯滤波
QImage QImageAPI::applyGaussianFilter(const QImage &inputImage, double sigma)
{
    QImage outputImage = inputImage;

    double kernel[GSfilterSize][GSfilterSize];
    generateGaussianFilter(kernel, sigma);

    int halfSize = GSfilterSize / 2;

    for (int y = halfSize; y < inputImage.height() - halfSize; ++y) {
        for (int x = halfSize; x < inputImage.width() - halfSize; ++x) {
            double sumRed = 0.0, sumGreen = 0.0, sumBlue = 0.0;

            for (int i = 0; i < GSfilterSize; ++i) {
                for (int j = 0; j < GSfilterSize; ++j) {
                    QColor color(inputImage.pixel(x + i - halfSize, y + j - halfSize));
                    sumRed += color.red() * kernel[i][j];
                    sumGreen += color.green() * kernel[i][j];
                    sumBlue += color.blue() * kernel[i][j];
                }
            }

            // 更新输出图像中的像素值
            QColor newColor(qBound(0, static_cast<int>(sumRed), 255),
                            qBound(0, static_cast<int>(sumGreen), 255),
                            qBound(0, static_cast<int>(sumBlue), 255));
            outputImage.setPixel(x, y, newColor.rgb());
        }
    }

    return outputImage;
}


// 对图像进行中值滤波
QImage QImageAPI::applyMedianFilter(const QImage &inputImage)
{
    // 中值滤波核大小
    const int filterSize = 3;

    QImage outputImage = inputImage;

    int halfSize = filterSize / 2;

    for (int y = halfSize; y < inputImage.height() - halfSize; ++y) {
        for (int x = halfSize; x < inputImage.width() - halfSize; ++x) {
            // 收集周围像素的颜色值
            QVector<int> redValues, greenValues, blueValues;

            for (int i = 0; i < filterSize; ++i) {
                for (int j = 0; j < filterSize; ++j) {
                    QColor color(inputImage.pixel(x + i - halfSize, y + j - halfSize));
                    redValues.append(color.red());
                    greenValues.append(color.green());
                    blueValues.append(color.blue());
                }
            }

            // 对颜色值进行排序
            std::sort(redValues.begin(), redValues.end());
            std::sort(greenValues.begin(), greenValues.end());
            std::sort(blueValues.begin(), blueValues.end());

            // 选择排序后的中间值作为新的像素值
            int medianRed = redValues.at(redValues.size() / 2);
            int medianGreen = greenValues.at(greenValues.size() / 2);
            int medianBlue = blueValues.at(blueValues.size() / 2);

            // 更新输出图像中的像素值
            QColor newColor(medianRed, medianGreen, medianBlue);
            outputImage.setPixel(x, y, newColor.rgb());
        }
    }

    return outputImage;
}
