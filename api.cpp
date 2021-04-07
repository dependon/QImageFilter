#include "api.h"


int api::Bound(int range_left, int data, int range_right)
{
    int index = data;
    if (data > range_right) {
        index = range_right;
    } else if (data < range_left) {
        index = range_left;
    }
    return index;
}

void api::RGBToYUV(int Red, int Green, int Blue, int *Y, int *U, int *V)

{

    *Y = ((Red << 6) + (Red << 3) + (Red << 2) + Red + (Green << 7) + (Green << 4) + (Green << 2) + (Green << 1) + (Blue << 4) + (Blue << 3) + (Blue << 2) + Blue) >> 8;

    *U = (-((Red << 5) + (Red << 2) + (Red << 1)) - ((Green << 6) + (Green << 3) + (Green << 1)) + ((Blue << 6) + (Blue << 5) + (Blue << 4))) >> 8;

    *V = ((Red << 7) + (Red << 4) + (Red << 3) + (Red << 2) + (Red << 1) - ((Green << 7) + (Green << 2)) - ((Blue << 4) + (Blue << 3) + (Blue << 1))) >> 8;

}

void api::YUVToRGB(int Y, int U, int V, int *Red, int *Green, int *Blue)

{

    *Red   = ((Y << 8) + ((V << 8) + (V << 5) + (V << 2))) >> 8;

    *Green = ((Y << 8) - ((U << 6) + (U << 5) + (U << 2)) - ((V << 7) + (V << 4) + (V << 2) + V)) >> 8;

    *Blue = ((Y << 8) + (U << 9) + (U << 3)) >> 8;

}

void api::RunBEEPSVerticalHorizontal(double *data, int width, int height, double spatialDecay, double *exp_table, double *g_table)
{
    int length0 = height * width;
    double *g = new double[length0];
    int m = 0;
    for (int k2 = 0; k2 < height; ++k2) {
        int n = k2;
        for (int k1 = 0; k1 < width; ++k1) {
            g[n] = data[m++];
            n += height;
        }
    }
    double *p = new double[length0];
    double *r = new double[length0];
    memcpy(p, g, sizeof(double) * length0);
    memcpy(r, g, sizeof(double) * length0);
    for (int k1 = 0; k1 < width; ++k1) {
        int startIndex = k1 * height;
        double mu = 0.0;
        for (int k = startIndex + 1, K = startIndex + height; k < K; ++k) {
            int div0 = fabs(p[k] - p[k - 1]);
            mu = exp_table[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);//文献中的公式1，这里做了一下修改，效果影响不大
        }
        for (int k = startIndex + height - 2; startIndex <= k; --k) {
            int div0 = fabs(r[k] - r[k + 1]);
            mu = exp_table[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu) ; //文献公式3
        }
    }
    double rho0 = 1.0 / (2 - spatialDecay);
    for (int k = 0; k < length0; ++k) {
        r[k] = (r[k] + p[k]) * rho0 - g_table[(int)g[k]];
    }
    m = 0;
    for (int k1 = 0; k1 < width; ++k1) {
        int n = k1;
        for (int k2 = 0; k2 < height; ++k2) {
            data[n] = r[m++];
            n += width;
        }
    }

    memcpy(p, data, sizeof(double) * length0);
    memcpy(r, data, sizeof(double) * length0);
    for (int k2 = 0; k2 < height; ++k2) {
        int startIndex = k2 * width;
        double mu = 0.0;
        for (int k = startIndex + 1, K = startIndex + width; k < K; ++k) {
            int div0 = fabs(p[k] - p[k - 1]);
            mu = exp_table[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);
        }
        for (int k = startIndex + width - 2; startIndex <= k; --k) {
            int div0 = fabs(r[k] - r[k + 1]);
            mu = exp_table[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu) ;
        }
    }

    double init_gain_mu = spatialDecay / (2 - spatialDecay);
    for (int k = 0; k < length0; k++) {
        data[k] = (p[k] + r[k]) * rho0 - data[k] * init_gain_mu; //文献中的公式5
    }
    delete[] p;
    delete[] r;
    delete[] g;
}

void api::RunBEEPSHorizontalVertical(double *data, int width, int height, double spatialDecay, double *exptable, double *g_table)
{
    int length = width * height;
    double *g = new double[length];
    double *p = new double[length];
    double *r = new double[length];
    memcpy(p, data, sizeof(double) * length);
    memcpy(r, data, sizeof(double) * length);
    double rho0 = 1.0 / (2 - spatialDecay);
    for (int k2 = 0; k2 < height; ++k2) {
        int startIndex = k2 * width;
        for (int k = startIndex + 1, K = startIndex + width; k < K; ++k) {
            int div0 = fabs(p[k] - p[k - 1]);
            double mu = exptable[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);//文献公式1

        }

        for (int k = startIndex + width - 2; startIndex <= k; --k) {
            int div0 = fabs(r[k] - r[k + 1]);
            double mu = exptable[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu);//文献公式3
        }
        for (int k = startIndex, K = startIndex + width; k < K; k++) {
            r[k] = (r[k] + p[k]) * rho0 - g_table[(int)data[k]];
        }
    }

    int m = 0;
    for (int k2 = 0; k2 < height; k2++) {
        int n = k2;
        for (int k1 = 0; k1 < width; k1++) {
            g[n] = r[m++];
            n += height;
        }
    }

    memcpy(p, g, sizeof(double) * height * width);
    memcpy(r, g, sizeof(double) * height * width);
    for (int k1 = 0; k1 < width; ++k1) {
        int startIndex = k1 * height;
        double mu = 0.0;
        for (int k = startIndex + 1, K = startIndex + height; k < K; ++k) {
            int div0 = fabs(p[k] - p[k - 1]);
            mu = exptable[div0];
            p[k] = p[k - 1] * mu + p[k] * (1.0 - mu);
        }
        for (int k = startIndex + height - 2; startIndex <= k; --k) {
            int div0 = fabs(r[k] - r[k + 1]);
            mu = exptable[div0];
            r[k] = r[k + 1] * mu + r[k] * (1.0 - mu);
        }
    }

    double init_gain_mu = spatialDecay / (2 - spatialDecay);
    for (int k = 0; k < length; ++k) {
        r[k] = (r[k] + p[k]) * rho0 - g[k] * init_gain_mu;
    }
    m = 0;
    for (int k1 = 0; k1 < width; ++k1) {
        int n = k1;
        for (int k2 = 0; k2 < height; ++k2) {
            data[n] = r[m++];
            n += width;
        }
    }

    delete p;
    delete r;
    delete g;
}

void api::QImageD_RunBEEPSHorizontalVertical(QImage *img, QImage *imgCopy, double spatialDecay, double photometricStandardDeviation)
{
    if (!img || !imgCopy) {
        return ;
    }

    double c = -0.5 / (photometricStandardDeviation * photometricStandardDeviation); //-1/2 *光度标准偏差的平方
    double mu = spatialDecay / (2 - spatialDecay);

    double *exptable = new double[256];;
    double *g_table = new double[256];;
    for (int i = 0; i <= 255; i++) {
        exptable[i] = (1 - spatialDecay) * exp(c * i * i);
        g_table[i] = mu * i;
    }
    int width = img->width();
    int height = img->height();
    int length = width * height;
    double *data2Red = new double[length];
    double *data2Green = new double[length];
    double *data2Blue = new double[length];

    int i = 0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            QRgb rgb = imgCopy->pixel(x, y);
            data2Red[i] = qRed(rgb);
            data2Green[i] = qGreen(rgb);
            data2Blue[i] = qBlue(rgb);
            i++;
        }
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
            pRed[k] = pRed[k - 1] * mu + pRed[k] * (1.0 - mu);//公式1

            int div0Green = fabs(pGreen[k] - pGreen[k - 1]);
            mu = exptable[div0Green];
            pGreen[k] = pGreen[k - 1] * mu + pGreen[k] * (1.0 - mu);//公式1

            int div0Blue = fabs(pBlue[k] - pBlue[k - 1]);
            mu = exptable[div0Blue];
            pBlue[k] = pBlue[k - 1] * mu + pBlue[k] * (1.0 - mu);//公式1

        }

        for (int k = startIndex + width - 2; startIndex <= k; --k) {
            int div0Red = fabs(rRed[k] - rRed[k + 1]);
            double mu = exptable[div0Red];
            rRed[k] = rRed[k + 1] * mu + rRed[k] * (1.0 - mu);//公式3

            int div0Green = fabs(rGreen[k] - rGreen[k + 1]);
            mu = exptable[div0Green];
            rGreen[k] = rGreen[k + 1] * mu + rGreen[k] * (1.0 - mu);//公式3

            int div0Blue = fabs(rBlue[k] - rBlue[k + 1]);
            mu = exptable[div0Blue];
            rBlue[k] = rBlue[k + 1] * mu + rBlue[k] * (1.0 - mu);//公式3
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
    for (int k1 = 0; k1 < width; ++k1) {
        int n = k1;
        for (int k2 = 0; k2 < height; ++k2) {

            data2Red[n] = rRed[m];
            data2Green[n] = rGreen[m];
            data2Blue[n] = rBlue[m];
            imgCopy->setPixel(k1, k2, qRgb(data2Red[n], data2Green[n], data2Blue[n]));
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

void api::warnImage(QImage *img, QImage *imgCopy, int index)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    QColor frontColor;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            frontColor = line[x];
            float r = frontColor.red() + index;
            float g = frontColor.green() + index;
            float b = frontColor.blue();
            r = Bound(0, r, 255);
            g = Bound(0, g, 255);
            imgCopy->setPixel(x, y, qRgb(r, g, b));
        }

    }
}

void api::coolImage(QImage *img, QImage *imgCopy, int index)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    QColor frontColor;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            frontColor = line[x];
            float r = frontColor.red();
            float g = frontColor.green();
            float b = frontColor.blue() + index;
            b = Bound(0, b, 255);
            imgCopy->setPixel(x, y, qRgb(r, g, b));
        }

    }
}

void api::GrayScaleImage(QImage *img, QImage *imgCopy)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            int average = (qRed(line[x]) + qGreen(line[x]) + qBlue(line[x])) / 3;
            imgCopy->setPixel(x, y, qRgb(average, average, average));
        }

    }
}

void api::lightContrastImage(QImage *img, QImage *imgCopy, int light, int Contrast)
{
    if (!img || !imgCopy) {
        return ;
    }
    //        *imgCopy= (imgCopy->scaled(1920,1080));
    QRgb *line;
    float lightf = light * 0.01 * 256;
    int Blight = lightf;
    int Bli = light * 0.01;
    int BContrast = Contrast - 150;
    int B1 = -1, B2 = -1, B3 = -1, B4 = -1, B5 = -1, B6 = -1, B7 = -1, B8 = -1;
    int n = 0;
    int num = 0;
    while (Blight > 0) {
        if (Blight & 1) {
            if (B1 < 0) {
                B1 = n;
            } else if (B2 < 0) {
                B2 = n;
            } else if (B3 < 0) {
                B3 = n;
            } else if (B4 < 0) {
                B4 = n;
            } else if (B5 < 0) {
                B5 = n;
            } else if (B6 < 0) {
                B6 = n;
            } else if (B7 < 0) {
                B7 = n;
            } else if (B8 < 0) {
                B8 = n;
            }
        }
        n++;
        Blight = Blight >> 1;

    }
    if (B8 < 0) {
        B8 = 0;
    }
    if (B7 < 0) {
        B7 = 0;
    }
    if (B6 < 0) {
        B6 = 0;
    }
    if (B5 < 0) {
        B5 = 0;
    }
    if (B4 < 0) {
        B4 = 0;
    }
    if (B3 < 0) {
        B3 = 0;
    }
    if (B2 < 0) {
        B2 = 0;
    }
    if (B1 < 0) {
        B1 = 0;
    }
    int index = 0;

    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            int lastR = qRed(line[x]);
            int lastG = qGreen(line[x]);
            int lastB = qBlue(line[x]);

            int r = (((lastR << B1) + (lastR << B2) + (lastR << B3) + (lastR << B4) + (lastR << B5) + (lastR << B6) + (lastR << B7) + (lastR << B8)) >> 8) + BContrast;
            int g = (((lastG << B1) + (lastG << B2) + (lastG << B3) + (lastG << B4) + (lastG << B5) + (lastG << B6) + (lastG << B7) + (lastG << B8)) >> 8) + BContrast;
            int b = (((lastB << B1) + (lastB << B2) + (lastB << B3) + (lastB << B4) + (lastB << B5) + (lastB << B6) + (lastB << B7) + (lastB << B8)) >> 8) + BContrast;

            r = Bound(0, r, 255);
            g = Bound(0, g, 255);
            b = Bound(0, b, 255);
            imgCopy->setPixel(x, y, qRgb(r, g, b));
        }

    }




}

void api::InverseColorImage(QImage *img, QImage *imgCopy)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {

            imgCopy->setPixel(x, y, qRgb(255 - qRed(line[x]), 255 - qGreen(line[x]), 255 - qBlue(line[x])));
        }

    }
}

void api::oldImage(QImage *img, QImage *imgCopy)
{
    if (!img || !imgCopy) {
        return ;
    }

    QRgb *line;
    for (int y = 0; y < imgCopy->height(); y++) {
        line = (QRgb *)imgCopy->scanLine(y);
        for (int x = 0; x < imgCopy->width(); x++) {
            float r = 0.393 * qRed(line[x]) + 0.769 * qGreen(line[x]) + 0.189 * qBlue(line[x]);
            float g = 0.349 * qRed(line[x]) + 0.686 * qGreen(line[x]) + 0.168 * qBlue(line[x]);
            float b = 0.272 * qRed(line[x]) + 0.534 * qGreen(line[x]) + 0.131 * qBlue(line[x]);
            r = Bound(0, r, 255);
            g = Bound(0, g, 255);
            b = Bound(0, b, 255);
            imgCopy->setPixel(x, y, qRgb(r, g, b));
        }

    }
}
