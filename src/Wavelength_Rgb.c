#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "Wavelength_Rgb.h"

static double clamp(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static double max3(double a, double b, double c) {
    double m = a > b ? a : b;
    return m > c ? m : c;
}

static double min3(double a, double b, double c) {
    double m = a < b ? a : b;
    return m < c ? m : c;
}


uint32_t wavelengthToHex(double wavelength) {
    double r, g, b;

    if (wavelength >= 380 && wavelength < 440) {
        r = -(wavelength - 440) / (440 - 380);
        g = 0;
        b = 1;
    } else if (wavelength >= 440 && wavelength < 490) {
        r = 0;
        g = (wavelength - 440) / (490 - 440);
        b = 1;
    } else if (wavelength >= 490 && wavelength < 510) {
        r = 0;
        g = 1;
        b = -(wavelength - 510) / (510 - 490);
    } else if (wavelength >= 510 && wavelength < 580) {
        r = (wavelength - 510) / (580 - 510);
        g = 1;
        b = 0;
    } else if (wavelength >= 580 && wavelength < 645) {
        r = 1;
        g = -(wavelength - 645) / (645 - 580);
        b = 0;
    } else if (wavelength >= 645 && wavelength <= 780) {
        r = 1;
        g = 0;
        b = 0;
    } else {
        r = g = b = 0;
    }

    /* Intensity falls off near the edges of the visible range */
    double factor;
    if (wavelength >= 380 && wavelength < 420) {
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
    } else if (wavelength >= 420 && wavelength <= 700) {
        factor = 1;
    } else if (wavelength > 700 && wavelength <= 780) {
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
    } else {
        factor = 0;
    }

    uint32_t R = (uint32_t)lround(clamp(r * factor * 255, 0, 255));
    uint32_t G = (uint32_t)lround(clamp(g * factor * 255, 0, 255));
    uint32_t B = (uint32_t)lround(clamp(b * factor * 255, 0, 255));

    return (R << 24) | (G << 16) | (B << 8) | 0xFFu;
}

void rgbToHsv(uint32_t rgb, int hsv[3]) {
    double r = ((rgb >> 16) & 0xFF) / 255.0;
    double g = ((rgb >> 8)  & 0xFF) / 255.0;
    double b = ( rgb        & 0xFF) / 255.0;

    double max = max3(r, g, b);
    double min = min3(r, g, b);
    double d = max - min;

    double h = 0;
    double s = (max == 0) ? 0 : d / max;
    double v = max;

    if (d != 0) {
        if (max == r)      h = 60 * fmod((g - b) / d, 6);
        else if (max == g) h = 60 * ((b - r) / d + 2);
        else               h = 60 * ((r - g) / d + 4);
    }
    if (h < 0) h += 360;

    hsv[0] = (int)lround(h);
    hsv[1] = (int)lround(s * 100);
    hsv[2] = (int)lround(v * 100);
}


int hexToWavelength(uint32_t hex) {
    int hsv[3];
    rgbToHsv(hex >> 8, hsv);           /* drop the alpha byte */

    if (hsv[1] < 10 || hsv[2] < 10) return -1;   /* grey or too dark: no wavelength */

    double h = hsv[0];
    double wl;

    if      (h <= 60)  wl = 645 - 65 * (h / 60);          /* red -> yellow    */
    else if (h <= 120) wl = 580 - 70 * ((h - 60)  / 60);  /* yellow -> green  */
    else if (h <= 180) wl = 510 - 20 * ((h - 120) / 60);  /* green -> cyan    */
    else if (h <= 240) wl = 490 - 50 * ((h - 180) / 60);  /* cyan -> blue     */
    else if (h <= 300) wl = 440 - 60 * ((h - 240) / 60);  /* blue -> violet   */
    else               wl = (h < 330) ? 380 : 645;        /* magenta/pink: not a spectral colour,
                                                             snap to the nearest end */

    return (int)lround(wl);
}
   