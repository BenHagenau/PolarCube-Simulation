/*
 *  vector3D.c
 *
 *  Created by Hanspeter Schaub on Sat Mar 01 2003.
 *  Copyright (c) 2003 Hanspeter Schaub. All rights reserved.
 *
 *  Provides a library for doing 3D matrix algebra.
 */

#include "vector3D.h"

void diag(double v[3], double m[3][3])
{
    m33SetZero(m);
    m[0][0] = v[0];
    m[1][1] = v[1];
    m[2][2] = v[2];
}

void set3(double v1, double v2, double v3, double result[3])
{
    result[0] = v1;
    result[1] = v2;
    result[2] = v3;
}

void setZero(double v[3])
{
    v[0] = 0.;
    v[1] = 0.;
    v[2] = 0.;
}

void equal(double v1[3], double v2[3])
{
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];

    return;
}

void add(double v1[3], double v2[3], double result[3])
{
    result[0] = v1[0] + v2[0];
    result[1] = v1[1] + v2[1];
    result[2] = v1[2] + v2[2];
}

void sub(double v1[3], double v2[3], double result[3])
{
    result[0] = v1[0] - v2[0];
    result[1] = v1[1] - v2[1];
    result[2] = v1[2] - v2[2];
}

void mult(double scaleFactor, double v[3], double result[3])
{
    result[0] = scaleFactor * v[0];
    result[1] = scaleFactor * v[1];
    result[2] = scaleFactor * v[2];
}

void cross(double v1[3], double v2[3], double result[3])
{
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double dot(double v1[3], double v2[3])
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void dotT(double v1[3], double v2[3], double result[3][3])
{
    mult(v1[0], v2, result[0]);
    mult(v1[1], v2, result[1]);
    mult(v1[2], v2, result[2]);
}

double norm(double v[3])
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

int equalCheck(double v1[3], double v2[3], double precisionExponent)
{
    int c = 1, i;

    for(i = 0; i <= 2; i++) {
        if(fabs(v1[i] - v2[i]) > pow(10., -precisionExponent)) {
            c = 0;
        }
    }

    return c;
}

int VequalCheck(double v1[3], double v2[3], double precisionExponent)
{
    int c = 1, i;

    for(i = 0; i <= 2; i++) {
        if(fabs(v1[i] - v2[i]) > pow(10., -precisionExponent)) {
            c = 0;
        }
    }

    return c;
}

void tilde(double v[3], double result[3][3])
{
    result[0][0] = result[1][1] = result[2][2] = 0.;
    result[0][1] = -v[2];
    result[1][0] =  v[2];
    result[0][2] =  v[1];
    result[2][0] = -v[1];
    result[1][2] = -v[0];
    result[2][1] =  v[0];

    return;
}

void printVector(const char *str, double v[3])
{
    printf("%s (%20.15g, %20.15g, %20.15g)\n", str, v[0], v[1], v[2]);

    return;
}

void set4(double v1, double v2, double v3, double v4, double result[5])
{
    result[0] = v1;
    result[1] = v2;
    result[2] = v3;
    result[3] = v4;
}

void setZero4(double v[5])
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 0.0;
}

void equal4(double v1[5], double v2[5])
{
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];
    v2[3] = v1[3];

    return;
}

double norm4(double v[5])
{
    double result;

    result = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);

    return result;
}

void printVector4(const char *str, double v[5])
{
    printf("%s (%20.15g, %20.15g, %20.15g, %20.15g)\n",
           str, v[0], v[1], v[2], v[3]);

    return;
}

int V4equalCheck(double v1[5], double v2[5], double precisionExponent)
{
    int c = 1, i;

    for(i = 0; i <= 2; i++) {
        if(fabs(v1[i] - v2[i]) > pow(10., -precisionExponent)) {
            c = 0;
        }
    }

    return c;
}

void setMatrix(double m11, double m12, double m13,
               double m21, double m22, double m23,
               double m31, double m32, double m33, double result[3][3])
{
    result[0][0] = m11;
    result[0][1] = m12;
    result[0][2] = m13;
    result[1][0] = m21;
    result[1][1] = m22;
    result[1][2] = m23;
    result[2][0] = m31;
    result[2][1] = m32;
    result[2][2] = m33;
}

void m33SetZero(double m[3][3])
{
    int i;
    int j;
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            m[i][j] = 0.0;
        }
    }
}

void eye(double m[3][3])
{
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;

    return;
}

void Madd(double m1[3][3], double m2[3][3], double m3[3][3])
{
    add(m1[0], m2[0], m3[0]);
    add(m1[1], m2[1], m3[1]);
    add(m1[2], m2[2], m3[2]);

    return;
}

void Msub(double m1[3][3], double m2[3][3], double m3[3][3])
{
    sub(m1[0], m2[0], m3[0]);
    sub(m1[1], m2[1], m3[1]);
    sub(m1[2], m2[2], m3[2]);

    return;
}

void MdotM(double m1[3][3], double m0[3][3], double ans[3][3])
{
    double m2[3][3];
    transpose(m0, m2);

    ans[0][0] = dot(m1[0], m2[0]);
    ans[0][1] = dot(m1[0], m2[1]);
    ans[0][2] = dot(m1[0], m2[2]);
    ans[1][0] = dot(m1[1], m2[0]);
    ans[1][1] = dot(m1[1], m2[1]);
    ans[1][2] = dot(m1[1], m2[2]);
    ans[2][0] = dot(m1[2], m2[0]);
    ans[2][1] = dot(m1[2], m2[1]);
    ans[2][2] = dot(m1[2], m2[2]);

    return;
}

void Mdot(double m[3][3], double v[3], double result[3])
{
    result[0] = dot(m[0], v);
    result[1] = dot(m[1], v);
    result[2] = dot(m[2], v);
}

void MdotMT(double m1[3][3], double m0[3][3], double ans[3][3])
{
    ans[0][0] = dot(m1[0], m0[0]);
    ans[0][1] = dot(m1[0], m0[1]);
    ans[0][2] = dot(m1[0], m0[2]);
    ans[1][0] = dot(m1[1], m0[0]);
    ans[1][1] = dot(m1[1], m0[1]);
    ans[1][2] = dot(m1[1], m0[2]);
    ans[2][0] = dot(m1[2], m0[0]);
    ans[2][1] = dot(m1[2], m0[1]);
    ans[2][2] = dot(m1[2], m0[2]);

    return;
}

void Mmult(double a, double m[3][3], double m2[3][3])
{
    mult(a, m[0], m2[0]);
    mult(a, m[1], m2[1]);
    mult(a, m[2], m2[2]);

    return;
}

void transpose(double m1[3][3], double ans[3][3])
{
    ans[0][0] = m1[0][0];
    ans[0][1] = m1[1][0];
    ans[0][2] = m1[2][0];
    ans[1][0] = m1[0][1];
    ans[1][1] = m1[1][1];
    ans[1][2] = m1[2][1];
    ans[2][0] = m1[0][2];
    ans[2][1] = m1[1][2];
    ans[2][2] = m1[2][2];

    return;
}
void inverse(double m[3][3], double ans[3][3])
{
    double det;

    det = detM(m);

    //if(fabs(det) > 1e-15) {
    ans[0][0] = (-m[1][2] * m[2][1] + m[1][1] * m[2][2]) / det;
    ans[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
    ans[0][2] = (-m[0][2] * m[1][1] + m[0][1] * m[1][2]) / det;
    ans[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
    ans[1][1] = (-m[0][2] * m[2][0] + m[0][0] * m[2][2]) / det;
    ans[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
    ans[2][0] = (-m[1][1] * m[2][0] + m[1][0] * m[2][1]) / det;
    ans[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
    ans[2][2] = (-m[0][1] * m[1][0] + m[0][0] * m[1][1]) / det;
    //} else {
    //    fprintf(stderr, "ERROR: singular 3x3 matrix inverse\n");
    //    ans[0][0] = ans[0][1] = ans[0][2] = ans[1][0] = ans[1][1] = ans[1][2] =
    //                ans[2][0] = ans[2][1] = ans[2][2] = NAN;
    //}

    return;
}

double detM(double m[3][3])
{
    return -m[0][2] * m[1][1] * m[2][0] + m[0][1] * m[1][2] * m[2][0]
            + m[0][2] * m[1][0] * m[2][1] - m[0][0] * m[1][2] * m[2][1]
            - m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][1] * m[2][2];
}

double trace(double m[3][3])
{
    return m[0][0] + m[1][1] + m[2][2];
}

int MequalCheck(double m1[3][3], double m2[3][3], double d)
{
    int c = 1, i, j;

    for(i = 0; i <= 2; i++) {
        for(j = 0; j <= 2; j++) {
            if(fabs(m1[i][j] - m2[i][j]) > pow(10., -d)) {
                c = 0;
            }
        }
    }

    return c;
}

void printMatrix(const char *str, double m[3][3])
{
    printf("%s:\n%20.15g, %20.15g, %20.15g\n", str, m[0][0], m[0][1], m[0][2]);
    printf("%20.15g, %20.15g, %20.15g\n",           m[1][0], m[1][1], m[1][2]);
    printf("%20.15g, %20.15g, %20.15g\n",           m[2][0], m[2][1], m[2][2]);

    return;
}

double arc_cosh(double z)
{
    return log(z + sqrt(z + 1.) * sqrt(z - 1.));
}

double arc_tanh(double z)
{
    return 0.5 * (log(1. + z) - log(1. - z));
}

void getColumn(float m[3][3], int col, float ans[3])
{
	ans[0] = m[0][col];
	ans[1] = m[1][col];
	ans[2] = m[2][col];
}
