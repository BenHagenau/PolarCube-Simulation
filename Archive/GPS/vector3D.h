/*
 *  vector3D.h
 *
 *  Created by Hanspeter Schaub on Sat Mar 01 2003.
 *  Copyright (c) 2003 Hanspeter Schaub. All rights reserved.
 */
#ifndef _VECTOR_3D_H_
#define _VECTOR_3D_H_

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

    void    diag(double v[3], double m[3][3]);
    void    set3(double v1, double v2, double v3, double result[3]);
    void    setZero(double v[3]);
    void    equal(double v1[3], double v2[3]);
    void    add(double v1[3], double v2[3], double result[3]);
    void    sub(double v1[3], double v2[3], double result[3]);
    void    mult(double scaleFactor, double v[3], double result[3]);
    void    cross(double v1[3], double v2[3], double result[3]);
    double  dot(double v1[3], double v2[3]);
    void    dotT(double v1[3], double v2[3], double result[3][3]);
    double  norm(double v[3]);
    int     equalCheck(double v1[3], double v2[3], double precisionExponent);
    int     VequalCheck(double v1[3], double v2[3], double precisionExponent);
    void    tilde(double v[3], double result[3][3]);
    void    printVector(const char *str, double v[3]);

    void    set4(double v1, double v2, double v3, double v4, double result[5]);
    void    setZero4(double v[5]);
    void    equal4(double v1[5], double v2[5]);
    double  norm4(double v[3]);
    void    printVector4(const char *str, double v[5]);
    int     V4equalCheck(double v1[5], double v2[5], double precisionExponent);

    void    setMatrix(double m11, double m12, double m13,
                      double m21, double m22, double m23,
                      double m31, double m32, double m33,
                      double result[3][3]);
    void    m33SetZero(double m[3][3]);
    void    eye(double m[3][3]);
    void    Madd(double m1[3][3], double m2[3][3], double m3[3][3]);
    void    Msub(double m1[3][3], double m2[3][3], double m3[3][3]);
    void    MdotM(double m1[3][3], double m2[3][3], double ans[3][3]);
    void    Mdot(double m[3][3], double v[3], double result[3]);
    void    MdotMT(double m1[3][3], double m2[3][3], double ans[3][3]);
    void    Mmult(double, double m[3][3], double m2[3][3]);
    void    transpose(double m[3][3], double ans[3][3]);
    void    inverse(double m[3][3], double ans[3][3]);
    double  detM(double a[3][3]);
    double  trace(double m[3][3]);
    int     MequalCheck(double m1[3][3], double m2[3][3], double d);
    void    printMatrix(const char *str, double m[3][3]);

    double  arc_cosh(double);
    double  arc_tanh(double);

	void	  getColumn(float m[3][3], int col, float ans[3]);

    /* common variables */
#ifndef M_PI
#define M_PI        3.141592653589793
#endif
#ifndef D2R
#define D2R         M_PI/180.
#endif
#ifndef R2D
#define R2D         180./M_PI
#endif
#ifndef NAN
#define NAN         sqrt((float)-1)
#endif

#ifdef __cplusplus
}
#endif

#endif
