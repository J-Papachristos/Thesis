#include "libs/FunctionDiscretize.h"
#include "libs/Sparse.h"

/// Utility Library from Potential Field Solution Code

#ifndef _LIB_POTENTIAL_SOLVE
#define _LIB_POTENTIAL_SOLVE

#define X_SIZE (sizeof(double) * N * N)
#define A_SIZE (X_SIZE * N * N)
#define deg2rad(X) ((X) * (M_PI / 180.0))

#define cosd(X) (cos(deg2rad(X)))
#define sind(X) (sin(deg2rad(X)))
#define tand(X) (tan(deg2rad(X)))

// Crack Line Function
// https://www.desmos.com/calculator/rthhj6yoh6

// Domain : x_c0 - (L/2)*cos(a) - (t/2)*sin(a) : x_c0 + (L/2)*cos(a) - (t/2)*sin(a)
double crackLineLeft(double x) {
    return (y_c0 + (crack_thick / cosd(crack_angle) / 2.0)) + tand(crack_angle) * (x - x_c0);
}
// Domain : x_c0 - (L/2)*cos(a) + (t/2)*sin(a) : x_c0 + (L/2)*cos(a) + (t/2)*sin(a)
double crackLineRight(double x) {
    return (y_c0 - (crack_thick / cosd(crack_angle) / 2.0)) + tand(crack_angle) * (x - x_c0);
}
// Domain : x_c0 + (L/2)*cos(a) - (t/2)*sin(a) : x_c0 + (L/2)*cos(a) + (t/2)*sin(a)
double crackLineTop(double x) {
    return y_c0 + (crack_len / 2.0 * sind(crack_angle)) + tand(crack_angle + 90) * (x - x_c0 - crack_len / 2.0 * cosd(crack_angle));
}
// Domain : x_c0 - (L/2)*cos(a) - (t/2)*sin(a) : x_c0 - (L/2)*cos(a) + (t/2)*sin(a)
double crackLineBottom(double x) {
    return y_c0 - (crack_len / 2.0 * sind(crack_angle)) + tand(crack_angle + 90) * (x - x_c0 + crack_len / 2.0 * cosd(crack_angle));
}

enum sideEnum {
    LEFT = 0,
    RIGHT = 1,
    TOP = 2,
    BOTTOM = 3,
    SIDE_SIZE
};

#endif