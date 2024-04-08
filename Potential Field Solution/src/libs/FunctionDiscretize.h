#ifndef _LIB_FUNC_DISCR_
#define _LIB_FUNC_DISCR_

#include <math.h>
#include <stdio.h>
#include <string.h>

/// @brief Struct for Multiple Point Representation
typedef struct point_set {
    int n_points; // Current # of Points
    int n_max;    // Max # of Points allocated in memory
    double *i;    // X/i Coordinates
    double *j;    // Y/j Coordinates

    /// @brief Constructor
    /// @param n Initial # of Points
    point_set(int n) {
        this->n_points = 0;
        this->n_max = n;
        this->i = (double *) malloc(sizeof(double) * n);
        this->j = (double *) malloc(sizeof(double) * n);
        memset(this->i, -1, sizeof(double) * n);
        memset(this->j, 0, sizeof(double) * n);
    }

    /// @brief Destructor
    void free_point_set() {
        free(this->i);
        free(this->j);
    }

    /// @brief Get the j Coordinate from the i Coordinate, of the (i,j) Point.
    /// Mimics the functionality of the [] operator
    /// @param i i Coordinate of Point (i,j)
    /// @return j Coordinate of Point (i,j)
    double operator[](int i) {
        if (i < 0)
            return 0;

        for (int k = 0; i < this->n_points; i++) {
            if (this->i[k] == i)
                return this->j[k];
        }
        return 0;
    }

    /// @brief Add Point (i,j) to the Point Set
    /// @param i i Coordinate of Point (i,j)
    /// @param j j Coordinate of Point (i,j)
    /// @return 0 for Success, -1 for Error
    int setPoint(double i, double j) {
        if (i < 0 || j < 0) {
            return -1;
        }
        if (this->n_points >= this->n_max) {
            this->n_max = this->n_points + 1;
            this->i = (double *) realloc(this->i, sizeof(double) * this->n_max);
            this->j = (double *) realloc(this->j, sizeof(double) * this->n_max);
        }
        this->i[this->n_points] = i;
        this->j[this->n_points] = j;
        this->n_points++;
        return 0;
    }
} pointSet;

/// @brief Struct containing information about
/// the Grid size
typedef struct grid_ {
    int L; // Length (X) of Grid (True Units)
    int H; // Height (Y) of Grid (True Units)
    int N; // # of Points in X/i - Direction
    int M; // # of Points in Y/j - Direction

    double Dx;
    double Dy;

    /// @brief Constructor for Square Grid
    /// @param L Side (True Units)
    /// @param N # of Points per Side
    grid_(int L, int N) {
        this->L = L;
        this->H = L;
        this->N = N;
        this->M = N;

        this->Dx = (double) L / (double) N;
        this->Dy = (double) L / (double) N;
    }

    /// @brief Constructor for Rectangular Grid
    /// @param L Length (X) of Grid (True Units)
    /// @param H Height (Y) of Grid (True Units)
    /// @param N # of Points in X/i - Direction
    /// @param M # of Points in Y/j - Direction
    grid_(int L, int H, int N, int M) {
        this->L = L;
        this->H = H;
        this->N = N;
        this->M = M;

        this->Dx = (double) L / (double) N;
        this->Dy = (double) H / (double) M;
    }
} grid;

/// @brief Struct containing information about
/// the function's Domain
typedef struct domain_ {
    // Real Range
    double x_min;
    double x_max;
    double y_min;
    double y_max;

    // Index Range
    int i_min;
    int i_max;
    int j_min;
    int j_max;

    /// @brief Constructor (Only X Coordinate Range needed)
    /// @param func Function pointer to a function of the form y = f(x)
    /// @warning Assumes function's y_min, y_max are at x_min, x_max
    domain_(double x_min, double x_max, double (*func)(double)) {
        this->x_min = x_min;
        this->x_max = x_max;
        this->y_min = func(x_min);
        this->y_max = func(x_max);
    }

    /// @brief Constructor (Full)
    domain_(double x_min, double x_max, double y_min, double y_max) {
        this->x_min = x_min;
        this->x_max = x_max;
        this->y_min = y_min;
        this->y_max = y_max;
    }

    /// @brief Calculates Index Range based on
    /// Real Range and given Grid G
    /// @param G Grid Struct
    void calcIndex(grid G) {
        this->i_min = floor(this->x_min / G.Dx);
        this->i_max = ceil(this->x_max / G.Dx);
        this->j_min = floor(this->y_min / G.Dy);
        this->j_max = ceil(this->y_max / G.Dy);
    }
} domain;

/// @brief Find Pointset P of Points that approximate the function func()
/// on a given grid G and domain D. X-Direction Search Algorithm
/// @param G Grid Struct
/// @param D Domain Struct
/// @param P Point Set
/// @param func Function Pointer
void funcDiscrete(grid G, domain D, pointSet *P, double (*func)(double)) {
    int j = D.j_min;
    for (int i = D.i_min; i < D.i_max; i++) {
        int counter = 0;
        double y = func(i * G.Dx);
        if (y > G.H)
            break;
        if (func((i + 1) * G.Dx) > y) {
            while (fabs((j + 1) * G.Dy - y) < fabs(y - j * G.Dy)) {
                j++, counter++;
                if (j + 1 > G.M)
                    break;
            }
        } else {
            while (fabs((j - 1) * G.Dy - y) < fabs(y - j * G.Dy)) {
                j--, counter++;
                if (j - 1 < 0)
                    break;
            }
        }

        if (counter > 1 && i + 1 > 1) {
            for (int k = 1; k <= floor(counter / 2.0); k++) {
                P->setPoint(i - 1, j - counter + k); // Filler Points at i-1
            }
            for (int k = floor(counter / 2.0) + 1; k <= counter - 1; k++) {
                P->setPoint(i, j - counter + k); // Filler Points at i+1
            }
        }
        P->setPoint(i, j); // Actual Point of function
    }
}

/// @brief Find Pointset P of Points that approximate the linear function func_lin()
/// on a given grid G and domain D.
/// @param G Grid Struct
/// @param D Domain Struct
/// @param P Point Set
/// @param func_lin Function Pointer
/// @param angle Angle of Linear Function. Determines direction of search algorith
void funcDiscrete(grid G, domain D, pointSet *P, double (*func_lin)(double), double angle) {
    if (angle <= 45) {
        funcDiscrete(G, D, P, func_lin);
    } else {
        int i = D.i_min;
        for (int j = D.j_min; j < D.j_max; j++) {
            int counter = 0;
            double x = D.x_min + (((j * G.Dy) - D.y_min) * (D.x_max - D.x_min)) / (func_lin(D.x_max) - D.y_min);
            printf("%lf %lf\n",x,(j * G.Dy));
            if (x > G.L)
                break;
            double x_next = D.x_min + ((((j + 1) * G.Dy) - D.y_min) * (D.x_max - D.x_min)) / (D.y_max - D.y_min);
            if (x_next > x) {
                while (fabs((i + 1) * G.Dx - x) < fabs(x - i * G.Dx)) {
                    i++, counter++;
                    if (i + 1 > G.N)
                        break;
                }
            } else {
                while (fabs((i - 1) * G.Dx - x) < fabs(x - i * G.Dx)) {
                    i--, counter++;
                    if (i - 1 < 0)
                        break;
                }
            }

            if (counter > 1 && j + 1 > 1) {
                for (int k = 1; k <= floor(counter / 2.0); k++) {
                    P->setPoint(i - counter + k, j - 1); // Filler Points at i-1
                }
                for (int k = floor(counter / 2.0) + 1; k <= counter - 1; k++) {
                    P->setPoint(i - counter + k, j); // Filler Points at i+1
                }
            }
            P->setPoint(i, j); // Actual Point of function
        }
    }
}

#endif // !_LIB_FUNC_DISCR_