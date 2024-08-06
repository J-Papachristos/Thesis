#include "libs/FunctionDiscretize.h"
#include "libs/Sparse.h"

#define X_SIZE (sizeof(double) * N * N)
#define A_SIZE (X_SIZE * N * N)
#define b2Mb(X) (((X) / 1024) / 1024)
#define deg2rad(X) ((X) * (M_PI / 180.0))

#define cosd(X) (cos(deg2rad(X)))
#define sind(X) (sin(deg2rad(X)))
#define tand(X) (tan(deg2rad(X)))

#define L (400.0 * 1e-3)  // Lengt              [m]
#define H (250.0 * 1e-3)  // Height             [m]
#define t (3.100 * 1e-3)  // Thickness          [m]
#define r (30.000 * 1e-3) // Source/Sink Radius [m]

// Crack Data (Global for Function use)
double x_c0, y_c0;  // [m]
double crack_angle; // [°]
double crack_len;   // [m]
double crack_thick; // [m]

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

int main(int argc, char const *argv[]) {
    /// Mesh Data :
    int N = atoi(argv[1]);
    double Dx = L / N, Dx2 = Dx * Dx;
    int M = argc > 2 ? atoi(argv[2]) : (int) ((H / L) * N);
    double Dy = H / M, Dy2 = Dy * Dy;
    double Dx2Dy2 = Dx2 * Dy2;

    // Create Grid Struct
    grid G(L, H, N, M);

    // Data Structure Initialization
    Sparse *A_sp = (Sparse *) malloc(N * M * sizeof(Sparse));
    for (int i = 0; i < N * M; i++) {
        A_sp[i] = {NNZ_INIT, N * M};
    }
    double *b = (double *) malloc(X_SIZE);
    double *x = (double *) malloc(X_SIZE);
    memset(b, 0, X_SIZE);
    memset(x, 0, X_SIZE);

    /// File Input :
    FILE *fp_init = fopen("input.dat", "r+");
    double I;
    fscanf(fp_init, "I %lf\n", &I);
    // Source/Sink Data :
    I = I / Dx / Dy;                      // [A/m^2]
    double sigma = 7476688.272623353 * t; // [siemens/m^2] * [m] = [A/Vm^2]
    double source = I / sigma;            // [A/m^2] / [A/Vm^2] = [V]

    int n_sources, n_sinks;
    double x_source[1], y_source[1];
    double x_sink[1], y_sink[1];
    int i_source[1], j_source[1];
    int i_sink[1], j_sink[1];

    fscanf(fp_init, "n_sources %d", &n_sources);
    if (!(n_sources > 0)) {
        return -1;
    }
    for (int i = 0; i < n_sources; i++) {
        fscanf(fp_init, "\tx %lf y %lf\n", &x_source[i], &y_source[i]);
        i_source[i] = x_source[i] / Dx, j_source[i] = y_source[i] / Dy;
    }
    fscanf(fp_init, "n_sinks %d", &n_sinks);
    if (!(n_sinks > 0)) {
        return -1;
    }
    for (int i = 0; i < n_sinks; i++) {
        fscanf(fp_init, "\tx %lf y %lf\n", &x_sink[i], &y_sink[i]);
        i_sink[i] = x_sink[i] / Dx, j_sink[i] = y_sink[i] / Dy;
    }

    /// Crack Data :
    int is_cracked;
    fscanf(fp_init, "is_cracked %d\n", &is_cracked);
    if (is_cracked) {
        fscanf(fp_init, "crack_pos %lf %lf\n", &x_c0, &y_c0);
        fscanf(fp_init, "crack_len %lf\n", &crack_len);
        fscanf(fp_init, "crack_thick %lf\n", &crack_thick);
        fscanf(fp_init, "crack_angle %lf\n", &crack_angle);
    } else {
        x_c0 = y_c0 = 0;
        crack_len = crack_thick = 0;
    }
    int i_c0 = x_c0 / Dx, j_c0 = y_c0 / Dy;

    int i_min_crack = i_c0 - (int) ceil(crack_thick / 2.0 / Dx);
    int i_max_crack = i_c0 + (int) ceil(crack_thick / 2.0 / Dx);

    int j_min_crack = j_c0 - (int) ceil(crack_len / 2.0 / Dy);
    int j_max_crack = j_c0 + (int) ceil(crack_len / 2.0 / Dy);

    domain D_left(x_c0 - (cosd(crack_angle) * crack_len / 2.0) - (sind(crack_angle) * crack_thick / 2.0),
                  x_c0 + (cosd(crack_angle) * crack_len / 2.0) - (sind(crack_angle) * crack_thick / 2.0),
                  &crackLineLeft, true);
    D_left.calcIndex(G);
    domain D_right(x_c0 - (cosd(crack_angle) * crack_len / 2.0) + (sind(crack_angle) * crack_thick / 2.0),
                   x_c0 + (cosd(crack_angle) * crack_len / 2.0) + (sind(crack_angle) * crack_thick / 2.0),
                   &crackLineRight, true);
    D_right.calcIndex(G);
    domain D_top(x_c0 + (cosd(crack_angle) * crack_len / 2.0) - (sind(crack_angle) * crack_thick / 2.0),
                 x_c0 + (cosd(crack_angle) * crack_len / 2.0) + (sind(crack_angle) * crack_thick / 2.0),
                 &crackLineTop, false);
    D_top.calcIndex(G);
    domain D_bottom(x_c0 - (cosd(crack_angle) * crack_len / 2.0) - (sind(crack_angle) * crack_thick / 2.0),
                    x_c0 - (cosd(crack_angle) * crack_len / 2.0) + (sind(crack_angle) * crack_thick / 2.0),
                    &crackLineBottom, false);
    D_bottom.calcIndex(G);

    int init_points_long, init_points_short;
    if (crack_angle > 45) {
        init_points_long = MAX(D_left.j_max - D_left.j_min, D_right.j_max - D_right.j_min);
        init_points_short = MAX(D_top.i_max - D_top.i_min, D_bottom.i_max - D_bottom.i_min);
    } else {
        init_points_long = MAX(D_left.i_max - D_left.i_min, D_right.i_max - D_right.i_min);
        init_points_short = MAX(D_top.j_max - D_top.j_min, D_bottom.j_max - D_bottom.j_min);
    }

    pointSet P_left(init_points_long), P_right(init_points_long);
    pointSet P_top(init_points_short), P_bottom(init_points_short);
    if (crack_angle != 90) { // If 90°, no need for PointSets
        funcDiscrete(G, D_left, &P_left, &crackLineLeft, crack_angle);
        funcDiscrete(G, D_right, &P_right, &crackLineRight, crack_angle);
        funcDiscrete(G, D_top, &P_top, &crackLineTop, crack_angle + 90);
        funcDiscrete(G, D_bottom, &P_bottom, &crackLineBottom, crack_angle + 90);
    }
    pointSet *P[SIDE_SIZE] = {&P_left, &P_right, &P_top, &P_bottom}; // Define List of PointSets

    // Set Source/Sink Terms
    double xr = (r / Dx);
    double yr = (r / Dy);
    for (int k = 0; k < n_sources; k++) {
        for (int i = i_source[k] - xr; i < i_source[k] + xr; i++) {
            for (int j = j_source[k] - yr; j < j_source[k] + yr; j++) {
                if (sqrt((i - i_source[k]) * (i - i_source[k]) * Dx2 + (j - j_source[k]) * (j - j_source[k]) * Dy2) <= r) {
                    b[i + j * N] = source / (M_PI * xr * yr);
                }
            }
        }
    }

    for (int k = 0; k < n_sources; k++) {
        for (int i = i_sink[k] - xr; i < i_sink[k] + xr; i++) {
            for (int j = j_sink[k] - yr; j < j_sink[k] + yr; j++) {
                if (sqrt((i - i_sink[k]) * (i - i_sink[k]) * Dx2 + (j - j_sink[k]) * (j - j_sink[k]) * Dy2) <= r) {
                    b[i + j * N] = -source / (M_PI * xr * yr);
                }
            }
        }
    }

    //// Neumann Boundary Conditions :
    // Surface 1 : (0,0) -> (L,0)
    // Surface 2 : (L,0) -> (L,M)
    // Surface 3 : (L,M) -> (0,M)
    // Surface 4 : (0,M) -> (0,0)
    // Surface 5 : Crack, centered around (x_c0,y_c0) with
    //             Length : crack_len
    //             Thickness : crack_thickness
    //             Angle : crack_angle
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            // S1 : -φ[i][2] + 4φ[i][1] - 3φ[i][0] = 0
            if (i < N - 1 && j == 0) {
                A_sp[i + N * j].set(i + 0 * N, -3);
                A_sp[i + N * j].set(i + 1 * N, +4);
                A_sp[i + N * j].set(i + 2 * N, -1);
                continue;
            }

            // S2 : φ[N-3][j] - 4φ[N-2][j] + 3φ[N-1][j] = 0
            if (j < M - 1 && i == N - 1) {
                A_sp[i + N * j].set((N - 1) + j * N, -3);
                A_sp[i + N * j].set((N - 2) + j * N, +4);
                A_sp[i + N * j].set((N - 3) + j * N, -1);
                continue;
            }

            // S3 : φ[i][N-3] - 4φ[i][N-2] + 3φ[i][N-1] = 0
            if (i > 0 && j == M - 1) {
                A_sp[i + N * j].set(i + (M - 1) * N, -3);
                A_sp[i + N * j].set(i + (M - 2) * N, +4);
                A_sp[i + N * j].set(i + (M - 3) * N, -1);
                continue;
            }

            // S4 : -φ[2][j] + 4φ[1][j] - 3φ[0][j] = 0
            if (j > 0 && i == 0) {
                A_sp[i + N * j].set(0 + j * N, -3);
                A_sp[i + N * j].set(1 + j * N, +4);
                A_sp[i + N * j].set(2 + j * N, -1);
                continue;
            }

            // S5 :
            if (crack_angle == 90) {    // Vertical Crack, don't use PointSets
                if (i == i_min_crack) { // Left Side of Crack
                    if (j > j_min_crack && j < j_max_crack) {
                        A_sp[i + N * j].set((i - 0) + j * N, -3);
                        A_sp[i + N * j].set((i - 1) + j * N, +4);
                        A_sp[i + N * j].set((i - 2) + j * N, -1);
                        continue;
                    }
                }
                if (i == i_max_crack) { // Right Side of Crack
                    if (j > j_min_crack && j < j_max_crack) {
                        A_sp[i + N * j].set((i + 0) + j * N, -3);
                        A_sp[i + N * j].set((i + 1) + j * N, +4);
                        A_sp[i + N * j].set((i + 2) + j * N, -1);
                        continue;
                    }
                }
                if (j == j_min_crack) { // Bottom Side of Crack
                    if (i > i_min_crack && i < i_max_crack) {
                        A_sp[i + N * j].set(i + N * (j - 0), -3);
                        A_sp[i + N * j].set(i + N * (j - 1), +4);
                        A_sp[i + N * j].set(i + N * (j - 2), -1);
                        continue;
                    }
                }
                if (j == j_max_crack) { // Top Side of Crack
                    if (i > i_min_crack && i < i_max_crack) {
                        A_sp[i + N * j].set(i + N * (j + 0), -3);
                        A_sp[i + N * j].set(i + N * (j + 1), +4);
                        A_sp[i + N * j].set(i + N * (j + 2), -1);
                        continue;
                    }
                }

                if (i == i_min_crack && j == j_min_crack) { // Bottom Left Corner
                    A_sp[i + N * j].set(i + N * j, 1);
                    continue;
                }
                if (i == i_max_crack && j == j_min_crack) { // Bottom Right Corner
                    A_sp[i + N * j].set(i + N * j, 1);
                    continue;
                }
                if (i == i_min_crack && j == j_max_crack) { // Top Left Corner
                    A_sp[i + N * j].set(i + N * j, 1);
                    continue;
                }
                if (i == i_max_crack && j == j_max_crack) { // Top Right Corner
                    A_sp[i + N * j].set(i + N * j, 1);
                    continue;
                }

                if ((i > i_min_crack && i < i_max_crack) &&
                    (j > j_min_crack && j < j_max_crack)) {
                    A_sp[i + N * j].set(i + N * j, 1);
                    continue;
                }
            }
            // Inner Nodes
            if (i > 0 && j > 0 && i < N - 1 && j < M - 1) {
                A_sp[i + N * j].set(i + j * N, -2 * (Dx2 + Dy2) / Dx2Dy2);
                A_sp[i + N * j].set((i + 1) + j * N, 1 / Dx2), A_sp[i + N * j].set((i - 1) + j * N, 1 / Dx2);
                A_sp[i + N * j].set(i + (j + 1) * N, 1 / Dy2), A_sp[i + N * j].set(i + (j - 1) * N, 1 / Dy2);
                continue;
            }
        }
    }

    if (crack_angle != 90) {
        for (int side = LEFT; side <= BOTTOM; side++) {
            int i = P[side]->i[0], j = P[side]->j[0];
            A_sp[i + N * j].clear();
            A_sp[i + N * j].set(i + N * j, 1);
            i = P[side]->i[P[side]->n_points - 1], j = P[side]->j[P[side]->n_points - 1];
            A_sp[i + N * j].clear();
            A_sp[i + N * j].set(i + N * j, 1);
            for (int point = 1; point < P[side]->n_points - 1; point++) {
                int i = P[side]->i[point];
                int j = P[side]->j[point];
                A_sp[i + N * j].clear();
                switch (side) {
                case LEFT:
                    A_sp[i + N * j].add((i - 0) + j * N, -3 * (sind(crack_angle)));
                    A_sp[i + N * j].add((i - 1) + j * N, +4 * (sind(crack_angle)));
                    A_sp[i + N * j].add((i - 2) + j * N, -1 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 0), -3 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 1), +4 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 2), -1 * (cosd(crack_angle)));
                    break;
                case RIGHT:
                    A_sp[i + N * j].add((i + 0) + j * N, -3 * (sind(crack_angle)));
                    A_sp[i + N * j].add((i + 1) + j * N, +4 * (sind(crack_angle)));
                    A_sp[i + N * j].add((i + 2) + j * N, -1 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 0), -3 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 1), +4 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 2), -1 * (cosd(crack_angle)));
                    break;
                case TOP:
                    A_sp[i + N * j].add((i + 0) + j * N, -3 * (cosd(crack_angle)));
                    A_sp[i + N * j].add((i + 1) + j * N, +4 * (cosd(crack_angle)));
                    A_sp[i + N * j].add((i + 2) + j * N, -1 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 0), -3 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 1), +4 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j + 2), -1 * (sind(crack_angle)));
                    break;
                case BOTTOM:
                    A_sp[i + N * j].add((i - 0) + j * N, -3 * (cosd(crack_angle)));
                    A_sp[i + N * j].add((i - 1) + j * N, +4 * (cosd(crack_angle)));
                    A_sp[i + N * j].add((i - 2) + j * N, -1 * (cosd(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 0), -3 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 1), +4 * (sind(crack_angle)));
                    A_sp[i + N * j].add(i + N * (j - 2), -1 * (sind(crack_angle)));
                    break;
                }
            }
        }
    }

    FILE *fp_A = fopen("A_Sparse.txt", "w+");
    fprintf(fp_A, "Row(int),Col(int),Data(float)\n");
    for (int i = 0; i < N * M; i++) {
        for (int j = 0; j < A_sp[i].n_nz - 1; j++) {
            if (A_sp[i].c[j] != -1) {
                fprintf(fp_A, "%d, %d, %.16lf\n", i + 1, A_sp[i].c[j] + 1, A_sp[i].v[j]);
            }
        }
    }
    fclose(fp_A);

    FILE *fp_b = fopen("b.txt", "w+");
    for (int i = 0; i < N * M; i++) {
        fprintf(fp_b, "%.16lf\n", b[i]);
    }
    fclose(fp_b);

    // Solve Sparse Linear System
    system("C:\\Miniforge3\\python.exe .\\SparseSolver.py");

    for (int i = 0; i < N * M; i++) {
        A_sp[i].free_sparse_row();
    }
    free(A_sp);
    free(b);
    free(x);
    return M != N ? M : 0;
}