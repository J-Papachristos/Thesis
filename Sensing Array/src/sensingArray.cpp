#include "libs/FunctionDiscretize.h"
#include "libs/Sparse.h"

#define X_SIZE (sizeof(double) * N * N)
#define A_SIZE (X_SIZE * N * N)
#define deg2rad(X) ((X) * (M_PI / 180.0))

#define L (2000.0 * 1e-3) // [m] Length
#define H (500.0 * 1e-3)  // [m] Height
#define t (3.0 * 1e-3)    // [m] Thickness

#define cosd(X) (cos(deg2rad(X)))
#define sind(X) (sin(deg2rad(X)))
#define tand(X) (tan(deg2rad(X)))

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
    // Mesh Data :
    int N = 1000;
    int M = (int) ((H / L) * N);
    double Dx = L / N, Dx2 = Dx * Dx;
    double Dy = H / M, Dy2 = Dy * Dy;
    double Dx2Dy2 = Dx2 * Dy2;

    // Create Grid Struct
    grid G(L, H, N, M);

    // Source/Sink Data :
    double I = 2.00;                      // [A]
    I = I / Dx / Dy;                      // [A/m^2]
    double sigma = 7476688.272623353 * t; // [siemens/m^2] * [m] = [A/Vm^2]
    double source = I / sigma;            // [A/m^2] / [A/Vm^2] = [V]

    // Data Structure Initialization
    Sparse *A_sp = (Sparse *) malloc(N * M * sizeof(Sparse));
    for (int i = 0; i < N * M; i++) {
        A_sp[i] = {NNZ_INIT, N * M};
    }

    // Crack Data :
    FILE *fp_init = fopen(argc > 1 ? argv[1] : "input.dat", "r+");
    if (!fp_init) {
        printf("Wrong File!\n");
        return -1;
    }

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
    if (crack_angle != 90 && is_cracked) {
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
    free(A_sp);

    // Fill up {b} matrix for Multi-Solve
    int N_sens = 26;
    int M_sens = 12;

    FILE *fp_b = fopen("b_Sparse.txt", "w+");
    FILE *fp_v = fopen("v_Points.txt", "w+");
    fprintf(fp_b, "Row(int),Col(int),Data(float)\n");
    fprintf(fp_v, "i_v1(int),i_v2(int),j_v1(int),j_v2(int)\n");
    // Horizontal Scan
    for (int j = 1; j <= M_sens; j++) {
        for (int i = 1; i <= N_sens - 3; i++) {
            int i_sens_Ii = (int) ((i + 0) * N / (N_sens + 1)); // Source i
            int i_sens_Io = (int) ((i + 3) * N / (N_sens + 1)); // Sink i

            int j_sens_Ii = (int) ((j + 0) * M / (M_sens + 1)); // Source j
            int j_sens_Io = (int) ((j + 0) * M / (M_sens + 1)); // Sink j

            int b_col_index = ((i - 1) + (j - 1) * (N_sens - 3)) + 1;
            fprintf(fp_b, "%d, %d, %.16lf\n", (i_sens_Ii + j_sens_Ii * N) + 1, b_col_index, source);
            fprintf(fp_b, "%d, %d, %.16lf\n", (i_sens_Io + j_sens_Io * N) + 1, b_col_index, -source);

            int i_sens_v1 = (int) ((i + 1) * N / (N_sens + 1)); // v1 i
            int i_sens_v2 = (int) ((i + 2) * N / (N_sens + 1)); // v2 i

            int j_sens_v1 = (int) ((j + 0) * M / (M_sens + 1)); // v1 j
            int j_sens_v2 = (int) ((j + 0) * M / (M_sens + 1)); // v2 j
            fprintf(fp_v, "%d, %d, %d, %d\n", i_sens_v1, i_sens_v2, j_sens_v1, j_sens_v2);
        }
    }
    // Vertical Scan
    int b_collumn_offset = (N_sens - 3) * (M_sens); // Offset the collumns to seperate Hor/Ver Scans
    for (int j = 1; j <= M_sens - 3; j++) {
        for (int i = 1; i <= N_sens; i++) {
            int i_sens_Ii = (int) ((i + 0) * N / (N_sens + 1)); // Source i
            int i_sens_Io = (int) ((i + 0) * N / (N_sens + 1)); // Sink i

            int j_sens_Ii = (int) ((j + 0) * M / (M_sens + 1)); // Source j
            int j_sens_Io = (int) ((j + 3) * M / (M_sens + 1)); // Sink j

            int b_col_index = ((i - 1) + (j - 1) * (N_sens)) + 1;
            fprintf(fp_b, "%d, %d, %.16lf\n", (i_sens_Ii + j_sens_Ii * N) + 1, b_collumn_offset + b_col_index, source);
            fprintf(fp_b, "%d, %d, %.16lf\n", (i_sens_Io + j_sens_Io * N) + 1, b_collumn_offset + b_col_index, -source);

            int i_sens_v1 = (int) ((i + 0) * N / (N_sens + 1)); // v1 i
            int i_sens_v2 = (int) ((i + 0) * N / (N_sens + 1)); // v2 i

            int j_sens_v1 = (int) ((j + 1) * M / (M_sens + 1)); // v1 j
            int j_sens_v2 = (int) ((j + 2) * M / (M_sens + 1)); // v2 j
            fprintf(fp_v, "%d, %d, %d, %d\n", i_sens_v1, i_sens_v2, j_sens_v1, j_sens_v2);
        }
    }
    fclose(fp_b);
    fclose(fp_v);

    // Output Simulation Info
    FILE *fp_info = fopen("info.txt", "w+");
    fprintf(fp_info, "%d\n%d\n", N, M);
    fprintf(fp_info, "%d\n%d\n", N_sens, M_sens);
    fclose(fp_info);

    // Solve Sparse Linear System for Multiple {b} Vectors
    system("C:\\Miniforge3\\python.exe .\\SparseMatrixSolveMultiple.py");
    return 0;
}