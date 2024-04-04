#include "Sparse.h"

#define X_SIZE (sizeof(double) * N * N)
#define A_SIZE (X_SIZE * N * N)
#define b2Mb(X) (((X) / 1024) / 1024)
#define deg2rad(X) ((X) * (M_PI / 180.0))

#define L 1.0   // [m]
#define H 1.0   // [m]
#define t 0.003 // [m]

#define N_SOURCES 1 // Number of Sources

double crackLine(int i, int i0, double crack_angle) {
    return (i - i0) * tan(deg2rad(crack_angle));
}

int main(int argc, char const *argv[]) {
    /// Mesh Data :
    int N = atoi(argv[1]);
    double Dx = L / N, Dx2 = Dx * Dx;
    int M = argc > 2 ? atoi(argv[2]) : (int) (H / L) * N;
    double Dy = H / M, Dy2 = Dy * Dy;
    double Dx2Dy2 = Dx2 * Dy2;

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
    I = I / Dx / Dy;           // [A/m^2]
    double sigma = 38e6 * t;   // [siemens/m^2] * [m] = [A/Vm^2]
    double source = I / sigma; // [A/m^2] / [A/Vm^2] = [V]

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
    double x_c0, y_c0;
    double crack_len, crack_thick, crack_angle;
    int is_cracked;
    fscanf(fp_init, "is_cracked %d\n", &is_cracked);
    if (is_cracked) {
        fscanf(fp_init, "crack_pos %lf %lf\n", &x_c0, &y_c0);
        fscanf(fp_init, "crack_len %lf\n", &crack_len);
        fscanf(fp_init, "crack_thick %lf\n", &crack_thick);
        fscanf(fp_init, "crack_angle %lf\n", &crack_angle);
    } else {
        x_c0 = y_c0 = 0;
        crack_angle = crack_len = crack_thick = 0;
    }

    int i_c0 = x_c0 / Dx, j_c0 = y_c0 / Dy;
    // int i_minL = i_c0 - (int) floor((crack_len / 2.0) / sin(deg2rad(crack_angle))); // Min i based on Centerline
    // int i_maxR = i_c0 + (int) ceil((crack_len / 2.0) / sin(deg2rad(crack_angle)));  // Max i based on Centerline

    int i_cL = i_c0 - (int) floor((crack_thick / Dx / 2.0) / sin(deg2rad(crack_angle)));
    int j_cL = j_c0 + (int) ceil((crack_thick / Dy / 2.0) * cos(deg2rad(crack_angle)));
    int i_min_L = i_cL - (int) floor((crack_len / Dx / 2.0) * cos(deg2rad(crack_angle)));
    int i_max_L = i_cL + (int) ceil((crack_len / Dx / 2.0) * cos(deg2rad(crack_angle)));

    int i_cR = i_c0 + (int) ceil((crack_thick / Dx / 2.0) / sin(deg2rad(crack_angle)));
    int j_cR = j_c0 - (int) floor((crack_thick / Dy / 2.0) * cos(deg2rad(crack_angle)));
    int i_min_R = i_cR - (int) floor((crack_len / Dx / 2.0) * cos(deg2rad(crack_angle)));
    int i_max_R = i_cR + (int) ceil((crack_len / Dx / 2.0) * cos(deg2rad(crack_angle)));

    // Set Source/Sink Terms
    for (int i = 0; i < n_sources; i++) {
        b[i_source[i] + j_source[i] * N] = source;
    }
    for (int i = 0; i < n_sinks; i++) {
        b[i_sink[i] + j_sink[i] * N] = -source;
    }

    //// Neumann Boundary Conditions :
    // Surface 1 : (0,0) -> (L,0)
    // Surface 2 : (L,0) -> (L,M)
    // Surface 3 : (L,M) -> (0,M)
    // Surface 4 : (0,M) -> (0,0)
    // Surface 5 : Crack, centered around (x_c0,y_c0) with
    //              Length : crack_len
    //              Thickness : crack_thickness
    //              Angle : crack_angle
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
                A_sp[i + N * j].set(i + (N - 1) * N, -3);
                A_sp[i + N * j].set(i + (N - 2) * N, +4);
                A_sp[i + N * j].set(i + (N - 3) * N, -1);
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
            // if (i == i_c0 && j > j_c0 && j < j_ce) {
            //     A_sp[i + N * j].set((i - 0) + j * N, -3);
            //     A_sp[i + N * j].set((i - 1) + j * N, +4);
            //     A_sp[i + N * j].set((i - 2) + j * N, -1);
            //     continue;
            // }
            // if (i == i_ce && j > j_c0 && j < j_ce) {
            //     A_sp[i + N * j].set((i + 0) + j * N, -3);
            //     A_sp[i + N * j].set((i + 1) + j * N, +4);
            //     A_sp[i + N * j].set((i + 2) + j * N, -1);
            //     continue;
            // }
            // if (j == j_c0 && i >= i_c0 && i <= i_ce) {
            //     A_sp[i + N * j].set(i + N * (j - 0), -3);
            //     A_sp[i + N * j].set(i + N * (j - 1), +4);
            //     A_sp[i + N * j].set(i + N * (j - 2), -1);
            //     continue;
            // }
            // if (j == j_ce && i >= i_c0 && i <= i_ce) {
            //     A_sp[i + N * j].set(i + N * (j + 0), -3);
            //     A_sp[i + N * j].set(i + N * (j + 1), +4);
            //     A_sp[i + N * j].set(i + N * (j + 2), -1);
            //     continue;
            // }

            if (i >= i_min_L && i <= i_max_L) { // Left Side of Crack
                int jc = j_cL + crackLine(i, i_cL, crack_angle);
                int jc_prev = i - 1 >= i_min_L ? j_cL + crackLine(i - 1, i_cL, crack_angle) : jc;
                int jc_next = i + 1 <= i_max_L ? j_cL + crackLine(i + 1, i_cL, crack_angle) : jc;
                if (j >= (int) ceil((jc_prev + jc) / 2) &&
                    j <= (int) floor((jc + jc_next) / 2)) {
                    A_sp[i + N * j].set((i - 0) + j * N, -3);
                    A_sp[i + N * j].set((i - 1) + j * N, +4);
                    A_sp[i + N * j].set((i - 2) + j * N, -1);
                    continue;
                }
            }
            if (i >= i_min_R && i <= i_max_R) { // Right Side of Crack
                int jc = j_cR + crackLine(i, i_cR, crack_angle);
                int jc_prev = i - 1 >= i_min_R ? j_cR + crackLine(i - 1, i_cR, crack_angle) : jc;
                int jc_next = i + 1 <= i_max_R ? j_cR + crackLine(i + 1, i_cR, crack_angle) : jc;
                if (j >= (int) ceil((jc_prev + jc) / 2) &&
                    j <= (int) floor((jc + jc_next) / 2)) {
                    A_sp[i + N * j].set((i + 0) + j * N, -3);
                    A_sp[i + N * j].set((i + 1) + j * N, +4);
                    A_sp[i + N * j].set((i + 2) + j * N, -1);
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
    system("C:\\Users\\John\\miniforge3\\python.exe .\\SparseSolver.py");

    for (int i = 0; i < N * M; i++) {
        A_sp[i].free_sparse_row();
    }
    free(A_sp);
    free(b);
    free(x);
    return 0;
}