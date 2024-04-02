#include "Sparse.h"

#define X_SIZE (sizeof(double) * N * N)
#define A_SIZE (X_SIZE * N * N)
#define b2Mb(X) (((X) / 1024) / 1024)

#define L 1.0
#define H 1.0

#define N_SOURCES 1 // Number of Sources

int main(int argc, char const *argv[]) {
    /// Mesh Data :
    int N = atoi(argv[1]);
    double Dx = L / N, Dx2 = Dx * Dx;
    int M = argc > 2 ? atoi(argv[2]) : (int) (H / L) * N;
    double Dy = H / M, Dy2 = Dy * Dy;
    double Dx2Dy2 = Dx2 * Dy2;

    /// Source/Sink Data :
    double I = 0.5 / Dx / Dy;    // [A/m^2]
    double sigma = 38e6 * 0.003; // [siemens/m^2] * [m] = [A/Vm^2]
    double source = I / sigma;   // [A/m^2] / [A/Vm^2] = [V]
    double x_source = (L / 2) - 137.5e-3, y_source = (H / 2);
    double x_sink = (L / 2) + 137.5e-3, y_sink = (H / 2);
    int i_source = x_source / Dx, j_source = y_source / Dy;
    int i_sink = x_sink / Dx, j_sink = y_sink / Dy;
    // printf("(%lf,%lf) -> (%d,%d)\n", x_source, y_source, i_source, j_source);
    // printf("(%lf,%lf) -> (%d,%d)\n", x_sink, y_sink, i_sink, j_sink);

    /// Crack Data :
    // double x_c0 = L / 2, y_c0 = M / 4;
    // double x_ce = L / 2, y_ce = 3 * M / 4;
    // int i_c0 = x_c0 / Dx - 1, j_c0 = y_c0 / Dy;
    // int i_ce = x_ce / Dx + 1, j_ce = y_ce / Dy;

    double *b = (double *) malloc(X_SIZE);
    double *x = (double *) malloc(X_SIZE);
    memset(b, 0, X_SIZE);
    memset(x, 0, X_SIZE);

    Sparse *A_sp = (Sparse *) malloc(N * M * sizeof(Sparse));
    for (int i = 0; i < N * M; i++) {
        A_sp[i] = {NNZ_INIT, N * M};
    }

    // Set Source/Sink Terms
    // todo Implement Multiple Source Terms
    b[i_source + j_source * N] = source;
    b[i_sink + j_sink * N] = -source;

    //// Neumann Boundary Conditions :
    // Surface 1 : (0,0) -> (L,0)
    // Surface 2 : (L,0) -> (L,M)
    // Surface 3 : (L,M) -> (0,M)
    // Surface 4 : (0,M) -> (0,0)
    // Surface 5 : Crack, (x_c0,y_c0) -> (x_ce,y_ce)
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
    system("C:\\Users\\John\\miniforge3\\python.exe .\\src\\SparseSolver.py");

    for (int i = 0; i < N * M; i++) {
        A_sp[i].free_sparse_row();
    }
    free(A_sp);
    free(b);
    free(x);
    return 0;
}