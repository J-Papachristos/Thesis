#ifndef _LIB_SPARSE_
#define _LIB_SPARSE_

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y) ((X) > (Y) ? (Y) : (X))

#define NNZ_INIT 5 // Typical # of Non Zero Elements in each row

typedef struct sparse_row {
    int max_cols; // Max Number of columns the row can have
    int n_cols;   // Initial # of Columns (Non-Zero Elements) to Allocate
    int n_nz = 1; // Number of Non Zero Elements currently stored (Index Variable)
    int *c;       // Vector of all indexes of Non Zero Elements
    double *v;    // Vector of all values of Non Zero Elements

    /// @brief Constructor
    /// @param n_cols Initial # of Columns (Non-Zero Elements) to Allocate
    sparse_row(int n_cols, int max_cols) {
        this->max_cols = max_cols;
        this->n_cols = n_cols;
        this->c = (int *) malloc(n_cols * sizeof(int));
        this->v = (double *) malloc(n_cols * sizeof(double));
        memset(this->c, -1, n_cols * sizeof(int)); // Set to -1 as indicator of NULL Value
        memset(this->v, 0, n_cols * sizeof(double));
    }

    /// @brief Destructor
    void free_sparse_row() {
        free(this->c);
        free(this->v);
    }

    /// @brief Get Value of Index j from Row.
    /// Mimics the functionality of the [] operator
    /// @param j Index of Value
    /// @return Value
    double operator[](int j) {
        if (j > this->max_cols || j < 0)
            return 0;

        for (int i = 0; i < this->n_nz; i++) {
            if (this->c[i] == j)
                return this->v[i];
        }
        return 0;
    }

    /// @brief Set Value of Index j in Row
    /// @param j Index of Value
    /// @param val Value
    /// @return 0 for Success, -1 for Error
    int set(int j, double val) {
        if (j > this->max_cols || j < 0) // If not in Bounds
            return -1;
        if (val == 0) // Only for Non-Zero Values
            return -1;

        if (!(*this)[j]) { // If pair (j,v) doesn't already exist in Row
            this->c[n_nz - 1] = j;
            this->v[n_nz - 1] = val;
            this->n_nz++; // Increment True Index

            if (n_nz > this->n_cols) { // If more memory needed
                this->c = (int *) realloc(this->c, (n_nz + 1) * sizeof(int));
                this->v = (double *) realloc(this->v, (n_nz + 1) * sizeof(double));
                this->c[n_nz - 1] = -1; // Init for get() safety
                this->v[n_nz - 1] = 0;  // Init
                this->n_cols++;         // Increment Check Value
            }
        }
        return 0;
    }

    void clear() {
        this->n_cols = 1;
        this->n_nz = 1;
        this->c = (int *) malloc(n_cols * sizeof(int));
        this->v = (double *) malloc(n_cols * sizeof(double));
        memset(this->c, -1, n_cols * sizeof(int)); // Set to -1 as indicator of NULL Value
        memset(this->v, 0, n_cols * sizeof(double));
    }
} Sparse;

#endif // !_LIB_SPARSE_