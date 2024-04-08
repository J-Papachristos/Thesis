#include <stdio.h>

#include "FunctionDiscretize.h"

double linear(double x) { // Example, remove later
    return tan(85 * 2 * M_PI / 360) * x;
}

int main(int argc, char const *argv[]) {
    grid G(2, 20);
    domain D(0, 2, 0, 2);
    D.calcIndex(G);

    pointSet P(1);
    funcDiscrete(G, D, &P, &linear, 85);
    printf("%d\n",P.n_points);
    for (int i = 0; i < P.n_points; i++) {
        printf("%lf %lf\n", P.i[i], P.j[i]);
    }

    return 0;
}
