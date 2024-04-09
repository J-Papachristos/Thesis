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
    P.printPointSet();
    P.free_point_set();
    return 0;
}
