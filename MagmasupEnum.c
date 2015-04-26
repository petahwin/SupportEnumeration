#include "magma.h"

int main() {
    int foo[3];
    magma_sgesv(1, 1, NULL, 1, foo, NULL, 1, NULL);
    return 0;
}

