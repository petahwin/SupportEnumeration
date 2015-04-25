#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cuda.h>

int main() {
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("No. Devices: %d\n", devCount);
    return 0;
}

