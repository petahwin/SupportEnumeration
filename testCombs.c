#include <stdio.h>
#include <stdlib.h>

int combinationsHelper(int * arr, int n, int k) {
    int finished = 0;
    int changed = 0;
    if (k > 0) {
        for (int i = k - 1; !finished && !changed; --i) {
            if (arr[i] < (n-1) - (k-1) + i) {
                ++arr[i];
                if (i < k - 1) {
                    for (int j = i + 1; j < k; ++j) {
                        arr[j] = arr[j-1] + 1;
                    }
                }
                changed = 1;
            }
            finished = i == 0;
        }
        if (!changed) {
            for (int i = 0; i < k; ++i) {
                arr[i] = i;
            }
        }
    }
    return changed;
}

void printVectorI(int * v, int size) {
    printf("[");
    for (int i=0; i<size; ++i) printf("%d ", v[i]);
    printf("]\n");
}

// Debug function to determine correctness of subset enumeration
void printPair(int * acc1, int * acc2, int k) {
    printf("Acc1: ");
    for (int i = 0; i < k; ++i) {
        printf("%d ", acc1[i]); 
    }
    printf("Acc2: ");
    for (int i = 0; i < k; ++i) {
        printf("%d ", acc2[i]); 
    }
    printf("\n");
}

void doNothing(int * i, int * j, int k) {
    return;
}

void combinations(int n, int k) {
    int arr1[k], arr2[k];
    for (int i = 0; i < k; ++i) arr1[i] = i, arr2[i] = i;
    do {
        do {
            printPair(arr1, arr2, k);
        } while (combinationsHelper(arr2, n, k));
    } while (combinationsHelper(arr1, n, k));
}

int main(int argc, char * argv[]) {
    if(argc != 3) {
        printf("Incorrect args: n k\n");
    } else {
        int n = atoi(argv[1]);
        // int k = atoi(argv[2]);    
        for (int i = 1; i <= n; ++i) {
            combinations(n,i);
        }
    }
}

