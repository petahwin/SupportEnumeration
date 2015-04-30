#include <stdbool.h>

extern int nActions1, nActions2;
extern int * allActions1, * allActions2;

void kSubsetsRecHelper(int k, int kCur, int * acc1, int * acc2, int index, 
                    bool proc, void (*f)(int *, int *, int)) 
{
    int * acc; 
    int * items; 
    int nActions = 0;
    if (proc) {
        acc = acc2; 
        items = allActions2; 
        nActions = nActions2;
    } else {
        acc = acc1; 
        items = allActions1; 
        nActions = nActions1;
    }

    if (kCur == 0) {
        if (proc) f(acc1, acc2, k);
        else kSubsetsRecHelper(k, k, acc1, acc2, 0, true, f);
    } else if (nActions - index < kCur || index >= nActions ) {
        return;
    } else {
        acc[k-kCur] = items[index];
        kSubsetsRecHelper(k, kCur-1, acc1, acc2, index+1, proc, f);
        kSubsetsRecHelper(k, kCur, acc1, acc2, index+1, proc, f);
    }
    
    return;
}

// Runs only for the side effects, I/O
void kSubsetsRec(int k, void (*f) (int *, int *, int)) {
    bool startProc = false;
    int acc1[k], acc2[k];
    
    for (int i = 0; i < k; ++i) {acc1[i] = 0; acc2[i] = 0;};
    kSubsetsRecHelper(k, k, acc1, acc2, 0, startProc, f);
}

int kSubsetsItHelper(int * arr, int n, int k) {
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

int kSubsetsChrItHelper(char * arr, int n, int k) {
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

void kSubsetsIt(int k, void (*f) (int *, int *, int)) {
    int arr1[k], arr2[k];
    for (int i = 0; i < k; ++i) arr1[i] = i, arr2[i] = i;
    do {
        do {
            f(arr1, arr2, k);
        } while (kSubsetsItHelper(arr2, nActions2, k));
    } while (kSubsetsItHelper(arr1, nActions1, k));
}

