#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define A 420

static double basic_uniform_distributed_value(unsigned int *seedp); // [0; 1]
static double* generate_m1(int N, unsigned int *seedp);
static double* generate_m2(int N, unsigned int *seedp);

static void map_m1(double* m1, int N);
static void map_m2(double* m2, int N);

static void merge(const double* m1, double* m2, int N);

static void heapSort(double* arr, int N);

static double reduce(double* m2, int N);


int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "You should provide exactly one command-line arg: N\n");
        return 1;
    }

    struct timeval T1, T2;

    int N = atoi(argv[1]);

    if (N < 2) {
        fprintf(stderr, "N should be equal or more than 2\n");
        return 1;
    }

    gettimeofday(&T1, NULL);

    double x;
    for (int i = 0; i < 100; i++) {
        srand(i);

        // ----------- STEP GENERATE --------------
        double* m1 = generate_m1(N, (unsigned int*) &i);
        double* m2 = generate_m2(N, (unsigned int*) &i);
        // ----------------------------------------

        // ----------- STEP MAP -------------------
        map_m1(m1, N);
        map_m2(m2, N);
        // ----------------------------------------

        // ----------- STEP MERGE------------------
        merge(m1, m2, N);
        // ----------------------------------------

        // ----------- STEP SORT ------------------
        // N / 2, because heapSort() (algo which has been found on the internet) expects second arg as
        // array size
        heapSort(m2, N / 2);
        // ----------------------------------------

        // ----------- STEP REDUCE-----------------
        x = reduce(m2, N);
        // ----------------------------------------

        // ----------- STEP FREE ------------------
        free(m1);
        free(m2);
        // ----------------------------------------

//        printf("\n\nX: %lf\n\n", x);
    }

    gettimeofday(&T2, NULL);
    long delta_ms = (T2.tv_sec - T1.tv_sec) * 1000 + (T2.tv_usec - T1.tv_usec) / 1000;
    printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

    x += 0.0; // just trick to prevent CLion of highlighting x = reduce(m2, N);

    return 0;
}

// ------------------------------ IMPLS --------------------------------
static double basic_uniform_distributed_value(unsigned int *seedp) {
    return ((double) rand_r(seedp)) / RAND_MAX;
}

static double* generate_m1(int N, unsigned int *seedp) {
    double* m1 = (double*) malloc(sizeof(double) * N);

    for (int i = 0; i < N; i++) {
        m1[i] = basic_uniform_distributed_value(seedp) * (A - 1) + 1;
    }

    return m1;
}

static double* generate_m2(int N, unsigned int *seedp) {
    int l = N / 2;
    double* m2 = (double*) malloc(sizeof(double) * l);

    for (int i = 0; i < l; i++) {
        m2[i] = basic_uniform_distributed_value(seedp) * (10 * A - A) + A;
    }

    return m2;
}

static void map_m1(double* m1, int N) {
    for (int i = 0; i < N; i++) {
        m1[i] = tanh(m1[i]) - 1;
    }
}

static void map_m2(double* m2, int N) {
    double* copy = (double*) malloc(sizeof(double) * (N / 2));
    for (int i = 0; i < N / 2; i++) {
        copy[i] = m2[i];
    }

    int prev_i = -1;
    for (int i = 0; i < N / 2; i++) {
        // no-op for first element (sum with zero)
        if (prev_i != -1) {
            m2[i] += copy[prev_i];
        }

        prev_i = i;
    }

    double tmp;
    for (int i = 0; i < N / 2; i++) {
        tmp = tan(m2[i]);
        if (tmp < 0) {
            tmp = -tmp;
        }

        m2[i] = log(tmp);
    }

    free(copy);
}

static void merge(const double* m1, double* m2, int N) {
    // only i until N / 2, because size(m2) == N / 2, size(m1) == N
    for (int i = 0; i < N / 2; i++) {
        m2[i] = m1[i] * m2[i];
    }
}

// Function to swap the position of two elements
static void swap(double* a, double* b) {

    double temp = *a;
    *a = *b;
    *b = temp;
}

// To heapify a subtree rooted with node i
// which is an index in arr[].
// n is size of heap
static void heapify(double* arr, int N, int i) {
    // Find largest among root,
    // left child and right child

    // Initialize largest as root
    int largest, left, right;

    // CHANGED INITIAL ALGO: Recursion => cycle
    while (1) {
        largest = i;
        left = 2 * i + 1;
        right = 2 * i + 2;

        // If left child is larger than root
        if (left < N && arr[left] > arr[largest]) {
            largest = left;
        }

        // If right child is larger than largest
        // so far
        if (right < N && arr[right] > arr[largest]) {
            largest = right;
        }

        // Swap and continue heapifying
        // if root is not largest
        // If largest is not root
        if (largest != i) {
            swap(&arr[i], &arr[largest]);

            // Heapify the affected sub-tree
            i = largest;
        } else {
            break;
        }
    }
}

// Main function to do heap sort
static void heapSort(double* arr, int N) {
    // Build max heap
    for (int i = N / 2 - 1; i >= 0; i--) {
        heapify(arr, N, i);
    }

    // Heap sort
    for (int i = N - 1; i >= 0; i--) {
        swap(&arr[0], &arr[i]);
        // Heapify root element
        // to get highest element at
        // root again
        heapify(arr, i, 0);
    }
}

static double reduce(double* m2, int N) {
    int first_not_zero_i = -1;
    for (int i = 0; i < N / 2; i++) {
        if (m2[i] > 0) {
            first_not_zero_i = i;
            break;
        }
    }

    if (first_not_zero_i == -1) {
        fprintf(stderr, "STEP REDUCE: no not-zero elements are found in m2 array\n");
        return 0.0;
    }

    double min_not_zero = m2[first_not_zero_i];
    for (int i = first_not_zero_i + 1; i < N / 2; i++) {
        if (m2[i] < min_not_zero) {
            min_not_zero = m2[i];
        }
    }

    double x = 0.0;
    double division;
    for (int i = 0; i < N / 2; i++) {
        division = m2[i] / min_not_zero;
        if (((int) division) % 2 == 0) {
            x += sin(m2[i]);
        }
    }

    return x;
}
// ---------------------------------------------------------------------