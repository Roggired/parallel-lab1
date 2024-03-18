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

#define EXPERIMENTS_NUMBER 5
#define THREADS_NUMBER 8
#define SCHEDULE schedule(guided, 1)

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

    double** m11 = (double**) malloc(sizeof(double*) * EXPERIMENTS_NUMBER);
    double** m12 = (double**) malloc(sizeof(double*) * EXPERIMENTS_NUMBER);

    int k;
    for (int i = 0; i < EXPERIMENTS_NUMBER; i++) {
        k = i;

        m11[i] = generate_m1(N, (unsigned int*) &k);
        m12[i] = generate_m2(N, (unsigned int*) &k);
    }

    gettimeofday(&T1, NULL);

    for (int i = 0; i < EXPERIMENTS_NUMBER; i++) {
        double* m1 = malloc(sizeof(double) * N);
        #pragma omp parallel for default(none) shared(m1, m11, N, i) num_threads(THREADS_NUMBER) SCHEDULE
        for (int j = 0; j < N; j++) {
            m1[j] = m11[i][j];
        }

        double* m2 = malloc(sizeof(double) * N / 2);
        #pragma omp parallel for default(none) shared(m2, m12, N, i) num_threads(THREADS_NUMBER) SCHEDULE
        for (int j = 0; j < N / 2; j++) {
            m2[j] = m12[i][j];
        }

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
        double x = reduce(m2, N);
        // ----------------------------------------

        // ----------- STEP FREE ------------------
        free(m1);
        free(m2);
        // ----------------------------------------
        printf("\n\nX: %lf\n\n", x);
        x = x + 0.0;
    }

    gettimeofday(&T2, NULL);
    long delta_ms = (T2.tv_sec - T1.tv_sec) * 1000 + (T2.tv_usec - T1.tv_usec) / 1000;
    printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

    for (int i = 0; i < EXPERIMENTS_NUMBER; i++) {
        free(m11[i]);
        free(m12[i]);
    }
    free(m11);
    free(m12);

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
    #pragma omp parallel for default(none) shared(m1, N) num_threads(THREADS_NUMBER) SCHEDULE
    for (int i = 0; i < N; i++) {
        m1[i] = tanh(m1[i]) - 1;
    }
}

static void map_m2(double* m2, int N) {
    double* copy = (double*) malloc(sizeof(double) * (N / 2));
    #pragma omp parallel for default(none) shared(copy, m2, N) num_threads(THREADS_NUMBER) SCHEDULE
    for (int i = 0; i < N / 2; i++) {
        copy[i] = m2[i];
    }

    int prev_i = -1;
    #pragma omp parallel for default(none) firstprivate(prev_i) shared(copy, m2, N) num_threads(THREADS_NUMBER) SCHEDULE
    for (int i = 0; i < N / 2; i++) {
        // no-op for first element (sum with zero)
        if (prev_i != -1) {
            m2[i] += copy[prev_i];
        }

        prev_i = i;
    }

    double tmp;
    #pragma omp parallel for default(none) private(tmp) shared(m2, N) num_threads(THREADS_NUMBER) SCHEDULE
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
    #pragma omp parallel for default(none) shared(m1, m2, N) num_threads(THREADS_NUMBER) SCHEDULE
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
    int a_start = 0, a_end = N / 2;
    int b_start = a_end, b_end = N;

    int a_size = N / 2;
    int b_size;
    if (N % 2 == 0) {
        b_size = a_size;
    } else {
        b_size = a_size + 1;
    }

    double* arr_a = (double*) malloc(sizeof(double*) * a_size);
    double* arr_b = (double*) malloc(sizeof(double*) * b_size);
    #pragma omp parallel sections shared(arr, N, a_start, a_end, b_start, b_end, arr_a, arr_b, a_size, b_size)
    {
        #pragma omp section
        {
            for (int i = a_start; i < a_end; i++) {
                arr_a[i] = arr[i];
            }

            // Build max heap
            for (int i = a_size / 2 - 1; i >= 0; i--) {
                heapify(arr_a, a_size, i);
            }

            // Heap sort
            for (int i = a_size - 1; i >= 0; i--) {
                swap(&arr_a[0], &arr_a[i]);
                // Heapify root element
                // to get highest element at
                // root again
                heapify(arr_a, i, 0);
            }
        }

        #pragma omp section
        {
            for (int i = b_start; i < b_end; i++) {
                arr_b[i - b_start] = arr[i];
            }

            // Build max heap
            for (int i = b_size / 2 - 1; i >= 0; i--) {
                heapify(arr_b, b_size, i);
            }

            // Heap sort
            for (int i = b_size - 1; i >= 0; i--) {
                swap(&arr_b[0], &arr_b[i]);
                // Heapify root element
                // to get highest element at
                // root again
                heapify(arr_b, i, 0);
            }
        }
    }

    int a_i = 0;
    int b_i = 0;
    for (int i = 0; i < N; i++) {
        if (a_i >= a_size) {
            arr[i] = arr_b[b_i];
            b_i++;
            continue;
        }

        if (b_i >= b_size) {
            arr[i] = arr_a[a_i];
            a_i++;
            continue;
        }

        if (arr_a[a_i] <= arr_b[b_i]) {
            arr[i] = arr_a[a_i];
            a_i++;
        } else {
            arr[i] = arr_b[b_i];
            b_i++;
        }
    }

    free(arr_a);
    free(arr_b);
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
    #pragma omp parallel for default(none) shared(m2, min_not_zero, N) reduction(+:x) num_threads(THREADS_NUMBER) SCHEDULE
    for (int i = 0; i < N / 2; i++) {
        double division = m2[i] / min_not_zero;
        if (((int) division) % 2 == 0) {
            x += sin(m2[i]);
        }
    }

    return x;
}
// ---------------------------------------------------------------------