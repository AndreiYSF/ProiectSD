#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
using namespace std;

const int MAX = 1000000 + 5;

// auxiliar memory init
double *Left = new double[MAX];
double *Right = new double[MAX];
double *aux = new double[MAX];

double *original = new double[MAX];
double *arr = new double[MAX];

void display(double *arr, int n) {
    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";
    cout << "\n";
}

bool isSorted(double *arr, int n) {
    for (int i = 1; i < n - 1; i++)
        if (arr[i] > arr[i + 1]) {
            return false;
        }
    return true;
}

void RadixSort1(double *arr, int n) {

    for (int i = 1; i <= n; i++) {


        if (floor(arr[i]) != arr[i]) {

            cout << "RadixSort1 nu a putut executa sortarea.\n";
            return;
        }
    }

    for (int k = 0; k <= 30; k++) {
        bool ok = true; // verifica daca toate elementele sunt < 2^k
        int p = (1 << k);

        int x = 0, y = 0;

        for (int i = 0; i < n; i++) {
            if (arr[i] >= p)
                ok = false;
            if (int(arr[i]) & p)
                Right[y++] = arr[i];
            else
                Left[x++] = arr[i];
        }

        for (int i = 0; i < x; i++) {
            arr[i] = Left[i];
        }
        for (int i = 0; i < y; i++) {
            arr[x + i] = Right[i];
        }

        if (ok)
            break;
    }
}

void RadixSort2(double *arr, int n) {

    for (int i = 1; i <= n; i++) {
        if (floor(arr[i]) != arr[i]) {
            cout << "RadixSort2 nu a putut executa sortarea.\n";
            return;
        }
    }

    int pow = 1;
    int auxSize = n;

    int **Q = new int *[10];
    for (int d = 0; d < 10; d++) {
        Q[d] = new int[n];
    }

    for (int k = 0; k <= 9; k++) {

        bool ok = true;
        int cnt[10] = {0};

        for (int i = 0; i < auxSize; i++) {
            if (arr[i] >= pow)
                ok = false;
            int digit = (int(arr[i]) / pow) % 10;
            Q[digit][cnt[digit]++] = arr[i];
        }

        int index = 0;
        for (int d = 0; d < 10; d++) {
            for (int j = 0; j < cnt[d]; j++) {
                arr[index++] = Q[d][j];
            }
        }

        if (ok)
            break;
        pow *= 10;
    }

    for (int d = 0; d < 10; d++) {
        delete[] Q[d];
    }
    delete[] Q;
}

// RadixSort3 (baza 2^16 = 65636)
void RadixSort3(double *arr, int n) {

    const int b = 65636;

    for (int i = 1; i <= n; i++) {

        if (floor(arr[i]) != arr[i]) {
            cout << "RadixSort3 nu a putut executa sortarea.\n";
            return;
        }
    }

    queue<double> Q[b];

    for (int pow: {1, b}) {

        for (int i = 1; i <= n; i++) {
            int cif = (int(arr[i]) / pow) % b;
            Q[cif].push(arr[i]);
        }

        int index = 0;
        for (int i = 0; i < b; i++) {
            while (!Q[i].empty()) {
                arr[++index] = Q[i].front();
                Q[i].pop();
            }
        }
    }
}

void MergeSortExt(double *arr, int l, int r) {

    if (l >= r)
        return;

    int m = (l + r) >> 1;

    MergeSortExt(arr, l, m);
    MergeSortExt(arr, m + 1, r);

    int i = l, j = m + 1, k = 0;

    while (i <= m && j <= r)
        if (arr[i] < arr[j])
            aux[++k] = arr[i++];
        else
            aux[++k] = arr[j++];

    while (i <= m)
        aux[++k] = arr[i++];
    while (j <= r)
        aux[++k] = arr[j++];

    for (i = l; i <= r; i++)
        arr[i] = aux[i - l + 1];
}

void MergeSort(double *arr, int n) {
    MergeSortExt(arr, 0, n - 1);
}

// ShellSort using arrays only
void ShellSort(double *arr, int n) {
    for (int gap = n / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i++) {
            double temp = arr[i];
            int j = i;
            while (j >= gap && arr[j - gap] > temp) {
                arr[j] = arr[j - gap];
                j -= gap;
            }
            arr[j] = temp;
        }
    }
}

// QuickSort extension (recursive) using arrays only
void QuickSortExt(double *arr, int l, int r) {

    if (l >= r) return;

    int m = (l + r) / 2;
    if (arr[l] > arr[m])
        swap(arr[l], arr[m]);
    if (arr[l] > arr[r])
        swap(arr[l], arr[r]);
    if (arr[m] > arr[r])
        swap(arr[m], arr[r]);
    swap(arr[m], arr[r]);

    double pivot = arr[r];
    int i = l;
    for (int j = l; j < r; j++) {
        if (arr[j] < pivot) {
            swap(arr[i], arr[j]);
            i++;
        }
    }
    swap(arr[i], arr[r]);
    QuickSortExt(arr, l, i - 1);
    QuickSortExt(arr, i + 1, r);
}

void QuickSort(double *arr, int n) {
    QuickSortExt(arr, 0, n - 1);
}

void heapify(double *arr, int n, int i) {

    int largest = i, l = 2 * i + 1, r = 2 * i + 2;
    if (l < n && arr[l] > arr[largest])
        largest = l;
    if (r < n && arr[r] > arr[largest])
        largest = r;
    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void HeapSort(double *arr, int n) {
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
    for (int i = n - 1; i > 0; i--) {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

class Set {
private:
    struct node {
        double key;
        int frg;
        node *left, *right;
        int height;

        node(double k) : key(k), frg(1), left(NULL), right(NULL), height(1) {}
    };

    node *root;

    int height(node *x) {
        return x == NULL ? 0 : x->height;
    }

    int getBalanced(node *x) {
        return x == NULL ? 0 : height(x->left) - height(x->right);
    }

    node *toRight(node *x) {
        node *y = x->left;
        node *z = y->right;
        y->right = x;
        x->left = z;
        x->height = max(height(x->left), height(x->right)) + 1;
        y->height = max(height(y->left), height(y->right)) + 1;
        return y;
    }

    node *toLeft(node *x) {
        node *y = x->right;
        node *z = y->left;
        y->left = x;
        x->right = z;
        x->height = max(height(x->left), height(x->right)) + 1;
        y->height = max(height(y->left), height(y->right)) + 1;
        return y;
    }

    node *addNode(node *x, double key) {
        if (x == NULL)
            return new node(key);
        if (key < x->key)
            x->left = addNode(x->left, key);
        else if (key > x->key)
            x->right = addNode(x->right, key);
        else {
            x->frg++;
            return x;
        }

        x->height = max(height(x->left), height(x->right)) + 1;
        int balance = getBalanced(x);

        if (balance > 1) {
            if (getBalanced(x->left) >= 0)
                return toRight(x);
            else {
                x->left = toLeft(x->left);
                return toRight(x);
            }
        }
        else if (balance < -1) {
            if (getBalanced(x->right) <= 0)
                return toLeft(x);
            else {
                x->right = toRight(x->right);
                return toLeft(x);
            }
        }
        return x;
    }

    void clean(node *x) {
        if (x == NULL)
            return;
        clean(x->left);
        clean(x->right);
        delete x;
    }

    void rec(node *x, double *arr, int &idx) {
        if (x == NULL)
            return;
        rec(x->left, arr, idx);
        for (int i = 0; i < x->frg; i++) {
            arr[++idx] = x->key;
        }
        rec(x->right, arr, idx);
    }

public:
    Set() { root = NULL; }

    ~Set() { clean(root); }

    void add(double k) { root = addNode(root, k); }

    void reset() {
        clean(root);
        root = NULL;
    }

    void generateArray(double *arr) {

        int idx = 0;
        rec(root, arr, idx);
    }
};

void AVLSort(double *arr, int n) {

    Set S;
    for (int i = 0; i < n; i++) {
        S.add(arr[i]);
    }

    S.generateArray(arr);
}

void stdSort(double *arr, int n) {
    std::sort(arr, arr + n);
}

typedef void (*SortFunction)(double *, int);

struct Sort {
    const char *name;
    SortFunction f;
};

int main() {

    ifstream fin("number.in");

    int T;
    cin >> T;

    Sort sorts[] = {
            {"RadixSort1", RadixSort1},
            {"RadixSort2", RadixSort2},
            {"RadixSort3", RadixSort3},
            {"MergeSort",  MergeSort},
            {"ShellSort",  ShellSort},
            {"QuickSort",  QuickSort},
            {"HeapSort",   HeapSort},
            {"AVLSort",    AVLSort},
            {"std::sort",  stdSort}
    };
    int size = 9;

    for (int t = 0; t < T; t++) {
        int N, Max;
        char type; // i for integers / d for doubles
        cin >> N >> Max >> type;

        if (type == 'i') {

            mt19937 rng(random_device{}());
            uniform_int_distribution<int> dist(0, Max);

            for (int i = 1; i <= N; i++) {
                original[i] = dist(rng);
            }

        }
        else if (type == 'd') {

            mt19937 rng(random_device{}());
            uniform_real_distribution<double> dist(0, Max);

            for (int i = 1; i <= N; i++) {
                original[i] = dist(rng);
            }
        }

        // Run each sorting algorithm

        for (int i = 0; i < size; i++) {

            for (int j = 1; j <= N; j++) {

                arr[j] = original[j];
            }
            auto start = chrono::high_resolution_clock::now();
            sorts[i].f(arr, N);
            auto end = chrono::high_resolution_clock::now();

            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            bool sortedFlag = isSorted(arr, N);
            cout << sorts[i].name << ": " << duration << " ms, sorted: " << (sortedFlag ? "yes" : "no") << "\n";
        }
        cout << "\n";
    }
    return 0;
}