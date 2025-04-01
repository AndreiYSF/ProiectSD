#include <bits/stdc++.h>
using namespace std;

//void fc() {
//
//    auto start = chrono::high_resolution_clock::now();
//
//    this_thread::sleep_for(std::chrono::milliseconds(500));
//
//    auto end = chrono::high_resolution_clock::now();
//
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//
//    cout << duration.count();
//}

void display(vector <double> v) {

    for (auto &i : v)
        cout << i << ' ';
    cout << '\n';
}

// verifica daca un vector este sortat
bool isSorted(vector<double> &v) {

    for (int i = 0; i < v.size() - 1; i++)
        if (v[i] > v[i + 1])
            return false;
    return true;
}

// RadixSort clasic
void RadixSort1(vector<double> &v) {

    vector<int> aux;

    for (auto &i: v) {

        aux.push_back(i);
        if (floor(i) != i) {

            aux.clear();
            cout << "RadixSort1 nu a putut executa sortarea.\n";
            return;
        }
    }

    for (int k = 0; k <= 30; k++) {

        bool ok = true; // verifica daca totate elementele sunt mai mici strict decat 2^k

        int pow = (1 << k);

        queue<int> leftSide, rightSide;

        for (int i = 0; i < aux.size(); i++) {

            if (aux[i] >= pow)
                ok = false;

            if (aux[i] & pow)
                rightSide.push(aux[i]);
            else
                leftSide.push(aux[i]);
        }

        aux.clear();
        while (!leftSide.empty())
            aux.push_back(leftSide.front()), leftSide.pop();
        while (!rightSide.empty())
            aux.push_back(rightSide.front()), rightSide.pop();

        if (ok)
            break;
    }

    v.clear();
    for (auto &i: aux)
        v.push_back(i);
}

// RadixSort baza 10
void RadixSort2(vector<double> &v) {

    vector<int> aux;

    for (auto &i: v) {

        aux.push_back(i);
        if (floor(i) != i) {

            aux.clear();
            cout << "RadixSort1 nu a putut executa sortarea.\n";
            return;
        }
    }

    queue<int> Q[10];

    for (int k = 0, pow = 1; k <= 9; k++, pow *= 10) {

        bool ok = true; // verifica daca totate elementele sunt mai mici strict decat 10^k

        for (int i: aux) {

            if (i >= pow)
                ok = false;

            int cif = (i / pow) % 10;

            Q[cif].push(i);
        }

        aux.clear();
        for (auto &i: Q)
            while (!i.empty())
                aux.push_back(i.front()), i.pop();

        if (ok)
            break;
    }

    v.clear();
    for (auto &i: aux)
        v.push_back(i);
}

// RadixSort baza 2^16 = 65636
void RadixSort3(vector<double> &v) {

    const int b = 65636;

    vector<int> aux;

    for (auto &i: v) {

        aux.push_back(i);
        if (floor(i) != i) {

            aux.clear();
            cout << "RadixSort1 nu a putut executa sortarea.\n";
            return;
        }
    }

    queue<int> Q[b];

    for (int pow: {1, b}) {

        for (int i: aux) {

            int cif = (i / pow) % b;

            Q[cif].push(i);
        }

        aux.clear();
        for (auto &i: Q)
            while (!i.empty())
                aux.push_back(i.front()), i.pop();
    }

    v.clear();
    for (auto &i: aux)
        v.push_back(i);
}

void MergeSort(vector<double> &v) {

    if (v.size() == 1)
        return;

    vector<double> left, right;

    int i = 0;

    while (left.size() < (v.size() / 2))
        left.push_back(v[i++]);
    while (i < v.size())
        right.push_back(v[i++]);

    MergeSort(left);
    MergeSort(right);

    i = 0;
    int j = 0, k = 0;

    while (i < left.size() && j < right.size())
        if (left[i] < right[j])
            v[k++] = left[i++];
        else
            v[k++] = right[j++];

    while (i < left.size())
        v[k++] = left[i++];
    while (j < right.size())
        v[k++] = right[j++];
}

void ShellSort(vector<double> &v) {

    int n = v.size();

    for (int gap = n / 2; gap > 0; gap /= 2) {

        for (int i = gap; i < n; i++) {

            double aux = v[i];
            int j = i;

            while (j >= gap && v[j - gap] > aux) {

                v[j] = v[j - gap];
                j -= gap;
            }
            v[j] = aux;
        }
    }
}

void QuickSortExt(vector <double>& v, int l, int r) {

    if(l >= r) return;

    int m = (l + r) / 2;

    if (v[l] > v[m])
        swap(v[l], v[m]);
    if (v[l] > v[r])
        swap(v[l], v[r]);
    if (v[m] > v[r])
        swap(v[m], v[r]);
    swap(v[m], v[r]);

    double pivot = v[r];
    int i = l;

    for(int j = l; j < r; j++) {
        if(v[j] < pivot) {
            swap(v[i], v[j]);
            i++;
        }
    }

    swap(v[i], v[r]);
    QuickSortExt(v, l, i - 1);
    QuickSortExt(v, i + 1, r);
}

void QuickSort(vector <double>& v) {

    if (v.size() <= 1)
        return;
    QuickSortExt(v, 0, v.size() - 1);
}

void heapify(vector<double>& v, int n, int i) {

    int largest = i, l = 2 * i + 1, r = 2 * i + 2;

    if(l < n && v[l] > v[largest])
        largest = l;
    if(r < n && v[r] > v[largest])
        largest = r;

    if(largest != i) {

        swap(v[i], v[largest]);
        heapify(v, n, largest);
    }
}

void HeapSort(vector<double>& v) {

    int n = v.size();
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(v, n, i);
    for (int i = n - 1; i > 0; i--) {
        swap(v[0], v[i]);
        heapify(v, i, 0);
    }
}

int main() {

    vector<double> v = {1, 8, 2, 55, 112, 12, 342, 2};

    HeapSort(v);

    cout << isSorted(v);

    return 0;
}