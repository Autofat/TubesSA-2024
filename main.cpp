#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <vector>

using namespace std;

// BRUTEFORCE
pair<int, vector<bool>> knapsackBruteForce(const vector<int> &weights, const vector<int> &values, int capacity) {
  int n = weights.size();
  int max_value = 0;
  vector<bool> bestSelected(n, false);

  // Fungsi rekursif untuk memeriksa semua kombinasi
  function<void(int, int, int, vector<bool>)> knapsackRecursive =
      [&](int index, int current_weight, int current_value,
          vector<bool> selected) {
        if (index == n) {
          if (current_weight <= capacity && current_value > max_value) {
            max_value = current_value;
            bestSelected = selected;
          }
          return;
        }

        // Tidak memasukkan item saat ini
        knapsackRecursive(index + 1, current_weight, current_value, selected);

        // Memasukkan item saat ini
        if (current_weight + weights[index] <= capacity) {
          selected[index] = true;
          knapsackRecursive(index + 1, current_weight + weights[index],
                            current_value + values[index], selected);
          selected[index] = false;
        }
      };

  knapsackRecursive(0, 0, 0, vector<bool>(n, false));
  return make_pair(max_value, bestSelected);
}

void printSelectedItemsBruteForce(const vector<int> &weights,
                                  const vector<int> &values, int capacity) {
  auto result = knapsackBruteForce(weights, values, capacity);
  vector<bool> bestSelected = result.second;

  cout << "Kombinasi BruteForce: {";
  bool first = true;
  for (int i = 0; i < bestSelected.size(); ++i) {
    if (bestSelected[i]) {
      if (!first) {
        cout << ", ";
      }
      cout << i + 1;
      first = false;
    }
  }
  cout << "}" << endl;

  int total_weight = 0;
  for (int i = 0; i < bestSelected.size(); ++i) {
    if (bestSelected[i]) {
      total_weight += weights[i];
    }
  }
  cout << "Total Weight: " << total_weight << endl;
}

// BRANCH N BOUND
struct Node {
  int level;  // level
  int value;  // profit
  int weight; // price
  float bound;
  vector<bool> selected;

  Node(int l, int v, int w, float b, const vector<bool> &s)
      : level(l), value(v), weight(w), bound(b), selected(s) {}
};

// Fungsi untuk menghitung bound node
float calculateBound(Node u, int n, int capacity, const vector<int> &weights,
                     const vector<int> &values) {
  float profit_bound = u.value;
  int j = u.level + 1;
  int total_weight = u.weight;

  while (j < n && total_weight + weights[j] <= capacity) {
    total_weight += weights[j];
    profit_bound += values[j];
    j++;
  }

  if (j < n) {
    profit_bound += (capacity - total_weight) * values[j] / weights[j];
  }

  return profit_bound;
}

pair<int, vector<bool>> knapsackBranchAndBound(const vector<int> &weights,
                                               const vector<int> &values,
                                               int capacity) {
  int n = weights.size();
  queue<Node> Q;
  vector<bool> initial_selected(n, false);
  Node u(-1, 0, 0, 0, initial_selected);
  u.bound = calculateBound(u, n, capacity, weights, values);
  Q.push(u);

  int max_value = 0;
  vector<bool> bestSelected(n, false);

  while (!Q.empty()) {
    Node v = Q.front();
    Q.pop();

    if (v.bound > max_value) {
      u.level = v.level + 1;

      // Memasukkan item saat ini
      if (u.level < n) {
        u.weight = v.weight + weights[u.level];
        u.value = v.value + values[u.level];
        u.selected = v.selected;
        u.selected[u.level] = true;

        if (u.weight <= capacity && u.value > max_value) {
          max_value = u.value;
          bestSelected = u.selected;
        }

        u.bound = calculateBound(u, n, capacity, weights, values);

        if (u.bound > max_value) {
          Q.push(u);
        }

        // Tidak memasukkan item saat ini
        u.weight = v.weight;
        u.value = v.value;
        u.selected = v.selected;
        u.selected[u.level] = false;
        u.bound = calculateBound(u, n, capacity, weights, values);

        if (u.bound > max_value) {
          Q.push(u);
        }
      }
    }
  }

  return {max_value, bestSelected};
}

void printSelectedItemsBranchAndBound(const vector<int> &weights,
                                      const vector<int> &values, int capacity) {
  auto result = knapsackBranchAndBound(weights, values, capacity);
  vector<bool> selected = result.second;

  cout << "Kombinasi Branch and Bound: {";
  bool first = true;
  for (int i = 0; i < selected.size(); ++i) {
    if (selected[i]) {
      if (!first) {
        cout << ", ";
      }
      cout << i + 1;
      first = false;
    }
  }
  cout << "}" << endl;

  int total_weight = 0;
  for (int i = 0; i < selected.size(); ++i) {
    if (selected[i]) {
      total_weight += weights[i];
    }
  }
  cout << "Total Weight: " << total_weight << endl;
}

int main() {
  // Inisialisasi seed untuk menghasilkan angka acak
  std::srand(std::time(nullptr));
  // excel
  ofstream csvFile("execution_times.csv");

  csvFile << "N,Brute Force (ms), BnB (ms)\n";
  int i = 1;
  // int times = 0;
  while (i <= 20) {
    // Mendefinisikan vektor weights dan values
    int N = 5 * i;
    vector<int> weights(N); // Berat Barang
    vector<int> values(N);  // Nilai Keba3hagiaan
    int capacity = 150;     // Gaji Budi

    // Mengisi vektor weights dan values dengan data acak
    for (int i = 0; i < N; ++i) {
      // Menghasilkan data acak antara 20 dan 100 untuk weights
      weights[i] = std::rand() % 81 + 20;
      // Menghasilkan data acak antara 50 dan 150 untuk values
      values[i] = std::rand() % 101 + 50; // (Range: 50-150)
    }

    // Menampilkan header tabel
    printf("--------------------------------------------------------------\n");
    printf("| %-6s | %-20s | %-26s |\n", "Barang", "Berat (Harga Barang)",
           "Nilai Kebahagiaan (Profit)");
    printf("--------------------------------------------------------------\n");

    // Menampilkan isi tabel
    for (size_t i = 0; i < weights.size(); i++) {
      printf("| %-6zu | %-20d | %-26d |\n", i + 1, weights[i], values[i]);
    }

    printf("--------------------------------------------------------------\n");

    // Mengukur waktu eksekusi knapsackBruteForce
    auto start = chrono::high_resolution_clock::now();
    auto bruteForceResult = knapsackBruteForce(weights, values, capacity);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> bruteForceDuration = end - start;

    // Mengukur waktu eksekusi knapsackBranchAndBound
    start = chrono::high_resolution_clock::now();
    auto branchAndBoundResult =
        knapsackBranchAndBound(weights, values, capacity);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> branchAndBoundDuration = end - start;

    printSelectedItemsBruteForce(weights, values, capacity);
    cout << "Brute Force Profit: $" << bruteForceResult.first << endl;
    cout << "Brute Force Execution Time: " << bruteForceDuration.count()
         << " ms" << endl;
    cout << " " << endl;
    printSelectedItemsBranchAndBound(weights, values, capacity);
    cout << "Branch and Bound Profit: $" << branchAndBoundResult.first << endl;
    cout << "Branch and Bound Execution Time: "
         << branchAndBoundDuration.count() << " ms" << endl;

    // export to csv

    csvFile << N << "," << bruteForceDuration.count() << ","
            << branchAndBoundDuration.count() << "\n";
    // times++;
    // if (times == 4) {
    //   times = 0;
    //   i++;
    // }
    i++;
  }
  csvFile.close();

  return 0;
}
