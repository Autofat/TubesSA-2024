#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <vector>

using namespace std;

struct BruteForceResult {
  int max_value;
  vector<bool> best_selected;
};

void knapsack_recursive(int index, int current_weight, int current_value,
                        const vector<int> &weights, const vector<int> &values,
                        int capacity, vector<bool> &selected,
                        BruteForceResult &result) {
  if (current_weight > capacity)
    return;
  if (index == weights.size()) {
    if (current_value > result.max_value) {
      result.max_value = current_value;
      result.best_selected = selected;
    }
    return;
  }

  // Do not include the current item
  knapsack_recursive(index + 1, current_weight, current_value, weights, values,
                     capacity, selected, result);

  // Include the current item
  selected[index] = true;
  knapsack_recursive(index + 1, current_weight + weights[index],
                     current_value + values[index], weights, values, capacity,
                     selected, result);
  selected[index] = false;
}

BruteForceResult knapsackBruteForce(const vector<int> &weights,
                                    const vector<int> &values, int capacity) {
  BruteForceResult result{0, vector<bool>(weights.size(), false)};
  vector<bool> selected(weights.size(), false);
  knapsack_recursive(0, 0, 0, weights, values, capacity, selected, result);
  return result;
}

struct Node {
  int level;
  int value;
  int weight;
  double bound;
  vector<bool> selected;

  Node(int l, int v, int w, double b, vector<bool> sel)
      : level(l), value(v), weight(w), bound(b), selected(sel) {}
};

double calculate_bound(const Node &u, int n, int capacity,
                       const vector<tuple<double, int, int>> &items) {
  if (u.weight >= capacity)
    return 0;
  double profit_bound = u.value;
  int j = u.level + 1;
  int total_weight = u.weight;

  while (j < n && total_weight + get<1>(items[j]) <= capacity) {
    total_weight += get<1>(items[j]);
    profit_bound += get<2>(items[j]);
    j++;
  }

  if (j < n) {
    profit_bound += (capacity - total_weight) * get<0>(items[j]);
  }
  return profit_bound;
}

struct Compare {
  bool operator()(const pair<double, Node> &a, const pair<double, Node> &b) {
    return a.first > b.first;
  }
};

pair<int, vector<bool>> knapsackBranchAndBound(const vector<int> &weights,
                                               const vector<int> &values,
                                               int capacity) {
  vector<tuple<double, int, int>> items;
  for (int i = 0; i < weights.size(); i++) {
    items.emplace_back((double)values[i] / weights[i], weights[i], values[i]);
  }
  sort(items.begin(), items.end(), greater<>());

  int n = items.size();
  priority_queue<pair<double, Node>, vector<pair<double, Node>>, Compare> Q;
  Node u(-1, 0, 0, 0.0, vector<bool>(n, false));
  u.bound = calculate_bound(u, n, capacity, items);
  Q.emplace(1 - u.bound, u);

  int max_value = 0;
  vector<bool> best_selected(n, false);

  while (!Q.empty()) {
    u = Q.top().second;
    Q.pop();

    if (u.bound > max_value) {
      Node v(u.level + 1, u.value, u.weight, 0.0, u.selected);
      v.weight = u.weight + get<1>(items[v.level]);
      v.value = u.value + get<2>(items[v.level]);

      if (v.weight <= capacity && v.value > max_value) {
        max_value = v.value;
        best_selected = v.selected;
      }

      v.bound = calculate_bound(v, n, capacity, items);
      if (v.bound > max_value) {
        v.selected[v.level] = true;
        Q.emplace(1 - v.bound, v);
      }

      v = Node(u.level + 1, u.value, u.weight, 0.0, u.selected);
      v.bound = calculate_bound(v, n, capacity, items);
      if (v.bound > max_value) {
        Q.emplace(1 - v.bound, v);
      }
    }
  }
  return {max_value, best_selected};
}

void printSelectedItemsBruteForce(const vector<int> &weights,
                                  const vector<int> &values, int capacity) {
  auto result = knapsackBruteForce(weights, values, capacity);
  cout << "Kombinasi BruteForce: {";
  bool first = true;
  for (int i = 0; i < result.best_selected.size(); i++) {
    if (result.best_selected[i]) {
      if (!first)
        cout << ", ";
      cout << i + 1;
      first = false;
    }
  }
  cout << "}\n";

  int total_weight = 0;
  for (int i = 0; i < result.best_selected.size(); i++) {
    if (result.best_selected[i]) {
      total_weight += weights[i];
    }
  }
  cout << "Total Weight: " << total_weight << "\n";
}

void printSelectedItemsBranchAndBound( vector<int> &weights,
                                       vector<int> &values, int capacity) {
  int n = weights.size();
  // Hitung rasio P/W dan urutkan weights dan values berdasarkan rasio tersebut
  vector<pair<double, pair<int, int>>> items(n);
  for (int i = 0; i < n; i++) {
    items[i] = {(double)values[i] / weights[i], {weights[i], values[i]}};
  }
  sort(items.begin(), items.end(), greater<>());

  for (int i = 0; i < n; i++) {
    weights[i] = items[i].second.first;
    values[i] = items[i].second.second;
  }

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
  srand(time(0));

  ofstream csv_file("execution_times.csv");
  csv_file << "N,Brute Force (ms),BnB (ms)\n";

  for (int i = 1; i <= 30; i++) {
    int N = 10 * i;
    vector<int> weights(N);
    vector<int> values(N);
    int capacity = 150;

    for (int j = 0; j < N; j++) {
      weights[j] = rand() % 81 + 20;
      values[j] = rand() % 101 + 50;
    }

    cout << "--------------------------------------------------------------\n";
    cout << "| " << left << setw(6) << "Barang"
         << " | " << setw(20) << "Berat (Harga Barang)"
         << " | " << setw(26) << "Nilai Kebahagiaan (Profit)"
         << " |\n";
    cout << "--------------------------------------------------------------\n";
    for (int j = 0; j < N; j++) {
      cout << "| " << left << setw(6) << j + 1 << " | " << setw(20)
           << weights[j] << " | " << setw(26) << values[j] << " |\n";
    }
    cout << "--------------------------------------------------------------\n";

    clock_t start_time, end_time;

    start_time = clock();
    auto brute_force_result = knapsackBruteForce(weights, values, capacity);
    end_time = clock();
    double brute_force_duration =
        1000.0 * (end_time - start_time) / CLOCKS_PER_SEC;

    start_time = clock();
    auto branch_and_bound_result =
        knapsackBranchAndBound(weights, values, capacity);
    end_time = clock();
    double branch_and_bound_duration =
        1000.0 * (end_time - start_time) / CLOCKS_PER_SEC;

    printSelectedItemsBruteForce(weights, values, capacity);
    cout << "Brute Force Profit: $" << brute_force_result.max_value << "\n";
    cout << "Brute Force Execution Time: " << fixed << setprecision(4)
         << brute_force_duration << " ms\n\n";
    printSelectedItemsBranchAndBound(weights, values, capacity);
    cout << "Branch and Bound Profit: $" << branch_and_bound_result.first
         << "\n";
    cout << "Branch and Bound Execution Time: " << fixed << setprecision(4)
         << branch_and_bound_duration << " ms\n";

    csv_file << N << "," << brute_force_duration << ","
             << branch_and_bound_duration << "\n";
  }

  csv_file.close();
  return 0;
}
