#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double EPS = 0.0001;

bool equal(double v1, double v2) {
  return ((v1 - EPS < v2) && (v1 + EPS > v2));
}

class RandomVariable {
private:
  vector<pair<double, double>> values;
  
public:
  RandomVariable(vector<pair<double, double>> values = {}) {
    double P = 0;
    for (int i = 0; i < values.size(); ++i) {
      P += values[i].second;
    }
    if (!equal(P, 1.0)) {
      throw exception();
    }

    this->values = values;
  }
  
  double mean() const {
    double m = 0;
    for (int i = 0; i < values.size(); ++i) {
      m += values[i].first * values[i].second;
    }

    return m;
  }
  
  double variance() const {
    double v = 0;
    for (int i = 0; i < values.size(); ++i) {
      v += pow(values[i].first, 2) * values[i].second;
    }

    return v - mean();
  }
  
  double standardDeviation() const {
    return sqrt(variance());
  }
  
  void print() {
    for (int i = 0; i < values.size(); ++i) {
      cout << "(" << values[i].first << ", " << values[i].second << ") ";
    }
    cout << endl;
  }
};

int factorial(int n) {
  int f = 1;
  for (int i = n; i > 1; --i) {
    f *= i;
  }

  return f;
}

double binomialDensity(int n, int k, double p) {
  return (factorial(n) / (factorial(k) * factorial(n - k))) * pow(p, k) * pow((1 - p), (n - k));
}

RandomVariable binomialRandomVariable(int n, double p) {
  vector<pair<double, double>> values;
  for (int i = 0; i <= n; ++i) {
    values.push_back(make_pair(i, binomialDensity(n, i, p)));
  }

  return RandomVariable(values);
}

double binomialRandom(int n, double p) {
  int k = 0;
  int t = 1.0 / p;
  for (int i = 0; i < n; ++i) {
    if (!(rand() % t)) {
      ++k;
    }
  }

  return k;
}

int main() {
  srand(time(nullptr));
  
  RandomVariable rv(vector<pair<double, double>>({ 
    make_pair(0, 0.5), 
    make_pair(1, 0.5)
  }));
  rv.print();
  
  cout << rv.mean() << " " << rv.variance() << endl;
  cout << binomialDensity(3, 2, 0.5) << endl;
  
  RandomVariable b(binomialRandomVariable(3, 0.5));
  b.print();
  
  /*
  // 7.
  cout << binomialDensity(234, 40, 1.0 / 6.0) << endl;
  cout << binomialDensity(234, 40, 1.0 / 6.0) + binomialDensity(234, 39, 1.0 / 6.0) + [...] + binomialDensity(234, 1, 1.0 / 6.0) << endl;
  */

  cout << binomialRandom(200, 0.5) << endl;

  return 0;
}
