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
    vector<pair<double, double>> values;
    
public:
    RandomVariable(vector<pair<double, double>> values = {}) {}
    
    double mean() const { return 0; }
    
    double variance() const { return 0; }
    
    double standardDeviation() const { return 0; }
    
    void print() {
        for(int i = 0; i< values.size(); ++i) {
            cout<< "("<< values[i].first<< ", "<< values[i].second<< ") ";
        }
        cout<< endl;
    }
};

double binomialDensity(int n, int k, double p) { return 0; }

RandomVariable binomialRandomVariable(int n, double p) { return RandomVariable(); }

double binomialRandom(int n, double p) { return 0; }

int main() {
    std::srand(std::time(nullptr));
    
    RandomVariable rv(vector<pair<double, double>>({ 
        make_pair(0, 0.5), 
        make_pair(1, 0.5) 
    }));
    rv.print();
    
    cout<< rv.mean()<< " "<< rv.variance()<< endl;
    cout<< binomialDensity(3, 2, 0.5)<< endl;
    
    RandomVariable b(binomialRandomVariable(3, 0.5));
    b.print();
    
    cout<< binomialRandom(200, 0.5);

    return 0;
}
