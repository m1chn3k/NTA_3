#include "IndexCalculus.h"
#include <chrono>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    uint64_t a = 0, b = 0, p = 0;
    try {
        if (argc == 4) {
            a = stoull(argv[1]);
            b = stoull(argv[2]);
            p = stoull(argv[3]);
        } else {
            cout << "Prime p (default 1009): ";
            string line;
            getline(cin, line);
            p = line.empty() ? 1009 : stoull(line);
            cout << "Generator a (0 = find automatically): ";
            getline(cin, line);
            a = line.empty() ? 0 : stoull(line);
            if (!a) a = find_generator(p);
            cout << "Target b (0 = create from x = 123): ";
            getline(cin, line);
            b = line.empty() ? 0 : stoull(line);
            if (!b) b = power_mod(a, 123, p);
        }
    } catch (...) { cerr << "Invalid input\n"; return 1; }
    if (p < 5 || p > 1000000000ULL || !is_prime(p)) {
        cerr << "p must be prime and <= 1000000000\n";
        return 1;
    }
    cout << "a = " << a << ", b = " << b << ", p = " << p << "\n";
    for (int mode = 0; mode < 2; ++mode) {
        uint64_t x = 0;
        int relations = 0;
        auto start = chrono::steady_clock::now();
        bool ok = index_calculus(a, b, p, mode == 1, x, relations);
        double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - start).count();
        cout << (mode ? "Parallel" : "Sequential") << ": ";
        if (ok) cout << "x = " << x << ", verified = " << (power_mod(a, x, p) == b % p);
        else cout << "no result";
        cout << ", relations = " << relations << ", time = " << ms << " ms\n";
    }
}
