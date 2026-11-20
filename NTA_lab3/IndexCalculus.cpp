#include "IndexCalculus.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

using namespace std;
using Row = vector<int64_t>;

static uint64_t mul(uint64_t a, uint64_t b, uint64_t p) {
    return a * b % p;
}

uint64_t power_mod(uint64_t a, uint64_t x, uint64_t p) {
    uint64_t result = 1;
    for (; x; x >>= 1) {
        if (x & 1) result = mul(result, a, p);
        a = mul(a, a, p);
    }
    return result;
}

bool is_prime(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t d = 2; d <= n / d; ++d)
        if (n % d == 0) return false;
    return true;
}

static vector<pair<uint64_t, int>> factor(uint64_t n) {
    vector<pair<uint64_t, int>> result;
    for (uint64_t d = 2; d <= n / d; ++d) {
        if (n % d) continue;
        int count = 0;
        do { n /= d; ++count; } while (n % d == 0);
        result.push_back({d, count});
    }
    if (n > 1) result.push_back({n, 1});
    return result;
}

uint64_t find_generator(uint64_t p) {
    auto factors = factor(p - 1);
    for (uint64_t a = 2; a < p; ++a) {
        bool good = true;
        for (auto [q, count] : factors)
            if (power_mod(a, (p - 1) / q, p) == 1) { good = false; break; }
        if (good) return a;
    }
    return 0;
}

static int64_t mod(int64_t a, int64_t p) {
    a %= p;
    return a < 0 ? a + p : a;
}

static int64_t inverse(int64_t a, int64_t m) {
    int64_t b = m, x = 1, y = 0;
    while (b) {
        int64_t q = a / b;
        int64_t t = a % b; a = b; b = t;
        t = x - q * y; x = y; y = t;
    }
    return a == 1 ? mod(x, m) : -1;
}

static vector<int> prime_base(uint64_t p) {
    double logn = log((double)(p - 1));
    int bound = (int)(3.38 * exp(0.5 * sqrt(logn * log(logn))));
    bound = max(11, min(bound, 500));
    vector<int> primes;
    for (int i = 2; i < bound && (uint64_t)i < p; ++i) {
        bool prime = true;
        for (int j = 2; j * j <= i; ++j)
            if (i % j == 0) { prime = false; break; }
        if (prime) primes.push_back(i);
    }
    return primes;
}

static bool smooth(uint64_t x, const vector<int>& base, Row& exponents) {
    exponents.assign(base.size(), 0);
    for (size_t i = 0; i < base.size(); ++i) {
        while (x % base[i] == 0) {
            x /= base[i];
            ++exponents[i];
        }
        if (x == 1) break;
    }
    return x == 1;
}

static bool solve_mod_prime(const vector<Row>& rows, int columns, int64_t q,
                            const vector<int64_t>& right, vector<int64_t>& answer) {
    int n = (int)rows.size();
    vector<Row> a(n, Row(columns + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < columns; ++j) a[i][j] = mod(rows[i][j], q);
        a[i][columns] = mod(right[i], q);
    }
    vector<int> pivot(columns, -1);
    int rank = 0;
    for (int col = 0; col < columns; ++col) {
        int best = rank;
        while (best < n && !a[best][col]) ++best;
        if (best == n) return false;
        swap(a[best], a[rank]);
        int64_t inv = inverse(a[rank][col], q);
        if (inv < 0) return false;
        for (int j = col; j <= columns; ++j) a[rank][j] = a[rank][j] * inv % q;
        for (int i = 0; i < n; ++i) {
            if (i == rank || !a[i][col]) continue;
            int64_t k = a[i][col];
            for (int j = col; j <= columns; ++j)
                a[i][j] = mod(a[i][j] - k * a[rank][j], q);
        }
        pivot[col] = rank++;
    }
    for (int i = rank; i < n; ++i)
        if (a[i][columns]) return false;
    answer.resize(columns);
    for (int j = 0; j < columns; ++j) answer[j] = a[pivot[j]][columns];
    return true;
}

static bool get_logs(const vector<Row>& rows, int columns, uint64_t groupOrder,
                     vector<int64_t>& logs) {
    auto factors = factor(groupOrder);
    vector<int64_t> residues(columns, 0);
    uint64_t currentModulus = 1;
    for (auto [q, exponent] : factors) {
        vector<int64_t> right(rows.size());
        for (size_t i = 0; i < rows.size(); ++i) right[i] = rows[i].back();
        vector<int64_t> x;
        if (!solve_mod_prime(rows, columns, (int64_t)q, right, x)) return false;
        uint64_t power = q;
        for (int digit = 1; digit < exponent; ++digit) {
            vector<int64_t> corrections(rows.size());
            uint64_t next = power * q;
            for (size_t i = 0; i < rows.size(); ++i) {
                int64_t sum = 0;
                for (int j = 0; j < columns; ++j)
                    sum = (sum + rows[i][j] * x[j]) % (int64_t)next;
                int64_t difference = mod(rows[i].back() - sum, (int64_t)next);
                if (difference % (int64_t)power) return false;
                corrections[i] = (difference / (int64_t)power) % (int64_t)q;
            }
            vector<int64_t> delta;
            if (!solve_mod_prime(rows, columns, (int64_t)q, corrections, delta)) return false;
            for (int j = 0; j < columns; ++j) x[j] += (int64_t)power * delta[j];
            power = next;
        }
        int64_t inv = inverse((int64_t)(currentModulus % power), (int64_t)power);
        if (inv < 0) return false;
        for (int j = 0; j < columns; ++j) {
            int64_t diff = mod(x[j] - residues[j], (int64_t)power);
            int64_t t = diff * inv % (int64_t)power;
            residues[j] += (int64_t)currentModulus * t;
        }
        currentModulus *= power;
    }
    logs = move(residues);
    return true;
}

bool index_calculus(uint64_t a, uint64_t b, uint64_t p, bool parallel,
                    uint64_t& answer, int& relations) {
    if (p < 5 || p > 1000000000ULL || !is_prime(p)) return false;
    a %= p;
    b %= p;
    if (!a || !b) return false;
    uint64_t n = p - 1;
    for (auto [q, exponent] : factor(n))
        if (power_mod(a, n / q, p) == 1) return false;

    auto base = prime_base(p);
    int columns = (int)base.size();
    vector<Row> rows;
    rows.reserve(columns + 160);
    uint64_t maxAttempts = 350000;
    uint64_t offset = 1;
    uint64_t step = n > 2 ? n - 1 : 1;
    atomic<uint64_t> next(0);
    atomic<int> accepted(0);
    int target = columns + 130;
    int threadCount = parallel ? (int)thread::hardware_concurrency() : 1;
    threadCount = max(1, min(threadCount, 8));
    vector<vector<Row>> found(threadCount);
    auto worker = [&](int id) {
        while (true) {
            uint64_t i = next.fetch_add(1);
            if (i >= maxAttempts || accepted.load() >= target) break;
            uint64_t k = (offset + i * step) % n;
            Row powers;
            if (smooth(power_mod(a, k, p), base, powers)) {
                powers.push_back((int64_t)k);
                if (accepted.fetch_add(1) < target) found[id].push_back(move(powers));
            }
        }
    };
    vector<thread> threads;
    for (int i = 1; i < threadCount; ++i) threads.emplace_back(worker, i);
    worker(0);
    for (auto& t : threads) t.join();
    for (auto& part : found) {
        for (auto& row : part) rows.push_back(move(row));
    }
    relations = (int)rows.size();
    if ((int)rows.size() < columns) return false;
    vector<int64_t> logs;
    if (!get_logs(rows, columns, n, logs)) return false;
    for (int j = 0; j < columns; ++j)
        if (power_mod(a, logs[j], p) != (uint64_t)base[j]) return false;

    uint64_t current = b;
    for (uint64_t l = 0; l < 1000000; ++l) {
        Row powers;
        if (smooth(current, base, powers)) {
            int64_t x = -(int64_t)(l % n);
            for (int j = 0; j < columns; ++j)
                x = (x + powers[j] * logs[j]) % (int64_t)n;
            answer = (uint64_t)mod(x, (int64_t)n);
            if (power_mod(a, answer, p) == b) return true;
        }
        current = mul(current, a, p);
    }
    return false;
}
