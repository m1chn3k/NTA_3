#pragma once
#include <cstdint>

bool index_calculus(uint64_t a, uint64_t b, uint64_t p, bool parallel, uint64_t& answer, int& relations);
uint64_t find_generator(uint64_t p);
uint64_t power_mod(uint64_t a, uint64_t x, uint64_t p);
bool is_prime(uint64_t n);
